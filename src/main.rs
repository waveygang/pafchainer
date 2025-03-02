use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::error::Error;
use clap::Parser;
use log::{info, warn};
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus};
use rust_htslib::faidx::Reader as FastaReader;

#[derive(Parser, Debug)]
#[clap(author, version, about = "Process PAF files and connect chains using WFA2 alignment")]
struct Args {
    /// Input PAF file
    #[clap(short, long)]
    paf: PathBuf,
    
    /// Query FASTA file
    #[clap(short, long)]
    query: PathBuf,
    
    /// Target FASTA file
    #[clap(short, long)]
    target: PathBuf,
    
    /// Size of boundary erosion in base pairs
    #[clap(short, long, default_value = "1000")]
    erosion_size: usize,
    
    /// Output PAF file
    #[clap(short, long)]
    output: Option<PathBuf>,
    
    /// Number of threads to use
    #[clap(short, long, default_value = "4")]
    threads: usize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[arg(short, long, default_value = "0")]
    verbose: u8,
}

#[derive(Debug, Clone)]
struct PafEntry {
    query_name: String,
    query_length: u64,
    query_start: u64,
    query_end: u64,
    strand: char,
    target_name: String,
    target_length: u64,
    target_start: u64,
    target_end: u64,
    num_matches: u64,
    alignment_length: u64,
    mapping_quality: u64,
    tags: Vec<(String, String)>, // Use Vec to preserve tag order
    chain_id: u64,
    chain_length: u64,
    chain_pos: u64,
    cigar: String,
}

impl PafEntry {
    fn parse_from_line(line: &str) -> Result<Self, Box<dyn Error>> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err(format!("PAF line has fewer than 12 required fields: {}", line).into());
        }

        let mut entry = PafEntry {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse()?,
            query_start: fields[2].parse()?,
            query_end: fields[3].parse()?,
            strand: if fields[4] == "+" { '+' } else { '-' },
            target_name: fields[5].to_string(),
            target_length: fields[6].parse()?,
            target_start: fields[7].parse()?,
            target_end: fields[8].parse()?,
            num_matches: fields[9].parse()?,
            alignment_length: fields[10].parse()?,
            mapping_quality: fields[11].parse()?,
            tags: Vec::new(),
            chain_id: 0,
            chain_length: 0,
            chain_pos: 0,
            cigar: String::new(),
        };

        let mut has_chain_tag = false;
        let mut has_cigar_tag = false;

        // Parse optional tags
        for i in 12..fields.len() {
            let tag_parts: Vec<&str> = fields[i].split(':').collect();
            if tag_parts.len() >= 3 {
                let tag_name = format!("{}:{}", tag_parts[0], tag_parts[1]);
                let tag_value = tag_parts[2..].join(":");
                
                // Parse chain information
                if tag_name == "chain:i" {
                    has_chain_tag = true;
                    let chain_parts: Vec<&str> = tag_value.split('.').collect();
                    if chain_parts.len() >= 3 {
                        entry.chain_id = chain_parts[0].parse()?;
                        entry.chain_length = chain_parts[1].parse()?;
                        entry.chain_pos = chain_parts[2].parse()?;
                    } else {
                        panic!("Invalid chain:i tag format: {}", tag_value);
                    }
                }
                
                // Parse CIGAR string
                if tag_name == "cg:Z" {
                    has_cigar_tag = true;
                    entry.cigar = tag_value.clone();
                }
                
                entry.tags.push((tag_name, tag_value));
            }
        }

        // Panic if required tags are missing
        if !has_chain_tag {
            panic!("Missing required 'chain:i' tag in PAF line: {}", line);
        }
        if !has_cigar_tag {
            panic!("Missing required 'cg:Z' tag in PAF line: {}", line);
        }

        Ok(entry)
    }
    
    fn to_string(&self) -> String {
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            if self.strand == '+' { "+" } else { "-" },
            self.target_name,
            self.target_length,
            self.target_start,
            self.target_end,
            self.num_matches,
            self.alignment_length,
            self.mapping_quality
        );
        
        // Add tags in the expected order
        for (tag, value) in &self.tags {
            line.push_str(&format!("\t{}:{}", tag, value));
        }
        
        line
    }
}

/// Sequence database for holding FASTA sequences and readers
struct SequenceDB {
    reader: FastaReader,
}

impl SequenceDB {
    fn new<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        // Create a FASTA reader with index
        let reader = FastaReader::from_path(path)?;
        
        Ok(Self {
            reader,
        })
    }
    
    fn get_subsequence(&self, name: &str, start: u64, end: u64, reverse_complement: bool) -> Result<Vec<u8>, Box<dyn Error>> {
        // Convert to usize for rust-htslib
        let start_usize = start as usize;
        let end_usize = std::cmp::min(end as usize, std::usize::MAX);
        
        // Fetch sequence - rust-htslib is zero-based for start and inclusive for end
        let seq = match self.reader.fetch_seq_string(name, start_usize, end_usize - 1) {
            Ok(seq) => seq.into_bytes(),
            Err(e) => return Err(format!("Failed to fetch sequence for {}: {}", name, e).into()),
        };
        
        if reverse_complement {
            // Implement reverse complement
            Ok(seq.iter().rev().map(|&b| match b {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'G' | b'g' => b'C',
                b'C' | b'c' => b'G',
                x => x,
            }).collect())
        } else {
            Ok(seq)
        }
    }
}

/// Represents a CIGAR operation and its count
#[derive(Debug, Clone, Copy)]
struct CigarOp(char, usize);

/// Parse a CIGAR string into a vector of operations
fn parse_cigar(cigar: &str) -> Result<Vec<CigarOp>, Box<dyn Error>> {
    let mut ops = Vec::new();
    let mut count = 0;
    
    for c in cigar.chars() {
        if c.is_digit(10) {
            count = count * 10 + c.to_digit(10).unwrap() as usize;
        } else {
            ops.push(CigarOp(c, count));
            count = 0;
        }
    }
    
    if count > 0 {
        return Err("CIGAR string ended with a number".into());
    }
    
    Ok(ops)
}

/// Calculate statistics from a CIGAR string
/// Returns (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len)
fn calculate_cigar_stats(cigar: &str) -> Result<(u64, u64, u64, u64, u64, u64, u64), Box<dyn Error>> {
    let ops = parse_cigar(cigar)?;
    let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) = 
        ops.iter().fold((0, 0, 0, 0, 0, 0, 0), |(m, mm, i, i_bp, d, d_bp, bl), &CigarOp(op, len)| {
            let len = len as u64;
            match op {
                'M' => (m + len, mm, i, i_bp, d, d_bp, bl + len), // Treat M as a match for statistics
                '=' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                'X' => (m, mm + len, i, i_bp, d, d_bp, bl + len),
                'I' => (m, mm, i + 1, i_bp + len, d, d_bp, bl + len),
                'D' => (m, mm, i, i_bp, d + 1, d_bp + len, bl + len),
                _ => (m, mm, i, i_bp, d, d_bp, bl),
            }
        });
    
    Ok((matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len))
}

/// Calculate identity scores from CIGAR statistics
fn calculate_identity_scores(
    matches: u64, 
    mismatches: u64, 
    insertions: u64, 
    inserted_bp: u64, 
    deletions: u64, 
    deleted_bp: u64
) -> (f64, f64) {
    // Gap-compressed identity
    let gap_compressed_identity = matches as f64 / (matches + mismatches + insertions + deletions) as f64;
    
    // Block identity
    let edit_distance = mismatches + inserted_bp + deleted_bp;
    let block_identity = matches as f64 / (matches + edit_distance) as f64;
    
    (gap_compressed_identity, block_identity)
}

/// Convert a vector of CIGAR operations back to a string
fn cigar_to_string(ops: &[CigarOp]) -> String {
    let mut result = String::new();
    for &CigarOp(op, count) in ops {
        result.push_str(&format!("{}{}", count, op));
    }
    result
}

/// Erode (remove) a specified number of base pairs from the end of a CIGAR string
/// Returns the eroded CIGAR string and the number of bases removed from query and target
fn erode_cigar_end(cigar: &str, bp_count: usize) -> Result<(String, usize, usize), Box<dyn Error>> {
    let ops = parse_cigar(cigar)?;
    let mut new_ops = Vec::new();
    let mut remaining_bp = bp_count;
    let mut query_bases_removed = 0;
    let mut target_bases_removed = 0;
    
    // Copy operations from the beginning, up to the point where we need to erode
    let mut i = 0;
    while i < ops.len() && remaining_bp > 0 {
        // Processing operations from the end, moving backward
        let idx = ops.len() - 1 - i;
        let CigarOp(op, count) = ops[idx];
        
        info!("erode_cigar_end - Processing operation: {}{}", count, op);
        
        match op {
            '=' | 'X' | 'M' => {
                // These operations advance both sequences
                if count > remaining_bp {
                    // Partial erosion of this operation
                    new_ops.insert(0, CigarOp(op, count - remaining_bp));
                    query_bases_removed += remaining_bp;
                    target_bases_removed += remaining_bp;
                    remaining_bp = 0;
                } else {
                    // Complete erosion of this operation
                    query_bases_removed += count;
                    target_bases_removed += count;
                    remaining_bp -= count;
                }
            },
            'I' => {
                // Insertion only advances query
                query_bases_removed += count;
                
                // If we still have remaining_bp to erode, we skip this insertion
                // Otherwise, we keep it
                if remaining_bp > 0 {
                    // Skip this operation (don't add to new_ops)
                } else {
                    new_ops.insert(0, CigarOp(op, count));
                }
            },
            'D' => {
                // Deletion only advances target
                target_bases_removed += count;
                
                // If we still have remaining_bp to erode, we skip this deletion
                // Otherwise, we keep it
                if remaining_bp > 0 {
                    // Skip this operation (don't add to new_ops)
                } else {
                    new_ops.insert(0, CigarOp(op, count));
                }
            },
            _ => new_ops.insert(0, CigarOp(op, count)), // Keep other operations
        }
        
        i += 1;
    }
    
    // Add remaining operations from the beginning of the CIGAR string
    for j in 0..(ops.len() - i) {
        new_ops.insert(j, ops[j]);
    }
    
    Ok((cigar_to_string(&new_ops), query_bases_removed, target_bases_removed))
}

/// Erode (remove) a specified number of base pairs from the start of a CIGAR string
/// Returns the eroded CIGAR string and the number of bases removed from query and target
fn erode_cigar_start(cigar: &str, bp_count: usize) -> Result<(String, usize, usize), Box<dyn Error>> {
    let ops = parse_cigar(cigar)?;
    let mut new_ops = Vec::new();
    let mut remaining_bp = bp_count;
    let mut query_bases_removed = 0;
    let mut target_bases_removed = 0;
    
    // Skip or partially include the first operations until we've eroded enough
    let mut i = 0;
    while i < ops.len() && remaining_bp > 0 {
        let CigarOp(op, count) = ops[i];
        info!("erode_cigar_start - Last operation: {}{}", count, op);
        
        match op {
            '=' | 'X' | 'M' => {
                // These operations advance both sequences
                if count > remaining_bp {
                    // Partial erosion of this operation
                    new_ops.push(CigarOp(op, count - remaining_bp));
                    query_bases_removed += remaining_bp;
                    target_bases_removed += remaining_bp;
                    remaining_bp = 0;
                } else {
                    // Complete erosion of this operation
                    query_bases_removed += count;
                    target_bases_removed += count;
                    remaining_bp -= count;
                }
            },
            'I' => {
                // Insertion only advances query
                query_bases_removed += count;
                
                // If we still have remaining_bp to erode, we skip this insertion
                // Otherwise, we keep it
                if remaining_bp > 0 {
                    // Skip this operation (don't add to new_ops)
                } else {
                    new_ops.push(CigarOp(op, count));
                }
            },
            'D' => {
                // Deletion only advances target
                target_bases_removed += count;
                
                // If we still have remaining_bp to erode, we skip this deletion
                // Otherwise, we keep it
                if remaining_bp > 0 {
                    // Skip this operation (don't add to new_ops)
                } else {
                    new_ops.push(CigarOp(op, count));
                }
            },
            _ => new_ops.push(CigarOp(op, count)), // Keep other operations
        }
        
        i += 1;
    }
    
    // Add all remaining operations
    while i < ops.len() {
        new_ops.push(ops[i]);
        i += 1;
    }
    
    Ok((cigar_to_string(&new_ops), query_bases_removed, target_bases_removed))
}

/// Convert a u8 CIGAR from WFA to a vector of operations with counts
fn cigar_u8_to_cigar_ops(cigar: &[u8]) -> Vec<(usize, char)> {
    if cigar.is_empty() {
        return Vec::new();
    }

    let mut ops = Vec::new();
    let mut count = 0;
    let mut last_op = cigar[0];

    for &op in cigar {
        if op == last_op {
            count += 1;
        } else {
            // Convert the WFA operation to a standard CIGAR operation
            let cigar_op = match last_op {
                b'M' => '=',  // WFA uses M for matches
                b'X' => 'X',  // WFA uses X for mismatches
                b'D' => 'D',  // WFA uses D for deletions
                b'I' => 'I',  // WFA uses I for insertions
                _ => panic!("Unknown CIGAR operation: {}", last_op),
            };
            ops.push((count, cigar_op));
            count = 1;
            last_op = op;
        }
    }

    // Add the last operation
    let cigar_op = match last_op {
        b'M' => '=',
        b'X' => 'X',
        b'D' => 'D',
        b'I' => 'I',
        _ => panic!("Unknown CIGAR operation: {}", last_op),
    };
    ops.push((count, cigar_op));

    ops
}

/// Align two sequence segments using WFA algorithm.
fn align_sequences_wfa(
    a: &[u8],
    b: &[u8],
    aligner: &mut AffineWavefronts
) -> Vec<(usize, char)> {   
    // Do the alignment b vs a (it is not a typo) to have insertions/deletions in the query as Is/Ds in the CIGAR string
    let status = aligner.align(b, a);
    
    match status {
        AlignmentStatus::Completed => {
            cigar_u8_to_cigar_ops(aligner.cigar())
        },
        s => {
            panic!("Alignment failed with status: {:?}", s);
        }
    }
}

/// Perform WFA alignment between query and target sequences
fn wfa_align(query: &[u8], target: &[u8]) -> Result<String, Box<dyn Error>> {
    // Initialize the WFA aligner with gap-affine penalties
    let match_score = 0;
    let mismatch = 3;
    let gap_open1 = 4;
    let gap_ext1 = 2;
    let gap_open2 = 24;
    let gap_ext2 = 1;
    
    // Create aligner and configure settings
    let mut aligner = AffineWavefronts::with_penalties_affine2p(
        match_score, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2
    );
    
    // Align the sequences
    let cigar_ops = align_sequences_wfa(target, query, &mut aligner);
    
    // Convert the CIGAR operations to a CIGAR string
    let cigar = cigar_ops.iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect::<String>();
    
    Ok(cigar)
}

/// Merge multiple CIGAR strings into one
fn merge_cigars(cigars: &[&str]) -> Result<String, Box<dyn Error>> {
    let mut all_ops = Vec::new();
    
    for &cigar in cigars {
        let ops = parse_cigar(cigar)?;
        all_ops.extend_from_slice(&ops);
    }
    
    // Optimize by combining adjacent operations of the same type
    let mut optimized_ops = Vec::new();
    let mut current_op = None;
    
    for op in all_ops {
        if let Some(CigarOp(prev_op, prev_count)) = current_op {
            if prev_op == op.0 {
                current_op = Some(CigarOp(prev_op, prev_count + op.1));
            } else {
                optimized_ops.push(CigarOp(prev_op, prev_count));
                current_op = Some(op);
            }
        } else {
            current_op = Some(op);
        }
    }
    
    if let Some(op) = current_op {
        optimized_ops.push(op);
    }
    
    Ok(cigar_to_string(&optimized_ops))
}

/// Compute coordinates based on CIGAR string
fn compute_paf_coordinates(
    query_start: u64,
    target_start: u64,
    cigar: &str
) -> Result<(u64, u64), Box<dyn Error>> {
    let cigar_ops = parse_cigar(cigar)?;
    
    let mut query_advance = 0u64;
    let mut target_advance = 0u64;
    
    // Calculate how much the query and target advance based on CIGAR operations
    for &CigarOp(op, len) in &cigar_ops {
        let len = len as u64;
        match op {
            '=' | 'X' | 'M' => {
                // These operations advance both sequences
                query_advance += len;
                target_advance += len;
            },
            'I' => {
                // Insertion advances only the query
                query_advance += len;
            },
            'D' => {
                // Deletion advances only the target
                target_advance += len;
            },
            _ => {
                // Other operations (soft clips, etc.) don't advance either sequence
            }
        }
    }
    
    // Compute expected coordinates
    let query_end = query_start + query_advance;
    let target_end = target_start + target_advance;
    
    Ok((query_end, target_end))
}

/// Validate PAF coordinates against CIGAR string
fn validate_paf_coordinates(entry: &PafEntry) -> Result<bool, Box<dyn Error>> {
    let (expected_query_end, expected_target_end) = 
        compute_paf_coordinates(entry.query_start, entry.target_start, &entry.cigar)?;
    
    // Validate coordinates
    let query_valid = expected_query_end == entry.query_end;
    let target_valid = expected_target_end == entry.target_end;
    
    if !query_valid || !target_valid {
        warn!(
            "Coordinate mismatch for {}: Query: expected {} got {}, Target: expected {} got {}", 
            entry.query_name,
            expected_query_end, entry.query_end,
            expected_target_end, entry.target_end
        );
    }
    
    Ok(query_valid && target_valid)
}

/// Ensure all PAF entries in a chain have valid CIGAR strings
fn verify_chain_cigars(chain: &[PafEntry], query_db: &SequenceDB, target_db: &SequenceDB) -> Result<bool, Box<dyn Error>> {
    let mut all_valid = true;

    for (i, entry) in chain.iter().enumerate() {
        // Validate that coordinates match CIGAR length
        let (computed_query_end, computed_target_end) = 
            compute_paf_coordinates(entry.query_start, entry.target_start, &entry.cigar)?;

        if computed_query_end != entry.query_end {
            warn!("Entry {}: Query end mismatch - computed {} but entry has {}",
                  i, computed_query_end, entry.query_end);
            all_valid = false;
        }

        if computed_target_end != entry.target_end {
            warn!("Entry {}: Target end mismatch - computed {} but entry has {}",
                  i, computed_target_end, entry.target_end);
            all_valid = false;
        }
    }

    Ok(all_valid)
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse command line arguments
    let args = Args::parse();
    
    // Initialize logging
    env_logger::Builder::new()
        .filter_level(match args.verbose {
            0 => log::LevelFilter::Warn,    // Errors and warnings
            1 => log::LevelFilter::Info,    // Errors, warnings, and info
            _ => log::LevelFilter::Debug,   // Errors, warnings, info, and debug
        })
        .init();

    // Set number of threads for potential parallel processing
    if args.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()?;
    }
    
    // Load sequence databases with rust-htslib
    info!("Loading query sequences from {}...", args.query.display());
    let query_db = SequenceDB::new(&args.query)?;
    
    info!("Loading target sequences from {}...", args.target.display());
    let target_db = SequenceDB::new(&args.target)?;
    
    // Read PAF file
    info!("Reading PAF file {}...", args.paf.display());
    let entries = read_paf_file(&args.paf)?;
    info!("Read {} PAF entries", entries.len());
    
    // Group entries by chain_id
    let chain_groups = group_by_chain(entries);
    info!("Found {} chains", chain_groups.len());
    
    // Process each chain
    let mut processed_entries = Vec::new();
    
    for (chain_id, mut group) in chain_groups {
        info!("Processing chain {} with {} entries", chain_id, group.len());
        
        // Verify CIGAR strings before processing
        info!("Verifying CIGAR strings for chain {}", chain_id);
        match verify_chain_cigars(&group, &query_db, &target_db) {
            Ok(true) => info!("All CIGAR strings valid for chain {}", chain_id),
            Ok(false) => warn!("Some CIGAR strings invalid for chain {}", chain_id),
            Err(e) => warn!("Error verifying CIGAR strings: {}", e),
        }
        
        // Sort by chain position
        group.sort_by_key(|entry| entry.chain_pos);
        
        // Process chain and add results to processed_entries
        let chain_result = process_chain(&group, &query_db, &target_db, args.erosion_size)?;
        
        // Verify the processed chain
        info!("Verifying processed chain {}", chain_id);
        match verify_chain_cigars(&chain_result, &query_db, &target_db) {
            Ok(true) => info!("Processed chain {} has valid CIGAR strings", chain_id),
            Ok(false) => warn!("Processed chain {} has invalid CIGAR strings", chain_id),
            Err(e) => warn!("Error verifying processed chain: {}", e),
        }
        
        processed_entries.extend(chain_result);
    }
    
    // Write output
    let output_path = args.output.unwrap_or_else(|| PathBuf::from("output.paf"));
    info!("Writing {} processed entries to {}", processed_entries.len(), output_path.display());
    write_paf_entries(&output_path, &processed_entries)?;
    
    info!("Done!");
    Ok(())
}

fn read_paf_file<P: AsRef<Path>>(path: P) -> Result<Vec<PafEntry>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();
    
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        match PafEntry::parse_from_line(&line) {
            Ok(entry) => entries.push(entry),
            Err(e) => warn!("Failed to parse PAF line {}: {}", i + 1, e),
        }
    }
    
    Ok(entries)
}

fn group_by_chain(entries: Vec<PafEntry>) -> HashMap<u64, Vec<PafEntry>> {
    let mut chain_groups: HashMap<u64, Vec<PafEntry>> = HashMap::new();
    
    for entry in entries {
        chain_groups.entry(entry.chain_id).or_insert_with(Vec::new).push(entry);
    }
    
    chain_groups
}

fn process_chain(
    chain: &[PafEntry], 
    query_db: &SequenceDB, 
    target_db: &SequenceDB,
    erosion_size: usize
) -> Result<Vec<PafEntry>, Box<dyn Error>> {
    if chain.len() <= 1 {
        // If only one entry, return it unchanged
        return Ok(chain.to_vec());
    }
    
    // Calculate the average md:f value from all entries in the chain
    let mut total_md = 0.0;
    let mut count_md = 0;
    
    for entry in chain {
        // Find the md:f tag if it exists
        for (tag, value) in &entry.tags {
            if tag == "md:f" {
                if let Ok(md_value) = value.parse::<f64>() {
                    total_md += md_value;
                    count_md += 1;
                }
                break;
            }
        }
    }
    
    // Calculate average md:f value
    let avg_md = if count_md > 0 {
        total_md / count_md as f64
    } else {
        // Default value if no md:f tags found
        0.0
    };
    
    // Format average md:f value
    let avg_md_str = format!("{:.6}", avg_md)
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string();
    
    info!("Average md:f value for chain: {}", avg_md_str);
    
    let mut merged_entry = chain[0].clone();
    
    for i in 0..(chain.len() - 1) {
        // Create a copy of the current entry to avoid borrowing issues
        let current = if i == 0 { chain[0].clone() } else { merged_entry };
        let next = &chain[i+1];
        
        info!("Processing chain pair {}/{}", i+1, chain.len()-1);
        
        // Get the CIGAR strings
        let current_cigar = &current.cigar;
        let next_cigar = &next.cigar;
        
        // Erode CIGARs at the boundaries and get the removed bases counts
        let (eroded_current_cigar, current_query_bases_removed, current_target_bases_removed) = 
            erode_cigar_end(current_cigar, erosion_size)?;
        
        let (eroded_next_cigar, next_query_bases_removed, next_target_bases_removed) = 
            erode_cigar_start(next_cigar, erosion_size)?;
        
        info!("Eroded current: query={}, target={} bases", 
              current_query_bases_removed, current_target_bases_removed);
        info!("Eroded next: query={}, target={} bases", 
              next_query_bases_removed, next_target_bases_removed);
        
        // Calculate correct coordinates for sequence extraction based on actual bases removed
        let current_query_extract_start = current.query_end - current_query_bases_removed as u64;
        let current_target_extract_start = current.target_end - current_target_bases_removed as u64;
        
        let next_query_extract_end = next.query_start + next_query_bases_removed as u64;
        let next_target_extract_end = next.target_start + next_target_bases_removed as u64;
        
        // Check if we have continuity in coordinates
        let query_gap = if next.query_start > current.query_end {
            next.query_start - current.query_end
        } else {
            0 // Overlapping or exactly contiguous
        };
        
        let target_gap = if next.target_start > current.target_end {
            next.target_start - current.target_end
        } else {
            0 // Overlapping or exactly contiguous
        };
        
        info!("Query gap: {}, Target gap: {}", query_gap, target_gap);
        
        // Get sequences for boundary regions using rust-htslib, with correct coordinates
        let current_query_seq = query_db.get_subsequence(
            &current.query_name,
            current_query_extract_start,
            current.query_end,
            current.strand == '-'
        ).map_err(|e| format!("Failed to get query sequence for {}: {}", current.query_name, e))?;
        
        let next_query_seq = query_db.get_subsequence(
            &next.query_name,
            next.query_start,
            next_query_extract_end,
            next.strand == '-'
        ).map_err(|e| format!("Failed to get query sequence for {}: {}", next.query_name, e))?;
        
        let current_target_seq = target_db.get_subsequence(
            &current.target_name,
            current_target_extract_start,
            current.target_end,
            false
        ).map_err(|e| format!("Failed to get target sequence for {}: {}", current.target_name, e))?;
        
        let next_target_seq = target_db.get_subsequence(
            &next.target_name,
            next.target_start,
            next_target_extract_end,
            false
        ).map_err(|e| format!("Failed to get target sequence for {}: {}", next.target_name, e))?;
        
        // Log sequence lengths for debugging
        info!("Current query seq length: {}, Current target seq length: {}", 
              current_query_seq.len(), current_target_seq.len());
        info!("Next query seq length: {}, Next target seq length: {}", 
              next_query_seq.len(), next_target_seq.len());
        
        // If there's a gap, also get the gap sequences
        let mut gap_query_seq = Vec::new();
        let mut gap_target_seq = Vec::new();
        
        if query_gap > 0 {
            let gap_query = query_db.get_subsequence(
                &current.query_name,
                current.query_end,
                next.query_start,
                current.strand == '-'
            ).map_err(|e| format!("Failed to get gap query sequence: {}", e))?;
            gap_query_seq = gap_query;
            info!("Got gap query sequence of length {}", gap_query_seq.len());
        }
        
        if target_gap > 0 {
            let gap_target = target_db.get_subsequence(
                &current.target_name,
                current.target_end,
                next.target_start,
                false
            ).map_err(|e| format!("Failed to get gap target sequence: {}", e))?;
            gap_target_seq = gap_target;
            info!("Got gap target sequence of length {}", gap_target_seq.len());
        }
        
        // Combine sequences to create a complete sequence
        let mut full_query = Vec::new();
        full_query.extend_from_slice(&current_query_seq);
        if !gap_query_seq.is_empty() {
            full_query.extend_from_slice(&gap_query_seq);
        }
        full_query.extend_from_slice(&next_query_seq);
        
        let mut full_target = Vec::new();
        full_target.extend_from_slice(&current_target_seq);
        if !gap_target_seq.is_empty() {
            full_target.extend_from_slice(&gap_target_seq);
        }
        full_target.extend_from_slice(&next_target_seq);
        
        // Perform WFA alignment on the combined sequences
        info!("Performing WFA alignment for combined regions...");
        info!("Full query length: {}, Full target length: {}", full_query.len(), full_target.len());
        let combined_cigar = wfa_align(&full_query, &full_target)?;
        
        // Merge the CIGARs
        let merged_cigar = merge_cigars(&[&eroded_current_cigar, &combined_cigar, &eroded_next_cigar])?;
        
        // Calculate alignment statistics from the merged CIGAR
        let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) = 
            calculate_cigar_stats(&merged_cigar)?;
        
        // Calculate identity scores
        let (gap_compressed_identity, block_identity) = 
            calculate_identity_scores(matches, mismatches, insertions, inserted_bp, deletions, deleted_bp);
        
        // Format identity scores
        let gi_str = format!("{:.6}", gap_compressed_identity)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
        let bi_str = format!("{:.6}", block_identity)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
        
        // Compute correct end coordinates based on CIGAR
        let (computed_query_end, computed_target_end) = 
            compute_paf_coordinates(current.query_start, current.target_start, &merged_cigar)?;
        
        // Log if there's a mismatch with the expected coordinates
        if computed_query_end != next.query_end {
            warn!("Query end coordinate mismatch: computed {} vs expected {} (from next entry)",
                  computed_query_end, next.query_end);
        }
        if computed_target_end != next.target_end {
            warn!("Target end coordinate mismatch: computed {} vs expected {} (from next entry)",
                  computed_target_end, next.target_end);
        }
        
        // Create a map of the current tag order by position
        let mut tag_positions = HashMap::new();
        for (i, (tag, _)) in current.tags.iter().enumerate() {
            tag_positions.insert(tag.as_str(), i);
        }
        
        // Create a map for new tag values
        let mut new_tag_values = HashMap::new();
        for (tag, value) in &current.tags {
            if !matches!(tag.as_str(), "gi:f" | "bi:f" | "cg:Z" | "md:f" | "chain:i") {
                new_tag_values.insert(tag.as_str(), value.clone());
            }
        }
        
        // Update with new values
        new_tag_values.insert("gi:f", gi_str);
        new_tag_values.insert("bi:f", bi_str);
        new_tag_values.insert("cg:Z", merged_cigar.clone());
        new_tag_values.insert("md:f", avg_md_str.clone());
        
        // Fix chain:i tag to have correct chain_len and chain_pos values
        // Format: chain_id.chain_len.chain_pos
        let chain_tag_value = format!("{}.1.1", current.chain_id);
        new_tag_values.insert("chain:i", chain_tag_value);
        
        // Create ordered tags vector preserving original order
        let mut new_tags = Vec::new();
        
        // First add any tags in their original order
        for (tag, _) in &current.tags {
            if let Some(value) = new_tag_values.get(tag.as_str()) {
                new_tags.push((tag.clone(), value.clone()));
                new_tag_values.remove(tag.as_str());
            }
        }
        
        // Add any remaining tags (should be none if properly preserved)
        for (tag, value) in new_tag_values {
            new_tags.push((tag.to_string(), value));
        }
        
        // Update the merged entry with computed coordinates
        merged_entry = PafEntry {
            query_name: current.query_name,
            query_length: current.query_length,
            query_start: current.query_start,
            query_end: computed_query_end,
            strand: current.strand,
            target_name: current.target_name,
            target_length: current.target_length,
            target_start: current.target_start,
            target_end: computed_target_end,
            num_matches: matches,
            alignment_length: block_len,
            mapping_quality: 255, // High mapping quality for merged alignments
            tags: new_tags,
            chain_id: current.chain_id,
            chain_length: 1,  // Set chain_length to 1
            chain_pos: 1,     // Set chain_pos to 1
            cigar: merged_cigar,
        };
    }
    
    // Final validation of the merged entry
    let (final_query_end, final_target_end) = 
        compute_paf_coordinates(merged_entry.query_start, merged_entry.target_start, &merged_entry.cigar)?;
    
    if final_query_end != merged_entry.query_end || final_target_end != merged_entry.target_end {
        warn!("Final coordinate validation failed: Query: computed {} vs stored {}, Target: computed {} vs stored {}",
                final_query_end, merged_entry.query_end, 
                final_target_end, merged_entry.target_end);
    } else {
        info!("Final coordinate validation successful");
    }
    
    Ok(vec![merged_entry])
}

fn write_paf_entries<P: AsRef<Path>>(path: P, entries: &[PafEntry]) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = io::BufWriter::new(file);
    
    for entry in entries {
        writeln!(writer, "{}", entry.to_string())?;
    }
    
    Ok(())
}