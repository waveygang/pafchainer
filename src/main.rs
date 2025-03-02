use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::error::Error;
use clap::Parser;
use bio::io::fasta;
use log::{info, warn, error};
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, HeuristicStrategy};

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
    chain_id: Option<u64>,
    chain_length: Option<u64>,
    chain_pos: Option<u64>,
    cigar: Option<String>,
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
            chain_id: None,
            chain_length: None,
            chain_pos: None,
            cigar: None,
        };

        // Parse optional tags
        for i in 12..fields.len() {
            let tag_parts: Vec<&str> = fields[i].split(':').collect();
            if tag_parts.len() >= 3 {
                let tag_name = format!("{}:{}", tag_parts[0], tag_parts[1]);
                let tag_value = tag_parts[2..].join(":");
                
                // Parse chain information
                if tag_name == "chain:i" {
                    let chain_parts: Vec<&str> = tag_value.split('.').collect();
                    if chain_parts.len() >= 3 {
                        entry.chain_id = Some(chain_parts[0].parse()?);
                        entry.chain_length = Some(chain_parts[1].parse()?);
                        entry.chain_pos = Some(chain_parts[2].parse()?);
                    }
                }
                
                // Parse CIGAR string
                if tag_name == "cg:Z" {
                    entry.cigar = Some(tag_value.clone());
                }
                
                entry.tags.push((tag_name, tag_value));
            }
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
        
        // Add tags
        for (tag, value) in &self.tags {
            line.push_str(&format!("\t{}:{}", tag, value));
        }
        
        line
    }
}

/// Sequence database for holding FASTA sequences
struct SequenceDB {
    sequences: HashMap<String, Vec<u8>>,
}

impl SequenceDB {
    fn new() -> Self {
        Self {
            sequences: HashMap::new(),
        }
    }
    
    fn load_fasta<P: AsRef<Path>>(&mut self, path: P) -> Result<(), Box<dyn Error>> {
        let reader = fasta::Reader::new(File::open(path)?);
        
        for record in reader.records() {
            let record = record?;
            self.sequences.insert(record.id().to_string(), record.seq().to_vec());
        }
        
        Ok(())
    }
    
    fn get_subsequence(&self, name: &str, start: u64, end: u64, reverse_complement: bool) -> Option<Vec<u8>> {
        self.sequences.get(name).map(|seq| {
            let start = start as usize;
            let end = std::cmp::min(end as usize, seq.len());
            
            if start >= seq.len() {
                return Vec::new();
            }
            
            let subseq = seq[start..end].to_vec();
            
            if reverse_complement {
                // Implement reverse complement
                subseq.iter().rev().map(|&b| match b {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'G' => b'C',
                    b'C' => b'G',
                    x => x,
                }).collect()
            } else {
                subseq
            }
        })
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

/// Convert a vector of CIGAR operations back to a string
fn cigar_to_string(ops: &[CigarOp]) -> String {
    let mut result = String::new();
    for &CigarOp(op, count) in ops {
        result.push_str(&format!("{}{}", count, op));
    }
    result
}

/// Erode (remove) a specified number of base pairs from the end of a CIGAR string
fn erode_cigar_end(cigar: &str, bp_count: usize) -> Result<String, Box<dyn Error>> {
    let ops = parse_cigar(cigar)?;
    let mut new_ops = Vec::new();
    let mut remaining_bp = bp_count;
    
    for i in 0..ops.len() {
        let CigarOp(op, count) = ops[i];
        
        if remaining_bp == 0 || i < ops.len() - 1 {
            // Keep operations that are not at the end
            new_ops.push(CigarOp(op, count));
            continue;
        }
        
        // Process the last operation
        match op {
            '=' | 'X' | 'M' | 'I' | 'D' => {
                if count > remaining_bp {
                    new_ops.push(CigarOp(op, count - remaining_bp));
                    remaining_bp = 0;
                } else {
                    remaining_bp -= count;
                }
            },
            _ => new_ops.push(CigarOp(op, count)), // Keep other operations
        }
    }
    
    Ok(cigar_to_string(&new_ops))
}

/// Erode (remove) a specified number of base pairs from the start of a CIGAR string
fn erode_cigar_start(cigar: &str, bp_count: usize) -> Result<String, Box<dyn Error>> {
    let ops = parse_cigar(cigar)?;
    let mut new_ops = Vec::new();
    let mut remaining_bp = bp_count;
    
    for i in (0..ops.len()).rev() {
        let CigarOp(op, count) = ops[i];
        
        if remaining_bp == 0 || i > 0 {
            // Keep operations that are not at the start
            new_ops.push(CigarOp(op, count));
            continue;
        }
        
        // Process the first operation
        match op {
            '=' | 'X' | 'M' | 'I' | 'D' => {
                if count > remaining_bp {
                    new_ops.push(CigarOp(op, count - remaining_bp));
                    remaining_bp = 0;
                } else {
                    remaining_bp -= count;
                }
            },
            _ => new_ops.push(CigarOp(op, count)), // Keep other operations
        }
    }
    
    new_ops.reverse(); // Reverse since we processed from the end
    Ok(cigar_to_string(&new_ops))
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
    let mismatch = 4;
    let gap_open1 = 6;
    let gap_ext1 = 2;
    let gap_open2 = 6;
    let gap_ext2 = 2;
    
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

fn main() -> Result<(), Box<dyn Error>> {
    // Initialize logging
    env_logger::init();
    
    // Parse command line arguments
    let args = Args::parse();
    
    // Load sequence databases
    info!("Loading query sequences from {}...", args.query.display());
    let mut query_db = SequenceDB::new();
    query_db.load_fasta(&args.query)?;
    
    info!("Loading target sequences from {}...", args.target.display());
    let mut target_db = SequenceDB::new();
    target_db.load_fasta(&args.target)?;
    
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
        
        // Sort by chain position
        group.sort_by_key(|entry| entry.chain_pos.unwrap_or(0));
        
        // Process chain and add results to processed_entries
        let chain_result = process_chain(&group, &query_db, &target_db, args.erosion_size)?;
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
        if let Some(chain_id) = entry.chain_id {
            chain_groups.entry(chain_id).or_insert_with(Vec::new).push(entry);
        } else {
            // For entries without chain_id, place each in its own group
            let unique_id = chain_groups.len() as u64 + 10000; // Use a high number to avoid collisions
            chain_groups.entry(unique_id).or_insert_with(Vec::new).push(entry);
        }
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
    
    let mut merged_entry = chain[0].clone();
    
    for i in 0..(chain.len() - 1) {
        let current = if i == 0 { &chain[0] } else { &merged_entry };
        let next = &chain[i+1];
        
        info!("Processing chain pair {}/{}", i+1, chain.len()-1);
        
        // Skip if no CIGAR
        if current.cigar.is_none() || next.cigar.is_none() {
            warn!("Missing CIGAR in entry, cannot merge");
            return Ok(chain.to_vec());
        }
        
        // Get the CIGAR strings
        let current_cigar = current.cigar.as_ref().unwrap();
        let next_cigar = next.cigar.as_ref().unwrap();
        
        // Erode CIGARs at the boundaries
        let eroded_current_cigar = erode_cigar_end(current_cigar, erosion_size)?;
        let eroded_next_cigar = erode_cigar_start(next_cigar, erosion_size)?;
        
        // Get sequences for boundary regions
        let current_query_seq = query_db.get_subsequence(
            &current.query_name,
            current.query_end - erosion_size as u64,
            current.query_end,
            current.strand == '-'
        ).ok_or_else(|| format!("Failed to get query sequence for {}", current.query_name))?;
        
        let next_query_seq = query_db.get_subsequence(
            &next.query_name,
            next.query_start,
            next.query_start + erosion_size as u64,
            next.strand == '-'
        ).ok_or_else(|| format!("Failed to get query sequence for {}", next.query_name))?;
        
        let current_target_seq = target_db.get_subsequence(
            &current.target_name,
            current.target_end - erosion_size as u64,
            current.target_end,
            false
        ).ok_or_else(|| format!("Failed to get target sequence for {}", current.target_name))?;
        
        let next_target_seq = target_db.get_subsequence(
            &next.target_name,
            next.target_start,
            next.target_start + erosion_size as u64,
            false
        ).ok_or_else(|| format!("Failed to get target sequence for {}", next.target_name))?;
        
        // Combine sequences to create a gap sequence
        let mut gap_query = Vec::new();
        gap_query.extend_from_slice(&current_query_seq);
        gap_query.extend_from_slice(&next_query_seq);
        
        let mut gap_target = Vec::new();
        gap_target.extend_from_slice(&current_target_seq);
        gap_target.extend_from_slice(&next_target_seq);
        
        // Perform WFA alignment on the gap
        info!("Performing WFA alignment for gap region...");
        let gap_cigar = wfa_align(&gap_query, &gap_target)?;
        
        // Merge the CIGARs
        let merged_cigar = merge_cigars(&[&eroded_current_cigar, &gap_cigar, &eroded_next_cigar])?;
        
        // Update the merged entry
        merged_entry = PafEntry {
            query_name: current.query_name.clone(),
            query_length: current.query_length,
            query_start: current.query_start,
            query_end: next.query_end,
            strand: current.strand,
            target_name: current.target_name.clone(),
            target_length: current.target_length,
            target_start: current.target_start,
            target_end: next.target_end,
            num_matches: current.num_matches + next.num_matches, // Approximate
            alignment_length: next.query_end - current.query_start,
            mapping_quality: current.mapping_quality,
            tags: current.tags.clone(),
            chain_id: current.chain_id,
            chain_length: current.chain_length,
            chain_pos: current.chain_pos,
            cigar: Some(merged_cigar.clone()),
        };
        
        // Update the CIGAR in tags
        for (tag, value) in &mut merged_entry.tags {
            if tag == "cg:Z" {
                *value = merged_cigar.clone();
                break;
            }
        }
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