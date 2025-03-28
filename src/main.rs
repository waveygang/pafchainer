use std::collections::{HashSet, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::error::Error;
use clap::Parser;
use log::{info, warn, debug};
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus};
use rust_htslib::faidx::Reader as FastaReader;
use libc;

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
    #[clap(short, long, default_value = "100")]
    erosion_size: usize,
    
    /// Output PAF file
    #[clap(short, long)]
    output: Option<PathBuf>,
    
    /// Output in SAM format instead of PAF
    #[clap(long)]
    sam: bool,

    /// Number of threads to use
    #[clap(short, long, default_value = "4")]
    threads: usize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[clap(short, long, default_value = "0")]
    verbose: u8,
}

// PAF entry representation
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
    tags: Vec<(String, String)>,
    chain_id: u64,
    chain_length: u64,
    chain_pos: u64,
    cigar: String,
}

// CIGAR operation
#[derive(Debug, Clone, Copy)]
struct CigarOp(char, usize);

// Sequence database to handle FASTA files
struct SequenceDB {
    reader: FastaReader,
}

impl SequenceDB {
    // Create a new sequence database
    fn new<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        Ok(Self {
            reader: FastaReader::from_path(path)?,
        })
    }
    
    // Get a subsequence from the database
    fn get_subsequence(&self, name: &str, start: u64, end: u64, reverse_complement: bool) -> Result<Vec<u8>, Box<dyn Error>> {
        let start_usize = start as usize;
        let end_usize = std::cmp::min(end as usize, std::usize::MAX);
        
        // Fetch sequence and properly handle memory
        let seq = match self.reader.fetch_seq(name, start_usize, end_usize - 1) {
            Ok(seq) => {
                let mut seq_vec = seq.to_vec();
                unsafe {libc::free(seq.as_ptr() as *mut std::ffi::c_void)}; // Free up memory to avoid memory leak (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                seq_vec
            },
            Err(e) => return Err(format!("Failed to fetch sequence for {}: {}", name, e).into()),
        };
        
        if reverse_complement {
            // Reverse complement the sequence
            Ok(seq.iter().rev().map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                x => x,
            }).collect())
        } else {
            Ok(seq)
        }
    }
}

impl PafEntry {
    // Parse a PAF entry from a line
    fn parse_from_line(line: &str) -> Result<Self, Box<dyn Error>> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err("PAF line has fewer than 12 required fields".into());
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

        // Parse optional tags
        for i in 12..fields.len() {
            let tag_parts: Vec<&str> = fields[i].split(':').collect();
            if tag_parts.len() >= 3 {
                let tag_name = format!("{}:{}", tag_parts[0], tag_parts[1]);
                let tag_value = tag_parts[2..].join(":");
                
                // Parse chain information
                if tag_name == "ch:Z" {
                    let chain_parts: Vec<&str> = tag_value.split('.').collect();
                    if chain_parts.len() >= 3 {
                        entry.chain_id = chain_parts[0].parse()?;
                        entry.chain_length = chain_parts[1].parse()?;
                        entry.chain_pos = chain_parts[2].parse()?;
                    } else {
                        return Err(format!("Invalid ch:Z tag format: {}", tag_value).into());
                    }
                }
                
                // Parse CIGAR string
                if tag_name == "cg:Z" {
                    entry.cigar = tag_value.clone();
                }
                
                entry.tags.push((tag_name, tag_value));
            }
        }

        // Verify required tags
        if entry.chain_id == 0 || entry.cigar.is_empty() {
            return Err("Missing required 'ch:Z' or 'cg:Z' tag".into());
        }

        Ok(entry)
    }
    
    // Convert PAF entry to string
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

// Parse CIGAR string into operations
fn cigar_to_ops(cigar: &str) -> Result<Vec<CigarOp>, Box<dyn Error>> {
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

// Convert WFA alignment CIGAR to standard operations
fn cigar_u8_to_ops(cigar: &[u8]) -> Vec<(usize, char)> {
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
            // Convert the operation
            let cigar_op = match last_op {
                b'M' => '=',  // Match
                b'X' => 'X',  // Mismatch
                b'D' => 'D',  // Deletion
                b'I' => 'I',  // Insertion
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

// Convert CIGAR operations to string
fn ops_to_cigar(ops: &[CigarOp]) -> String {
    ops.iter().map(|&CigarOp(op, count)| format!("{}{}", count, op)).collect()
}

// Erode CIGAR from the end with separate query and target erosion limits
fn erode_cigar_end(cigar: &str, query_bp_limit: usize, target_bp_limit: usize) -> Result<(String, usize, usize), Box<dyn Error>> {
    let ops = cigar_to_ops(cigar)?;
    let mut new_ops = Vec::new();
    let mut remaining_query_bp = query_bp_limit;
    let mut remaining_target_bp = target_bp_limit;
    let mut query_bases_removed = 0;
    let mut target_bases_removed = 0;
    
    debug!("\tEroding from end of CIGAR - Query limit: {} bp, Target limit: {} bp", 
          query_bp_limit, target_bp_limit);
    
    // Copy operations from beginning until we need to erode
    let mut i = 0;
    while i < ops.len() && (remaining_query_bp > 0 || remaining_target_bp > 0) {
        // Process operations from the end
        let idx = ops.len() - 1 - i;
        let CigarOp(op, count) = ops[idx];
        
        // debug!("\t\terode_cigar_end - Processing operation: {}{}", count, op);
        
        match op {
            '=' | 'X' | 'M' => {
                // Advances both sequences
                if remaining_query_bp > 0 || remaining_target_bp > 0 {
                    let query_erosion = std::cmp::min(count, remaining_query_bp);
                    let target_erosion = std::cmp::min(count, remaining_target_bp);
                    
                    // Determine actual erosion
                    let actual_erosion = std::cmp::max(query_erosion, target_erosion);
                    
                    if actual_erosion < count {
                        // Partial erosion
                        new_ops.insert(0, CigarOp(op, count - actual_erosion));
                    }
                    
                    query_bases_removed += actual_erosion;
                    target_bases_removed += actual_erosion;
                    remaining_query_bp = remaining_query_bp.saturating_sub(actual_erosion);
                    remaining_target_bp = remaining_target_bp.saturating_sub(actual_erosion);
                } else {
                    // No more erosion needed
                    new_ops.insert(0, CigarOp(op, count));
                }
            },
            'I' => {
                // Insertion only advances query
                if remaining_query_bp > 0 {
                    let actual_erosion = count;
                    
                    query_bases_removed += actual_erosion;
                    remaining_query_bp = remaining_query_bp.saturating_sub(actual_erosion);
                } else {
                    // No more query erosion needed
                    new_ops.insert(0, CigarOp(op, count));
                }
            },
            'D' => {
                // Deletion only advances target
                if remaining_target_bp > 0 {
                    let actual_erosion = count;

                    target_bases_removed += actual_erosion;
                    remaining_target_bp = remaining_target_bp.saturating_sub(actual_erosion);
                } else {
                    // No more target erosion needed
                    new_ops.insert(0, CigarOp(op, count));
                }
            },
            _ => new_ops.insert(0, CigarOp(op, count)), // Keep other operations
        }
        
        i += 1;
        
        // debug!("\t\tProgress - Eroded from query: {}/{}, target: {}/{}", 
        //       query_bases_removed, query_bp_limit, 
        //       target_bases_removed, target_bp_limit);
    }
    
    // Add remaining operations
    for j in 0..(ops.len() - i) {
        new_ops.insert(j, ops[j]);
    }
    
    Ok((ops_to_cigar(&new_ops), query_bases_removed, target_bases_removed))
}

// Erode CIGAR from the start with separate query and target erosion limits
fn erode_cigar_start(cigar: &str, query_bp_limit: usize, target_bp_limit: usize) -> Result<(String, usize, usize), Box<dyn Error>> {
    let ops = cigar_to_ops(cigar)?;
    let mut new_ops = Vec::new();
    let mut remaining_query_bp = query_bp_limit;
    let mut remaining_target_bp = target_bp_limit;
    let mut query_bases_removed = 0;
    let mut target_bases_removed = 0;
    
    debug!("\tEroding from start of CIGAR - Query limit: {} bp, Target limit: {} bp", 
          query_bp_limit, target_bp_limit);
    
    // Skip or partially include first operations until eroded enough
    let mut i = 0;
    while i < ops.len() && (remaining_query_bp > 0 || remaining_target_bp > 0) {
        let CigarOp(op, count) = ops[i];
        
        // debug!("\t\terode_cigar_start - Processing operation: {}{}", count, op);
        
        match op {
            '=' | 'X' | 'M' => {
                // Advances both sequences
                if remaining_query_bp > 0 || remaining_target_bp > 0 {
                    let query_erosion = std::cmp::min(count, remaining_query_bp);
                    let target_erosion = std::cmp::min(count, remaining_target_bp);
                    
                    // Determine actual erosion
                    let actual_erosion = std::cmp::max(query_erosion, target_erosion);

                    if actual_erosion < count {
                        // Partial erosion
                        new_ops.push(CigarOp(op, count - actual_erosion));
                    }
                    
                    query_bases_removed += actual_erosion;
                    target_bases_removed += actual_erosion;
                    remaining_query_bp = remaining_query_bp.saturating_sub(actual_erosion);
                    remaining_target_bp = remaining_target_bp.saturating_sub(actual_erosion);
                } else {
                    // No more erosion needed
                    new_ops.push(CigarOp(op, count));
                }
            },
            'I' => {
                // Insertion only advances query
                if remaining_query_bp > 0 {
                    let actual_erosion = std::cmp::min(count, remaining_query_bp);
                    
                    if actual_erosion < count {
                        // Partial erosion
                        new_ops.push(CigarOp(op, count - actual_erosion));
                    }
                    
                    query_bases_removed += actual_erosion;
                    remaining_query_bp = remaining_query_bp.saturating_sub(actual_erosion);
                } else {
                    // No more query erosion needed
                    new_ops.push(CigarOp(op, count));
                }
            },
            'D' => {
                // Deletion only advances target
                if remaining_target_bp > 0 {
                    let actual_erosion = std::cmp::min(count, remaining_target_bp);
                    
                    if actual_erosion < count {
                        // Partial erosion
                        new_ops.push(CigarOp(op, count - actual_erosion));
                    }
                    
                    target_bases_removed += actual_erosion;
                    remaining_target_bp = remaining_target_bp.saturating_sub(actual_erosion);
                } else {
                    // No more target erosion needed
                    new_ops.push(CigarOp(op, count));
                }
            },
            _ => new_ops.push(CigarOp(op, count)), // Keep other operations
        }
        
        i += 1;
        
        // debug!("\t\tProgress - Eroded from query: {}/{}, target: {}/{}", 
        //       query_bases_removed, query_bp_limit, 
        //       target_bases_removed, target_bp_limit);
    }
    
    // Add all remaining operations
    while i < ops.len() {
        new_ops.push(ops[i]);
        i += 1;
    }
    
    Ok((ops_to_cigar(&new_ops), query_bases_removed, target_bases_removed))
}

// Align two sequences using WFA
fn align_sequences_wfa(query: &[u8], target: &[u8]) -> Result<String, Box<dyn Error>> {
    debug!("Performing WFA alignment between sequences of lengths {} and {}", query.len(), target.len());
    
    // Initialize WFA aligner with gap-affine penalties
    let match_score = 0;
    let mismatch = 3;
    let gap_open1 = 4;
    let gap_ext1 = 2;
    let gap_open2 = 24;
    let gap_ext2 = 1;
    
    let aligner = AffineWavefronts::with_penalties_affine2p(
        match_score, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2
    );
    
    // Perform alignment (note that WFA expects target, query order)
    let status = aligner.align(target, query);
    
    // Check alignment status (can't use != since AlignmentStatus doesn't implement PartialEq)
    match status {
        AlignmentStatus::Completed => {
            debug!("WFA alignment completed successfully");
        },
        s => {
            warn!("Alignment failed with status: {:?}", s);
            return Err(format!("Alignment failed with status: {:?}", s).into());
        }
    }
    
    // Convert CIGAR to string
    let cigar_ops = cigar_u8_to_ops(aligner.cigar());
    let cigar = cigar_ops.iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect::<String>();
    
    Ok(cigar)
}

// Merge multiple CIGAR strings
fn merge_cigars(cigars: &[&str]) -> Result<String, Box<dyn Error>> {
    let mut all_ops = Vec::new();
    
    // Collect all operations
    for &cigar in cigars {
        let ops = cigar_to_ops(cigar)?;
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
    
    Ok(ops_to_cigar(&optimized_ops))
}

// Calculate CIGAR statistics
fn calculate_cigar_stats(cigar: &str) -> Result<(u64, u64, u64, u64, u64, u64, u64), Box<dyn Error>> {
    let ops = cigar_to_ops(cigar)?;
    let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) = 
        ops.iter().fold((0, 0, 0, 0, 0, 0, 0), |(m, mm, i, i_bp, d, d_bp, bl), &CigarOp(op, len)| {
            let len = len as u64;
            match op {
                'M' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                '=' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                'X' => (m, mm + len, i, i_bp, d, d_bp, bl + len),
                'I' => (m, mm, i + 1, i_bp + len, d, d_bp, bl + len),
                'D' => (m, mm, i, i_bp, d + 1, d_bp + len, bl + len),
                _ => (m, mm, i, i_bp, d, d_bp, bl),
            }
        });
    
    Ok((matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len))
}

// Read PAF file
fn read_paf_file<P: AsRef<Path>>(path: P) -> Result<Vec<PafEntry>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();
    
    for line in reader.lines() {
        if let Ok(entry) = PafEntry::parse_from_line(&line?) {
            entries.push(entry);
        }
    }
    
    Ok(entries)
}

// Write PAF entries to a writer
fn write_paf_content<W: Write>(writer: &mut W, entries: &[PafEntry]) -> Result<(), Box<dyn Error>> {
    for entry in entries {
        writeln!(writer, "{}", entry.to_string())?;
    }
    
    Ok(())
}

// Generate SAM header
fn generate_sam_header(entries: &[PafEntry]) -> String {
    let mut header = String::new();
    
    // Add SAM version header
    header.push_str("@HD\tVN:1.6\tSO:unsorted\n");
    
    // Add reference sequence information
    let mut target_seen = HashSet::new();
    for entry in entries {
        if !target_seen.contains(&entry.target_name) {
            header.push_str(&format!(
                "@SQ\tSN:{}\tLN:{}\n", 
                entry.target_name, entry.target_length
            ));
            target_seen.insert(entry.target_name.clone());
        }
    }
    
    // Add program information
    header.push_str("@PG\tID:chain_connect\tPN:chain_connect\tVN:0.1\n");
    
    header
}

// Convert PafEntry to SAM format
fn paf_entry_to_sam(
    entry: &PafEntry, 
    query_db: &SequenceDB,
    reference_id_map: &HashMap<String, usize>
) -> Result<String, Box<dyn Error>> {
    // Calculate flag
    let flag = if entry.strand == '-' { 16 } else { 0 };
    
    // Get reference ID
    let ref_id = match reference_id_map.get(&entry.target_name) {
        Some(_) => entry.target_name.clone(),
        None => "*".to_string(),
    };

    let sam_cigar = &entry.cigar;
    
    // Fetch sequence
    let sequence = query_db.get_subsequence(
        &entry.query_name,
        entry.query_start,
        entry.query_end,
        entry.strand == '-'
    )?;
    
    // Convert sequence to string
    let seq_str = std::str::from_utf8(&sequence)
        .map_err(|e| format!("Invalid UTF-8 sequence: {}", e))?;
    
    // MAPQ (mapping quality)
    let mapq = entry.mapping_quality.min(255);
    
    // Calculate TLEN (observed template length)
    let tlen = entry.target_end as i64 - entry.target_start as i64;
    
    // Build optional fields
    let mut optional_fields = Vec::new();
    
    // Add chain info as optional field
    optional_fields.push(format!("ch:Z:{}.{}.{}", 
        entry.chain_id, entry.chain_length, entry.chain_pos));
    
    // Add NM tag (edit distance)
    let (matches, mismatches, _, ins_bp, _, del_bp, _) = 
        calculate_cigar_stats(&entry.cigar)?;
    let edit_distance = mismatches + ins_bp + del_bp;
    optional_fields.push(format!("NM:i:{}", edit_distance));
    
    // Add AS tag (alignment score)
    // Simple scoring: +1 for match, -1 for mismatch/indel
    let alignment_score = matches as i64 - edit_distance as i64;
    optional_fields.push(format!("AS:i:{}", alignment_score));
    
    // Format SAM line
    // QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
    let sam_line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t{}\t{}\t*\t{}",
        entry.query_name,
        flag,
        ref_id,
        entry.target_start + 1, // SAM is 1-based
        mapq,
        sam_cigar,
        tlen,
        seq_str,
        optional_fields.join("\t")
    );
    
    Ok(sam_line)
}

// Write SAM content to a writer
fn write_sam_content<W: Write>(writer: &mut W, entries: &[PafEntry], query_db: &SequenceDB) -> Result<(), Box<dyn Error>> {
    // Generate and write header
    let header = generate_sam_header(entries);
    write!(writer, "{}", header)?;
    
    // Create reference ID map for quick lookups
    let mut reference_id_map = HashMap::new();
    let mut ref_id = 0;
    for entry in entries {
        if !reference_id_map.contains_key(&entry.target_name) {
            reference_id_map.insert(entry.target_name.clone(), ref_id);
            ref_id += 1;
        }
    }
    
    // Write alignment lines
    for entry in entries {
        match paf_entry_to_sam(entry, query_db, &reference_id_map) {
            Ok(sam_line) => writeln!(writer, "{}", sam_line)?,
            Err(e) => warn!("Failed to convert entry to SAM: {}", e),
        }
    }
    
    Ok(())
}

// Group entries by chain_id
fn group_by_chain(entries: Vec<PafEntry>) -> HashMap<u64, Vec<PafEntry>> {
    let mut chain_groups = HashMap::new();
    
    for entry in entries {
        chain_groups.entry(entry.chain_id).or_insert_with(Vec::new).push(entry);
    }
    
    chain_groups
}

// Helper function to update tags in a PafEntry
fn update_tags(tags: &[(String, String)], cigar: &str, chain_id: u64) -> Vec<(String, String)> {
    let mut new_tags = Vec::with_capacity(tags.len());
    for (tag, value) in tags {
        if tag == "cg:Z" {
            new_tags.push((tag.clone(), cigar.to_string()));
        } else if tag == "ch:Z" {
            new_tags.push((tag.clone(), format!("{}.1.1", chain_id)));
        } else {
            new_tags.push((tag.clone(), value.clone()));
        }
    }
    new_tags
}

// Helper function to create a merged PafEntry
fn create_merged_entry(
    base_entry: &PafEntry,
    next_entry: &PafEntry,
    merged_cigar: &str
) -> PafEntry {
    // Calculate statistics from merged CIGAR
    let (matches, _, _, _, _, _, block_len) = calculate_cigar_stats(merged_cigar)
        .expect("Failed to calculate CIGAR stats");
    
    // Adjust target range based on strand
    let (adjusted_target_start, adjusted_target_end) = if base_entry.strand == '-' {
        (next_entry.target_start, base_entry.target_end)
    } else {
        (base_entry.target_start, next_entry.target_end)
    };

    // Create the merged entry
    PafEntry {
        query_name: base_entry.query_name.clone(),
        query_length: base_entry.query_length,
        query_start: base_entry.query_start,
        query_end: next_entry.query_end,
        strand: base_entry.strand,
        target_name: base_entry.target_name.clone(),
        target_length: base_entry.target_length,
        target_start: adjusted_target_start,
        target_end: adjusted_target_end,
        num_matches: matches,
        alignment_length: block_len,
        mapping_quality: base_entry.mapping_quality,
        tags: update_tags(&base_entry.tags, merged_cigar, base_entry.chain_id),
        chain_id: base_entry.chain_id,
        chain_length: 1,
        chain_pos: 1,
        cigar: merged_cigar.to_string(),
    }
}

// Process a chain
fn process_chain(
    chain: &[PafEntry], 
    query_db: &SequenceDB, 
    target_db: &SequenceDB,
    erosion_size: usize
) -> Result<Vec<PafEntry>, Box<dyn Error>> {
    // Handle simple cases
    if chain.len() <= 1 {
        debug!("Chain has only one entry, returning unchanged");
        return Ok(chain.to_vec());
    }
        
    // Sort chain by position
    let mut sorted_chain = chain.to_vec();
    sorted_chain.sort_by_key(|entry| entry.chain_pos);
    
    // Start with the first entry as our base for merging
    let mut merged_entry = sorted_chain[0].clone();
    
    // Process each pair of adjacent entries
    for i in 0..(sorted_chain.len() - 1) {      
        debug!("\tProcessing entries {} and {}", i, i + 1);

        // Get the current base entry and the next entry to merge
        let current = &merged_entry;
        let next = &sorted_chain[i+1];

        debug!("\tMerged query {}-{}; {}; target {}-{}", merged_entry.query_start, merged_entry.query_end, merged_entry.strand, merged_entry.target_start, merged_entry.target_end);
        debug!("\tNext   query {}-{}; {}; target {}-{}", next.query_start, next.query_end, next.strand, next.target_start, next.target_end);
        //debug!("\tCIGAR: {}", merged_entry.cigar);
        debug!("\tChain: {}.{}.{}", merged_entry.chain_id, merged_entry.chain_length, next.chain_pos);

        // Check for query overlap
        let query_overlap = if next.query_start < current.query_end {
            // Calculate actual overlap size in base pairs
            current.query_end.saturating_sub(next.query_start) as usize
        } else {
            0
        };

        // Check for target overlap with inverted logic for reverse strand
        let target_overlap = if current.strand == '-' {
            if next.target_end > current.target_start {
                // Calculate actual overlap size in base pairs
                next.target_end.saturating_sub(current.target_start) as usize
            } else {
                0
            }
        } else {
            if next.target_start < current.target_end {
                // Calculate actual overlap size in base pairs
                current.target_end.saturating_sub(next.target_start) as usize
            } else {
                0
            }
        };

        // Determine minimum erosion sizes for query and target separately
        let min_query_erosion_size = if query_overlap > 0 {
            // Use at least the overlap size for erosion if an overlap exists
            let min_size = std::cmp::max(erosion_size, query_overlap);
            debug!("\tUsing min. query erosion size of {} bp due to query overlap of {} bp", 
                  min_size, query_overlap);
            min_size
        } else {
            // No overlap, use the specified erosion size
            debug!("\tUsing specified erosion size of {} bp due to no query overlap detected", erosion_size);
            erosion_size
        };
        
        let min_target_erosion_size = if target_overlap > 0 {
            // Use at least the overlap size for erosion if an overlap exists
            let min_size = std::cmp::max(erosion_size, target_overlap);
            debug!("\tUsing min. target erosion size of {} bp due to target overlap of {} bp", 
                  min_size, target_overlap);
            min_size
        } else {
            // No overlap, use the specified erosion size
            debug!("\tUsing specified erosion size of {} bp due to no target overlap detected", erosion_size);
            erosion_size
        };
        
        // Erode CIGARs at boundaries with the effective erosion sizes
        let (eroded_current_cigar, current_query_removed, current_target_removed) = 
        if min_query_erosion_size > 0 || min_target_erosion_size > 0 {
            if current.strand == '-' {
                // For reverse strand, erode the start of the current entry
                erode_cigar_start(&current.cigar, min_query_erosion_size, min_target_erosion_size)?
            } else {
                // For forward strand, erode the end of the current entry
                erode_cigar_end(&current.cigar, min_query_erosion_size, min_target_erosion_size)?
            }
        } else {
            (current.cigar.clone(), 0, 0)
        };
        debug!("\tEroded CIGARs - Current: query: {}, target: {}", current_query_removed, current_target_removed);

        let (eroded_next_cigar, next_query_removed, next_target_removed) = 
        // Only erode the start of next entry if there's no overlap
        // (otherwise we've already handled the overlap by eroding the current entry)
        if (min_query_erosion_size > 0 && query_overlap == 0) || 
        (min_target_erosion_size > 0 && target_overlap == 0) {
            if current.strand == '-' {
                // For reverse strand, erode the end of the next entry
                erode_cigar_end(&next.cigar, min_query_erosion_size, min_target_erosion_size)?
            } else {
                // For forward strand, erode the start of the next entry
                erode_cigar_start(&next.cigar, min_query_erosion_size, min_target_erosion_size)?
            }
        } else {
            (next.cigar.clone(), 0, 0)
        };
        debug!("\tEroded CIGARs - Next: query: {}, target: {}", next_query_removed, next_target_removed);

        // Calculate gap coordinates accounting for erosion
        let query_gap_start = current.query_end - current_query_removed as u64;
        let query_gap_end = next.query_start + next_query_removed as u64;
        // Check the strand to determine the target gap coordinates
        let (target_gap_start, target_gap_end) = if current.strand == '-' {
            // For reverse strand, subtract the erosion from the start of the next entry
            (next.target_end - next_target_removed as u64, current.target_start + current_target_removed as u64)
        } else {
            // For forward strand, add the erosion to the start of the next entry
            (current.target_end - current_target_removed as u64, next.target_start + next_target_removed as u64)
        };
        debug!("\tGap coordinates - Query: {}-{}, Target: {}-{}", query_gap_start, query_gap_end, target_gap_start, target_gap_end);

        // Calculate gap sizes
        let query_gap = query_gap_end - query_gap_start;
        let target_gap = target_gap_end - target_gap_start;
        debug!("\tAfter erosion - Query gap: {}, Target gap: {}", query_gap, target_gap);
        
        // Fetch sequences for alignment
        let gap_cigar = if query_gap == 0 && target_gap == 0 {
            // No gaps after erosion, nothing to align
            String::new()
        } else if query_gap > 0 && target_gap == 0 {
            // Query gap only - add insertion
            format!("{}I", query_gap)
        } else if query_gap == 0 && target_gap > 0 {
            // Target gap only - add deletion
            format!("{}D", target_gap)
        } else {
            // Both query and target have gaps - perform alignment
            let gap_query_seq = query_db.get_subsequence(
                &current.query_name,
                query_gap_start,
                query_gap_end,
                current.strand == '-'
            )?;
            
            let gap_target_seq = target_db.get_subsequence(
                &current.target_name,
                target_gap_start,
                target_gap_end,
                false
            )?;
            
            // Align the gap sequences
            align_sequences_wfa(&gap_query_seq, &gap_target_seq)?
        };
        
        // Merge CIGARs
        let merged_cigar = if current.strand == '-' {
            merge_cigars(&[&eroded_next_cigar, &gap_cigar, &eroded_current_cigar])?
        } else {
            merge_cigars(&[&eroded_current_cigar, &gap_cigar, &eroded_next_cigar])?
        };
        
        // Create the merged entry
        merged_entry = create_merged_entry(
            current, next, &merged_cigar
        );
    }
    
    debug!("\tMerged query {}-{}; {}; target {}-{}", merged_entry.query_start, merged_entry.query_end, merged_entry.strand, merged_entry.target_start, merged_entry.target_end);

    Ok(vec![merged_entry])
}


fn main() -> Result<(), Box<dyn Error>> {
    // Parse command-line arguments
    let args = Args::parse();
    
    // Initialize logging
    env_logger::Builder::new()
        .filter_level(match args.verbose {
            0 => log::LevelFilter::Warn,    // Errors and warnings
            1 => log::LevelFilter::Info,    // Errors, warnings, and info
            _ => log::LevelFilter::Debug,   // Errors, warnings, info, and debug
        })
        .init();
    
    // Set number of threads if needed
    if args.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()?;
    }
    
    // Load sequence databases
    info!("Loading query sequences from {}...", args.query.display());
    let query_db = SequenceDB::new(&args.query)?;
    
    info!("Loading target sequences from {}...", args.target.display());
    let target_db = SequenceDB::new(&args.target)?;
    
    // Read and process PAF file
    info!("Reading PAF file {}...", args.paf.display());
    let entries = read_paf_file(&args.paf)?;
    info!("Read {} PAF entries", entries.len());
    
    // Group entries by chain
    let chain_groups = group_by_chain(entries);
    info!("Found {} chains", chain_groups.len());
    
    // Process each chain
    let mut processed_entries = Vec::new();
    
    for (chain_id, group) in chain_groups {
        debug!("Processing chain {} with {} entries", chain_id, group.len());
        
        // Process chain and add results
        let chain_result = process_chain(&group, &query_db, &target_db, args.erosion_size)?;
        processed_entries.extend(chain_result);
    }
    
    // Determine output path
    let mut writer = match args.output {
        Some(file_path) => {
            let file = File::create(file_path)?;
            Box::new(io::BufWriter::new(file)) as Box<dyn Write>
        },
        None => {
            Box::new(io::BufWriter::new(io::stdout())) as Box<dyn Write>
        }
    };

    // Write output in appropriate format
    if args.sam {
        write_sam_content(&mut writer, &processed_entries, &query_db)?;
    } else {
        write_paf_content(&mut writer, &processed_entries)?;
    }
    
    info!("Done!");
    Ok(())
}