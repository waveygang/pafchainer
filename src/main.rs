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
        
        let seq = match self.reader.fetch_seq_string(name, start_usize, end_usize - 1) {
            Ok(seq) => seq.into_bytes(),
            Err(e) => return Err(format!("Failed to fetch sequence for {}: {}", name, e).into()),
        };
        
        if reverse_complement {
            // Reverse complement the sequence
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
                if tag_name == "chain:i" {
                    let chain_parts: Vec<&str> = tag_value.split('.').collect();
                    if chain_parts.len() >= 3 {
                        entry.chain_id = chain_parts[0].parse()?;
                        entry.chain_length = chain_parts[1].parse()?;
                        entry.chain_pos = chain_parts[2].parse()?;
                    } else {
                        return Err(format!("Invalid chain:i tag format: {}", tag_value).into());
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
            return Err("Missing required 'chain:i' or 'cg:Z' tag".into());
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

// Convert CIGAR operations to string
fn cigar_to_string(ops: &[CigarOp]) -> String {
    ops.iter().map(|&CigarOp(op, count)| format!("{}{}", count, op)).collect()
}

    // Erode CIGAR from the end
fn erode_cigar_end(cigar: &str, bp_count: usize) -> Result<(String, usize, usize), Box<dyn Error>> {
    let ops = parse_cigar(cigar)?;
    let mut new_ops = Vec::new();
    let mut remaining_bp = bp_count;
    let mut query_bases_removed = 0;
    let mut target_bases_removed = 0;
    
    info!("Eroding {} base pairs from end of CIGAR", bp_count);
    
    // Copy operations from beginning until we need to erode
    let mut i = 0;
    while i < ops.len() && remaining_bp > 0 {
        // Process operations from the end
        let idx = ops.len() - 1 - i;
        let CigarOp(op, count) = ops[idx];
        
        info!("erode_cigar_end - Processing operation: {}{}", count, op);
        
        match op {
            '=' | 'X' | 'M' => {
                // Advances both sequences
                if count > remaining_bp {
                    // Partial erosion
                    new_ops.insert(0, CigarOp(op, count - remaining_bp));
                    query_bases_removed += remaining_bp;
                    target_bases_removed += remaining_bp;
                    remaining_bp = 0;
                } else {
                    // Complete erosion
                    query_bases_removed += count;
                    target_bases_removed += count;
                    remaining_bp -= count;
                }
            },
            'I' => {
                // Insertion only advances query
                query_bases_removed += count;
                // Skip this operation if still eroding
                if remaining_bp == 0 {
                    new_ops.insert(0, CigarOp(op, count));
                }
            },
            'D' => {
                // Deletion only advances target
                target_bases_removed += count;
                // Skip this operation if still eroding
                if remaining_bp == 0 {
                    new_ops.insert(0, CigarOp(op, count));
                }
            },
            _ => new_ops.insert(0, CigarOp(op, count)), // Keep other operations
        }
        
        i += 1;
    }
    
    // Add remaining operations
    for j in 0..(ops.len() - i) {
        new_ops.insert(j, ops[j]);
    }
    
    Ok((cigar_to_string(&new_ops), query_bases_removed, target_bases_removed))
}

    // Erode CIGAR from the start
fn erode_cigar_start(cigar: &str, bp_count: usize) -> Result<(String, usize, usize), Box<dyn Error>> {
    let ops = parse_cigar(cigar)?;
    let mut new_ops = Vec::new();
    let mut remaining_bp = bp_count;
    let mut query_bases_removed = 0;
    let mut target_bases_removed = 0;
    
    info!("Eroding {} base pairs from start of CIGAR", bp_count);
    
    // Skip or partially include first operations until eroded enough
    let mut i = 0;
    while i < ops.len() && remaining_bp > 0 {
        let CigarOp(op, count) = ops[i];
        
        info!("erode_cigar_start - Processing operation: {}{}", count, op);
        
        match op {
            '=' | 'X' | 'M' => {
                // Advances both sequences
                if count > remaining_bp {
                    // Partial erosion
                    new_ops.push(CigarOp(op, count - remaining_bp));
                    query_bases_removed += remaining_bp;
                    target_bases_removed += remaining_bp;
                    remaining_bp = 0;
                } else {
                    // Complete erosion
                    query_bases_removed += count;
                    target_bases_removed += count;
                    remaining_bp -= count;
                }
            },
            'I' => {
                // Insertion only advances query
                query_bases_removed += count;
                // Skip this operation if still eroding
                if remaining_bp == 0 {
                    new_ops.push(CigarOp(op, count));
                }
            },
            'D' => {
                // Deletion only advances target
                target_bases_removed += count;
                // Skip this operation if still eroding
                if remaining_bp == 0 {
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

// Align two sequences using WFA
fn align_sequences_wfa(query: &[u8], target: &[u8]) -> Result<String, Box<dyn Error>> {
    info!("Performing WFA alignment between sequences of lengths {} and {}", query.len(), target.len());
    
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
            info!("WFA alignment completed successfully");
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

// Compute coordinates based on CIGAR string
fn compute_coordinates(
    query_start: u64,
    target_start: u64,
    cigar: &str
) -> Result<(u64, u64), Box<dyn Error>> {
    let cigar_ops = parse_cigar(cigar)?;
    
    let mut query_advance = 0u64;
    let mut target_advance = 0u64;
    
    // Calculate advancement based on CIGAR
    for &CigarOp(op, len) in &cigar_ops {
        let len = len as u64;
        match op {
            '=' | 'X' | 'M' => {
                // Advances both sequences
                query_advance += len;
                target_advance += len;
            },
            'I' => {
                // Insertion advances only query
                query_advance += len;
            },
            'D' => {
                // Deletion advances only target
                target_advance += len;
            },
            _ => {} // Other operations don't advance
        }
    }
    
    // Compute end coordinates
    let query_end = query_start + query_advance;
    let target_end = target_start + target_advance;
    
    Ok((query_end, target_end))
}

// Calculate CIGAR statistics
fn calculate_cigar_stats(cigar: &str) -> Result<(u64, u64, u64, u64, u64, u64, u64), Box<dyn Error>> {
    let ops = parse_cigar(cigar)?;
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

// Group entries by chain_id
fn group_by_chain(entries: Vec<PafEntry>) -> HashMap<u64, Vec<PafEntry>> {
    let mut chain_groups = HashMap::new();
    
    for entry in entries {
        chain_groups.entry(entry.chain_id).or_insert_with(Vec::new).push(entry);
    }
    
    chain_groups
}

// Process a chain
fn process_chain(
    chain: &[PafEntry], 
    query_db: &SequenceDB, 
    target_db: &SequenceDB,
    erosion_size: usize
) -> Result<Vec<PafEntry>, Box<dyn Error>> {
    if chain.len() <= 1 {
        // If only one entry, return it unchanged
        info!("Chain has only one entry, returning unchanged");
        return Ok(chain.to_vec());
    }
    
    info!("Processing chain with {} entries", chain.len());
    
    // Sort chain by position
    let mut sorted_chain = chain.to_vec();
    sorted_chain.sort_by_key(|entry| entry.chain_pos);
    
    // Start with the first entry
    let mut merged_entry = sorted_chain[0].clone();
    
    // Process each pair of adjacent entries
    for i in 0..(sorted_chain.len() - 1) {
        info!("Processing chain pair {}/{}", i+1, sorted_chain.len()-1);
        
        // Get the current merged entry and the next entry to merge
        let current = if i == 0 { sorted_chain[0].clone() } else { merged_entry };
        let next = &sorted_chain[i+1];
        
        // Erode CIGARs at boundaries only if erosion_size > 0
        let (eroded_current_cigar, current_query_removed, current_target_removed) = 
            if erosion_size > 0 {
                erode_cigar_end(&current.cigar, erosion_size)?
            } else {
                (current.cigar.clone(), 0, 0)
            };
        
        let (eroded_next_cigar, next_query_removed, next_target_removed) = 
            if erosion_size > 0 {
                erode_cigar_start(&next.cigar, erosion_size)?
            } else {
                (next.cigar.clone(), 0, 0)
            };
        
        // Check for gaps and calculate gap coordinates
        let query_gap_start = if erosion_size > 0 {
            current.query_end - current_query_removed as u64
        } else {
            current.query_end
        };
        
        let query_gap_end = if erosion_size > 0 {
            next.query_start + next_query_removed as u64
        } else {
            next.query_start
        };
        
        let target_gap_start = if erosion_size > 0 {
            current.target_end - current_target_removed as u64
        } else {
            current.target_end
        };
        
        let target_gap_end = if erosion_size > 0 {
            next.target_start + next_target_removed as u64
        } else {
            next.target_start
        };
        
        // Calculate gap sizes
        let query_gap = if query_gap_end > query_gap_start {
            query_gap_end - query_gap_start
        } else { 0 };
        
        let target_gap = if target_gap_end > target_gap_start {
            target_gap_end - target_gap_start
        } else { 0 };
        
        info!("Query gap: {}, Target gap: {}", query_gap, target_gap);
        
        // Special handling for no gaps
        if query_gap == 0 && target_gap == 0 && erosion_size == 0 {
            // If no gaps and no erosion, just merge the CIGAR strings directly
            info!("No erosion and no gaps, directly merging CIGAR strings");
            let merged_cigar = merge_cigars(&[&current.cigar, &next.cigar])?;
            
            // Calculate statistics from merged CIGAR
            let (matches, _mismatches, _insertions, _inserted_bp, _deletions, _deleted_bp, block_len) = 
                calculate_cigar_stats(&merged_cigar)?;
            
            // Create a map for new tag values
            let mut new_tag_values = HashMap::new();
            for (tag, _) in &current.tags {
                if matches!(tag.as_str(), "cg:Z" | "chain:i") {
                    // For tags we know we'll update, store their new values
                    if tag == "cg:Z" {
                        new_tag_values.insert(tag.clone(), merged_cigar.clone());
                    } else if tag == "chain:i" {
                        let chain_tag_value = format!("{}.1.1", current.chain_id);
                        new_tag_values.insert(tag.clone(), chain_tag_value);
                    }
                }
            }
            
            // Create new tags array preserving original order
            let mut new_tags = Vec::with_capacity(current.tags.len());
            for (tag, value) in &current.tags {
                if let Some(new_value) = new_tag_values.get(tag) {
                    // Use the new value for tags that need updating
                    new_tags.push((tag.clone(), new_value.clone()));
                } else {
                    // Keep the original value for other tags
                    new_tags.push((tag.clone(), value.clone()));
                }
            }
            
            // Create merged entry with updated CIGAR and end coordinates
            merged_entry = PafEntry {
                query_name: current.query_name,
                query_length: current.query_length,
                query_start: current.query_start,
                query_end: next.query_end,
                strand: current.strand,
                target_name: current.target_name,
                target_length: current.target_length,
                target_start: current.target_start,
                target_end: next.target_end,
                num_matches: matches,
                alignment_length: block_len,
                mapping_quality: current.mapping_quality,
                tags: new_tags,
                chain_id: current.chain_id,
                chain_length: 1,
                chain_pos: 1,
                cigar: merged_cigar,
            };
            
            continue;
        } else if query_gap > 0 && target_gap == 0 && erosion_size == 0 {
            // If gap in query only and no erosion, add insertions to the CIGAR
            info!("No erosion with query gap only, adding insertions to CIGAR");
            let insertion_cigar = format!("{}I", query_gap);
            let merged_cigar = merge_cigars(&[&current.cigar, &insertion_cigar, &next.cigar])?;
            
            // Calculate statistics from merged CIGAR
            let (matches, _mismatches, _insertions, _inserted_bp, _deletions, _deleted_bp, block_len) = 
                calculate_cigar_stats(&merged_cigar)?;
            
            // Create merged entry with insertion CIGAR
            let mut new_tag_values = HashMap::new();
            for (tag, _) in &current.tags {
                if matches!(tag.as_str(), "cg:Z" | "chain:i") {
                    // For tags we know we'll update, store their new values
                    if tag == "cg:Z" {
                        new_tag_values.insert(tag.clone(), merged_cigar.clone());
                    } else if tag == "chain:i" {
                        let chain_tag_value = format!("{}.1.1", current.chain_id);
                        new_tag_values.insert(tag.clone(), chain_tag_value);
                    }
                }
            }
            
            // Create new tags array preserving original order
            let mut new_tags = Vec::with_capacity(current.tags.len());
            for (tag, value) in &current.tags {
                if let Some(new_value) = new_tag_values.get(tag) {
                    // Use the new value for tags that need updating
                    new_tags.push((tag.clone(), new_value.clone()));
                } else {
                    // Keep the original value for other tags
                    new_tags.push((tag.clone(), value.clone()));
                }
            }
            
            merged_entry = PafEntry {
                query_name: current.query_name,
                query_length: current.query_length,
                query_start: current.query_start,
                query_end: next.query_end,
                strand: current.strand,
                target_name: current.target_name,
                target_length: current.target_length,
                target_start: current.target_start,
                target_end: next.target_end,
                num_matches: matches,
                alignment_length: block_len,
                mapping_quality: current.mapping_quality,
                tags: new_tags,
                chain_id: current.chain_id,
                chain_length: 1,
                chain_pos: 1,
                cigar: merged_cigar,
            };
            
            continue;
        } else if query_gap == 0 && target_gap > 0 && erosion_size == 0 {
            // If gap in target only and no erosion, add deletions to the CIGAR
            info!("No erosion with target gap only, adding deletions to CIGAR");
            let deletion_cigar = format!("{}D", target_gap);
            let merged_cigar = merge_cigars(&[&current.cigar, &deletion_cigar, &next.cigar])?;
            
            // Calculate statistics from merged CIGAR
            let (matches, _mismatches, _insertions, _inserted_bp, _deletions, _deleted_bp, block_len) = 
                calculate_cigar_stats(&merged_cigar)?;
            
            // Create merged entry with deletion CIGAR
            let mut new_tag_values = HashMap::new();
            for (tag, _) in &current.tags {
                if matches!(tag.as_str(), "cg:Z" | "chain:i") {
                    // For tags we know we'll update, store their new values
                    if tag == "cg:Z" {
                        new_tag_values.insert(tag.clone(), merged_cigar.clone());
                    } else if tag == "chain:i" {
                        let chain_tag_value = format!("{}.1.1", current.chain_id);
                        new_tag_values.insert(tag.clone(), chain_tag_value);
                    }
                }
            }
            
            // Create new tags array preserving original order
            let mut new_tags = Vec::with_capacity(current.tags.len());
            for (tag, value) in &current.tags {
                if let Some(new_value) = new_tag_values.get(tag) {
                    // Use the new value for tags that need updating
                    new_tags.push((tag.clone(), new_value.clone()));
                } else {
                    // Keep the original value for other tags
                    new_tags.push((tag.clone(), value.clone()));
                }
            }
            
            merged_entry = PafEntry {
                query_name: current.query_name,
                query_length: current.query_length,
                query_start: current.query_start,
                query_end: next.query_end,
                strand: current.strand,
                target_name: current.target_name,
                target_length: current.target_length,
                target_start: current.target_start,
                target_end: next.target_end,
                num_matches: matches,
                alignment_length: block_len,
                mapping_quality: current.mapping_quality,
                tags: new_tags,
                chain_id: current.chain_id,
                chain_length: 1,
                chain_pos: 1,
                cigar: merged_cigar,
            };
            
            continue;
        }
        
        // Handle gap alignment with or without erosion
        let mut gap_query_seq = Vec::new();
        let mut gap_target_seq = Vec::new();
        
        // Get gap sequences directly from computed gap positions
        if query_gap > 0 {
            info!("Getting gap query sequence of length {}", query_gap);
            gap_query_seq = query_db.get_subsequence(
                &current.query_name,
                query_gap_start,
                query_gap_end,
                current.strand == '-'
            )?;
            info!("Fetched sequence for {} from {} to {}", current.query_name, query_gap_start, query_gap_end);
            info!("Got gap query sequence of length {}", gap_query_seq.len());
        }
        
        if target_gap > 0 {
            info!("Getting gap target sequence of length {}", target_gap);
            gap_target_seq = target_db.get_subsequence(
                &current.target_name,
                target_gap_start,
                target_gap_end,
                false
            )?;
            info!("Fetched sequence for {} from {} to {}", current.target_name, target_gap_start, target_gap_end);
            info!("Got gap target sequence of length {}", gap_target_seq.len());
        }
        
        info!("Combined sequence lengths - Query: {}, Target: {}", gap_query_seq.len(), gap_target_seq.len());
        
        // Perform WFA alignment on gap sequences
        let gap_cigar = align_sequences_wfa(&gap_query_seq, &gap_target_seq)?;
        
        // Merge CIGARs
        let merged_cigar = merge_cigars(&[&eroded_current_cigar, &gap_cigar, &eroded_next_cigar])?;
        
        // Calculate statistics from merged CIGAR
        let (matches, _mismatches, _insertions, _inserted_bp, _deletions, _deleted_bp, block_len) = 
            calculate_cigar_stats(&merged_cigar)?;
        
        // Create a map for new tag values
        let mut new_tag_values = HashMap::new();
        for (tag, _) in &current.tags {
            if matches!(tag.as_str(), "cg:Z" | "chain:i") {
                // For tags we know we'll update, store their new values
                if tag == "cg:Z" {
                    new_tag_values.insert(tag.clone(), merged_cigar.clone());
                } else if tag == "chain:i" {
                    let chain_tag_value = format!("{}.1.1", current.chain_id);
                    new_tag_values.insert(tag.clone(), chain_tag_value);
                }
            }
        }
        
        // Create new tags array preserving original order
        let mut new_tags = Vec::with_capacity(current.tags.len());
        for (tag, value) in &current.tags {
            if let Some(new_value) = new_tag_values.get(tag) {
                // Use the new value for tags that need updating
                new_tags.push((tag.clone(), new_value.clone()));
            } else {
                // Keep the original value for other tags
                new_tags.push((tag.clone(), value.clone()));
            }
        }
        
        // Use the expected coordinates directly instead of computing them
        merged_entry = PafEntry {
            query_name: current.query_name,
            query_length: current.query_length,
            query_start: current.query_start,
            query_end: next.query_end,  // Use the end position from the next entry directly
            strand: current.strand,
            target_name: current.target_name,
            target_length: current.target_length,
            target_start: current.target_start,
            target_end: next.target_end,  // Use the end position from the next entry directly
            num_matches: matches,
            alignment_length: block_len,
            mapping_quality: 255,
            tags: new_tags,
            chain_id: current.chain_id,
            chain_length: 1,
            chain_pos: 1,
            cigar: merged_cigar,
        };
    }
    
    Ok(vec![merged_entry])
}

// Write PAF entries to file
fn write_paf_file<P: AsRef<Path>>(path: P, entries: &[PafEntry]) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = io::BufWriter::new(file);
    
    for entry in entries {
        writeln!(writer, "{}", entry.to_string())?;
    }
    
    Ok(())
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
        info!("Processing chain {} with {} entries", chain_id, group.len());
        
        // Verify CIGAR strings before processing
        info!("Verifying chain {}", chain_id);
        
        // Process chain and add results
        let chain_result = process_chain(&group, &query_db, &target_db, args.erosion_size)?;
        processed_entries.extend(chain_result);
    }
    
    // Write output
    let output_path = args.output.unwrap_or_else(|| PathBuf::from("output.paf"));
    info!("Writing {} processed entries to {}", processed_entries.len(), output_path.display());
    write_paf_file(&output_path, &processed_entries)?;
    
    info!("Done!");
    Ok(())
}