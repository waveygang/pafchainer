use clap::Parser;
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus};
use libc;
use log::{debug, info, warn};
use rust_htslib::faidx::Reader as FastaReader;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use pafchainer::chain_index::{ChainIndex, PafEntry};
use pafchainer::cigar::{
    calculate_cigar_stats, cigar_u8_to_ops, erode_cigar_end, erode_cigar_start, merge_cigar_ops,
    ops_to_cigar, CigarOp,
};

use rayon::prelude::*;
use std::sync::{Arc, Mutex};

#[derive(Parser, Debug)]
#[clap(
    author,
    version,
    about = "A tool for merging WFMASH's alignment chains using the WFA algorithm."
)]
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

    #[clap(
        long = "wfa-params",
        default_value = "5,8,2,24,1",
        help = "WFA alignment parameters: mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2"
    )]
    wfa_params: String,

    /// Size of boundary erosion in base pairs
    #[clap(short, long, default_value = "200")]
    erosion_size: usize,

    /// Output PAF file
    #[clap(short, long)]
    output: Option<PathBuf>,

    /// Output in SAM format instead of PAF
    #[clap(short = 'a', long)]
    sam: bool,

    /// Number of threads to use
    #[clap(long, default_value = "4")]
    threads: NonZeroUsize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[clap(short, long, default_value = "0")]
    verbose: u8,
}

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
    fn get_subsequence(
        &self,
        name: &str,
        start: u64,
        end: u64,
        reverse_complement: bool,
    ) -> Result<Vec<u8>, Box<dyn Error>> {
        let start_usize = start as usize;
        let end_usize = std::cmp::min(end as usize, std::usize::MAX);

        // Fetch sequence and properly handle memory
        let seq = match self.reader.fetch_seq(name, start_usize, end_usize - 1) {
            Ok(seq) => {
                let mut seq_vec = seq.to_vec();
                unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) }; // Free up memory to avoid memory leak (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                seq_vec
                    .iter_mut()
                    .for_each(|byte| *byte = byte.to_ascii_uppercase());
                seq_vec
            }
            Err(e) => return Err(format!("Failed to fetch sequence for {}: {}", name, e).into()),
        };

        if reverse_complement {
            // Reverse complement the sequence
            Ok(seq
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'G' => b'C',
                    b'C' => b'G',
                    x => x,
                })
                .collect())
        } else {
            Ok(seq)
        }
    }
}

// Align two sequences using WFA
fn align_sequences_wfa(
    query: &[u8],
    target: &[u8],
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
) -> Result<Vec<CigarOp>, Box<dyn Error>> {
    debug!(
        "Performing WFA alignment between sequences of lengths {} and {}",
        query.len(),
        target.len()
    );

    // Initialize WFA aligner with gap-affine penalties
    let aligner = AffineWavefronts::with_penalties_affine2p(
        0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2,
    );

    // Perform alignment (note that WFA expects target, query order)
    let status = aligner.align(target, query);

    // Check alignment status (can't use != since AlignmentStatus doesn't implement PartialEq)
    match status {
        AlignmentStatus::Completed => {
            debug!("WFA alignment completed successfully");
        }
        s => {
            warn!("Alignment failed with status: {:?}", s);
            return Err(format!("Alignment failed with status: {:?}", s).into());
        }
    }

    // Convert CIGAR to CigarOps
    let cigar_ops = cigar_u8_to_ops(aligner.cigar());
    Ok(cigar_ops)
}

// Convert PafEntry to SAM format
fn paf_entry_to_sam(
    entry: &PafEntry,
    query_db: &SequenceDB,
    reference_info: &HashMap<String, u64>,
) -> Result<String, Box<dyn Error>> {
    // Calculate flag
    let flag = if entry.strand == '-' { 16 } else { 0 };

    // Get reference ID
    let ref_id = match reference_info.get(&entry.target_name) {
        Some(_) => entry.target_name.clone(),
        None => "*".to_string(),
    };

    // Fetch sequence
    let sequence = query_db.get_subsequence(
        &entry.query_name,
        entry.query_start,
        entry.query_end,
        entry.strand == '-',
    )?;

    // Convert sequence to string
    let seq_str =
        std::str::from_utf8(&sequence).map_err(|e| format!("Invalid UTF-8 sequence: {}", e))?;

    // MAPQ (mapping quality)
    let mapq = entry.mapping_quality.min(255);

    // Calculate TLEN (observed template length)
    let tlen = entry.target_end as i64 - entry.target_start as i64;

    // Build optional fields
    let mut optional_fields = Vec::new();

    // Add chain info as optional field
    optional_fields.push(format!(
        "ch:Z:{}.{}.{}",
        entry.chain_id, entry.chain_length, entry.chain_pos
    ));

    let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, _) =
        calculate_cigar_stats(&entry.cigar_ops).expect("Failed to calculate CIGAR stats");

    // Add NM tag (edit distance)
    let edit_distance = mismatches + inserted_bp + deleted_bp;
    optional_fields.push(format!("NM:i:{}", edit_distance));

    // Format bi and gi fields without trailing zeros
    let gap_compressed_identity =
        (matches as f64) / (matches + mismatches + insertions + deletions) as f64;
    let block_identity = (matches as f64) / (matches + edit_distance) as f64;
    optional_fields.push(
        format!("gi:f:{:.6}", gap_compressed_identity)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string(),
    );
    optional_fields.push(
        format!("bi:f:{:.6}", block_identity)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string(),
    );

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
        ops_to_cigar(&entry.cigar_ops),
        tlen,
        seq_str,
        optional_fields.join("\t")
    );

    Ok(sam_line)
}

// Process a chain
fn process_chain(
    chain: &[PafEntry],
    query_db: &SequenceDB,
    target_db: &SequenceDB,
    erosion_size: usize,
    wfa_params: &(i32, i32, i32, i32, i32),
) -> Result<PafEntry, Box<dyn Error>> {
    // Handle simple cases
    if chain.len() <= 1 {
        debug!("Chain has only one entry, returning unchanged");
        return Ok(chain[0].clone());
    }

    // Sort chain by position
    let mut sorted_chain = chain.to_vec();
    sorted_chain.sort_by_key(|entry| entry.chain_pos);

    // Start with the first entry as our base for merging
    let mut merged_entry = sorted_chain[0].clone();

    // Unpack WFA parameters
    let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = *wfa_params;

    // Process each pair of adjacent entries
    for i in 0..(sorted_chain.len() - 1) {
        debug!("\tProcessing entries {} and {}", i, i + 1);

        // Get the next entry to merge
        let next = &sorted_chain[i + 1];

        debug!(
            "\tMerged query {}-{}; {}; target {}-{}",
            merged_entry.query_start,
            merged_entry.query_end,
            merged_entry.strand,
            merged_entry.target_start,
            merged_entry.target_end
        );
        debug!(
            "\tNext   query {}-{}; {}; target {}-{}",
            next.query_start, next.query_end, next.strand, next.target_start, next.target_end
        );
        debug!(
            "\tChain: {}.{}.{}",
            merged_entry.chain_id, merged_entry.chain_length, next.chain_pos
        );

        // Check for query overlap
        let query_overlap = if next.query_start < merged_entry.query_end {
            // Calculate actual overlap size in base pairs
            merged_entry.query_end.saturating_sub(next.query_start) as usize
        } else {
            0
        };

        // Check for target overlap with inverted logic for reverse strand
        let target_overlap = if merged_entry.strand == '-' {
            if next.target_end > merged_entry.target_start {
                // Calculate actual overlap size in base pairs
                next.target_end.saturating_sub(merged_entry.target_start) as usize
            } else {
                0
            }
        } else {
            if next.target_start < merged_entry.target_end {
                // Calculate actual overlap size in base pairs
                merged_entry.target_end.saturating_sub(next.target_start) as usize
            } else {
                0
            }
        };

        // Determine minimum erosion sizes for query and target separately
        let min_query_erosion_size = if query_overlap > 0 {
            // Use at least the overlap size for erosion if an overlap exists
            let min_size = std::cmp::max(erosion_size, query_overlap);
            debug!(
                "\tUsing min. query erosion size of {} bp due to query overlap of {} bp",
                min_size, query_overlap
            );
            min_size
        } else {
            // No overlap, use the specified erosion size
            debug!(
                "\tUsing specified erosion size of {} bp due to no query overlap detected",
                erosion_size
            );
            erosion_size
        };

        let min_target_erosion_size = if target_overlap > 0 {
            // Use at least the overlap size for erosion if an overlap exists
            let min_size = std::cmp::max(erosion_size, target_overlap);
            debug!(
                "\tUsing min. target erosion size of {} bp due to target overlap of {} bp",
                min_size, target_overlap
            );
            min_size
        } else {
            // No overlap, use the specified erosion size
            debug!(
                "\tUsing specified erosion size of {} bp due to no target overlap detected",
                erosion_size
            );
            erosion_size
        };

        // Erode CIGARs at boundaries with the effective erosion sizes
        let (eroded_current_cigar_ops, current_query_removed, current_target_removed) =
            if min_query_erosion_size > 0 || min_target_erosion_size > 0 {
                if merged_entry.strand == '-' {
                    // For reverse strand, erode the start of the current entry
                    erode_cigar_start(
                        &merged_entry.cigar_ops,
                        min_query_erosion_size,
                        min_target_erosion_size,
                    )?
                } else {
                    // For forward strand, erode the end of the current entry
                    erode_cigar_end(
                        &merged_entry.cigar_ops,
                        min_query_erosion_size,
                        min_target_erosion_size,
                    )?
                }
            } else {
                (merged_entry.cigar_ops.clone(), 0, 0)
            };
        debug!(
            "\tEroded CIGARs - Current: query: {}, target: {}",
            current_query_removed, current_target_removed
        );

        // Only erode the start of next entry if there's no overlap
        // (otherwise we've already handled the overlap by eroding the current entry)
        let (eroded_next_cigar_ops, next_query_removed, next_target_removed) =
            if (min_query_erosion_size > 0 && query_overlap == 0)
                || (min_target_erosion_size > 0 && target_overlap == 0)
            {
                if merged_entry.strand == '-' {
                    // For reverse strand, erode the end of the next entry
                    erode_cigar_end(
                        &next.cigar_ops,
                        min_query_erosion_size,
                        min_target_erosion_size,
                    )?
                } else {
                    // For forward strand, erode the start of the next entry
                    erode_cigar_start(
                        &next.cigar_ops,
                        min_query_erosion_size,
                        min_target_erosion_size,
                    )?
                }
            } else {
                (next.cigar_ops.clone(), 0, 0)
            };
        debug!(
            "\tEroded CIGARs - Next: query: {}, target: {}",
            next_query_removed, next_target_removed
        );

        // Calculate gap coordinates accounting for erosion
        let query_gap_start = merged_entry.query_end - current_query_removed as u64;
        let query_gap_end = next.query_start + next_query_removed as u64;

        // Check the strand to determine the target gap coordinates
        let (target_gap_start, target_gap_end) = if merged_entry.strand == '-' {
            // For reverse strand, subtract the erosion from the start of the next entry
            (
                next.target_end - next_target_removed as u64,
                merged_entry.target_start + current_target_removed as u64,
            )
        } else {
            // For forward strand, add the erosion to the start of the next entry
            (
                merged_entry.target_end - current_target_removed as u64,
                next.target_start + next_target_removed as u64,
            )
        };
        debug!(
            "\tGap coordinates - Query: {}-{}, Target: {}-{}",
            query_gap_start, query_gap_end, target_gap_start, target_gap_end
        );

        // Calculate gap sizes
        let query_gap = query_gap_end - query_gap_start;
        let target_gap = target_gap_end - target_gap_start;
        debug!(
            "\tAfter erosion - Query gap: {}, Target gap: {}",
            query_gap, target_gap
        );

        // Generate gap CIGAR operations
        let gap_ops = if query_gap == 0 && target_gap == 0 {
            // No gaps after erosion, nothing to align
            Vec::new()
        } else if query_gap > 0 && target_gap == 0 {
            // Query gap only - add insertion
            vec![CigarOp('I', query_gap as usize)]
        } else if query_gap == 0 && target_gap > 0 {
            // Target gap only - add deletion
            vec![CigarOp('D', target_gap as usize)]
        } else {
            // Both query and target have gaps - perform alignment
            let gap_query_seq = query_db.get_subsequence(
                &merged_entry.query_name,
                query_gap_start,
                query_gap_end,
                merged_entry.strand == '-',
            )?;

            let gap_target_seq = target_db.get_subsequence(
                &merged_entry.target_name,
                target_gap_start,
                target_gap_end,
                false,
            )?;

            // Align the gap sequences
            align_sequences_wfa(
                &gap_query_seq,
                &gap_target_seq,
                mismatch,
                gap_open1,
                gap_ext1,
                gap_open2,
                gap_ext2,
            )?
        };

        // Merge CIGARs
        merged_entry.cigar_ops = if merged_entry.strand == '-' {
            merge_cigar_ops(&[&eroded_next_cigar_ops, &gap_ops, &eroded_current_cigar_ops])?
        } else {
            merge_cigar_ops(&[&eroded_current_cigar_ops, &gap_ops, &eroded_next_cigar_ops])?
        };

        // Update only the necessary fields of the merged entry
        merged_entry.query_end = next.query_end;

        // Adjust target boundaries based on strand
        if merged_entry.strand == '-' {
            merged_entry.target_start = next.target_start;
        } else {
            merged_entry.target_end = next.target_end;
        }
    }

    // Now that all merging is complete, update the remaining fields based on the final merged CIGAR
    let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) =
        calculate_cigar_stats(&merged_entry.cigar_ops).expect("Failed to calculate CIGAR stats");

    // Update the merged entry with final statistics
    merged_entry.num_matches = matches;
    merged_entry.alignment_length = block_len;
    merged_entry.chain_length = 1;
    merged_entry.chain_pos = 1;

    // Update tags
    let mut new_tags = Vec::with_capacity(merged_entry.tags.len());

    let edit_distance = mismatches + inserted_bp + deleted_bp;
    let gap_compressed_identity =
        (matches as f64) / (matches + mismatches + insertions + deletions) as f64;
    let block_identity = (matches as f64) / (matches + edit_distance) as f64;

    // Format bi and gi fields without trailing zeros
    for (tag, value) in &merged_entry.tags {
        if tag == "cg:Z" {
            new_tags.push((tag.clone(), ops_to_cigar(&merged_entry.cigar_ops)));
        } else if tag == "ch:Z" {
            new_tags.push((tag.clone(), format!("{}.1.1", merged_entry.chain_id)));
        } else if tag == "gi:f" {
            new_tags.push((
                tag.clone(),
                format!("{:.6}", gap_compressed_identity)
                    .trim_end_matches('0')
                    .trim_end_matches('.')
                    .to_string(),
            ));
        } else if tag == "bi:f" {
            new_tags.push((
                tag.clone(),
                format!("{:.6}", block_identity)
                    .trim_end_matches('0')
                    .trim_end_matches('.')
                    .to_string(),
            ));
        } else if tag == "md:f" {
            // Drop it
        } else {
            new_tags.push((tag.clone(), value.clone()));
        }
    }
    merged_entry.tags = new_tags;

    debug!(
        "\tMerged query {}-{}; {}; target {}-{}",
        merged_entry.query_start,
        merged_entry.query_end,
        merged_entry.strand,
        merged_entry.target_start,
        merged_entry.target_end
    );

    Ok(merged_entry)
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse command-line arguments
    let args = Args::parse();

    // Initialize logging
    env_logger::Builder::new()
        .filter_level(match args.verbose {
            0 => log::LevelFilter::Warn,  // Errors and warnings
            1 => log::LevelFilter::Info,  // Errors, warnings, and info
            _ => log::LevelFilter::Debug, // Errors, warnings, info, and debug
        })
        .init();

    // Configure thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads.into())
        .build_global()?;

    // Parse WFA parameters
    let wfa_params = {
        let parts: Vec<i32> = args
            .wfa_params
            .split(',')
            .map(|s| s.trim().parse::<i32>())
            .collect::<Result<Vec<i32>, _>>()?;

        if parts.len() != 5 {
            return Err("WFA parameters must have exactly 5 values: mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2".into());
        }

        (parts[0], parts[1], parts[2], parts[3], parts[4])
    };

    // Check for or build chain index
    let index_path = PathBuf::from(format!("{}.cidx", args.paf.to_string_lossy()));
    let chain_index = if index_path.exists() {
        info!("Loading chain index from {}...", index_path.display());
        ChainIndex::load(&index_path)?
    } else {
        info!("Building chain index for {}...", args.paf.display());
        let index = ChainIndex::build(&args.paf)?;
        info!("Saving chain index to {}...", index_path.display());
        index.save(&index_path)?;
        index
    };
    info!("Found {} chains in PAF file", chain_index.num_chains());

    // Create a writer that is a Boxed trait object implementing Write + Send.
    let mut writer: Box<dyn Write + Send> = match args.output {
        Some(file_path) => {
            let file = File::create(file_path)?;
            Box::new(io::BufWriter::new(file))
        }
        None => Box::new(io::BufWriter::new(io::stdout())),
    };

    // For SAM output, we need to handle the header with reference sequences
    let mut reference_info: HashMap<String, u64> = HashMap::new();
    if args.sam {
        // First pass: collect reference information
        for chain_id in chain_index.get_chain_ids() {
            let chain_entries = chain_index.load_chain(chain_id)?;
            for entry in &chain_entries {
                reference_info.insert(entry.target_name.clone(), entry.target_length);
            }
        }

        // Write SAM header
        write!(writer, "@HD\tVN:1.6\tSO:unsorted\n")?;

        // Write reference sequences with lengths collected from PAF entries
        for (ref_name, ref_length) in &reference_info {
            write!(writer, "@SQ\tSN:{}\tLN:{}\n", ref_name, ref_length)?;
        }

        // Add program information
        write!(writer, "@PG\tID:pafchainer\tPN:pafchainer\tVN:0.1.0\n")?;
    }

    // Prepare shared writer
    let writer = Arc::new(Mutex::new(writer));

    // Process in parallel and write entries
    chain_index
        .get_chain_ids()
        .into_par_iter()
        .for_each(|chain_id| {
            // Each thread creates its own FASTA readers.
            let query_db =
                SequenceDB::new(&args.query).expect("Failed to open query FASTA in thread");
            let target_db =
                SequenceDB::new(&args.target).expect("Failed to open target FASTA in thread");

            // Load and process the chain
            let chain_entries = chain_index
                .load_chain(chain_id)
                .expect("Failed to load chain");
            debug!(
                "Loaded {} entries for chain {}",
                chain_entries.len(),
                chain_id
            );
            let merged_entry = process_chain(
                &chain_entries,
                &query_db,
                &target_db,
                args.erosion_size,
                &wfa_params,
            )
            .expect("Failed to process chain");

            // Convert to SAM or PAF format and write immediately
            if args.sam {
                let sam_line = paf_entry_to_sam(&merged_entry, &query_db, &reference_info)
                    .expect("Failed to convert to SAM");
                let mut writer_lock = writer.lock().expect("Failed to lock writer");
                writeln!(writer_lock, "{}", sam_line).expect("Failed to write output");
            } else {
                let paf_line = merged_entry.to_string();
                let mut writer_lock = writer.lock().expect("Failed to lock writer");
                writeln!(writer_lock, "{}", paf_line).expect("Failed to write output");
            }
        });

    info!("Done!");
    Ok(())
}
