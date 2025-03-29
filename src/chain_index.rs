use log::{debug, info};
use noodles::bgzf;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

// PAF entry representation
#[derive(Debug, Clone)]
pub struct PafEntry {
    pub query_name: String,
    pub query_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: char,
    pub target_name: String,
    pub target_length: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub num_matches: u64,
    pub alignment_length: u64,
    pub mapping_quality: u64,
    pub tags: Vec<(String, String)>,
    pub chain_id: u64,
    pub chain_length: u64,
    pub chain_pos: u64,
    pub cigar: String,
}

impl PafEntry {
    // Parse a PAF entry from a line
    pub fn parse_from_line(line: &str) -> Result<Self, Box<dyn Error>> {
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
    pub fn to_string(&self) -> String {
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

/// Represents a single entry within a chain
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ChainEntryInfo {
    /// File offset (uncompressed for BGZF files)
    pub offset: u64,
    /// Length of the entry in bytes (including newline)
    pub length: usize,
    /// Position in the chain
    pub chain_pos: u64,
}

/// Index structure for efficient processing of PAF chains
#[derive(Serialize, Deserialize, Debug)]
pub struct ChainIndex {
    /// Maps chain_id to a vector of file offsets
    pub chains: HashMap<u64, Vec<ChainEntryInfo>>,
    /// Path to the PAF file
    pub paf_file: String,
    /// Whether the PAF file is compressed
    pub is_compressed: bool,
}

impl ChainIndex {
    /// Build an index from a PAF file
    pub fn build<P: AsRef<Path>>(paf_file: P) -> Result<Self, Box<dyn Error>> {
        let paf_path = paf_file.as_ref();
        let paf_str = paf_path.to_string_lossy().to_string();
        let is_compressed = [".gz", ".bgz"].iter().any(|e| paf_str.ends_with(e));

        info!("Building chain index for {}...", paf_str);
        let mut chains = HashMap::new();

        if is_compressed {
            // Handle compressed files
            Self::build_index_for_compressed_file(paf_path, &mut chains)?;
        } else {
            // Handle uncompressed files
            Self::build_index_for_uncompressed_file(paf_path, &mut chains)?;
        }

        // Sort entries by chain position
        for entries in chains.values_mut() {
            entries.sort_by_key(|e| e.chain_pos);
        }

        info!("Found {} chains", chains.len());
        Ok(Self {
            chains,
            paf_file: paf_str,
            is_compressed,
        })
    }

    /// Build index for a compressed PAF file
    fn build_index_for_compressed_file(
        paf_path: &Path,
        chains: &mut HashMap<u64, Vec<ChainEntryInfo>>,
    ) -> Result<(), Box<dyn Error>> {
        let file = File::open(paf_path)?;
        let reader = bgzf::Reader::new(file);
        let mut buf_reader = BufReader::new(reader);

        let mut offset = 0;
        let mut line = String::new();

        loop {
            line.clear();
            let bytes_read = buf_reader.read_line(&mut line)?;
            if bytes_read == 0 {
                break;
            }

            if let Ok(entry) = PafEntry::parse_from_line(&line) {
                let chain_entry = ChainEntryInfo {
                    offset,
                    length: bytes_read,
                    chain_pos: entry.chain_pos,
                };
                chains
                    .entry(entry.chain_id)
                    .or_insert_with(Vec::new)
                    .push(chain_entry);
            }

            offset += bytes_read as u64;
        }

        Ok(())
    }

    /// Build index for an uncompressed PAF file
    fn build_index_for_uncompressed_file(
        paf_path: &Path,
        chains: &mut HashMap<u64, Vec<ChainEntryInfo>>,
    ) -> Result<(), Box<dyn Error>> {
        let file = File::open(paf_path)?;
        let reader = BufReader::new(file);

        let mut offset = 0;

        for line_result in reader.lines() {
            let line = line_result?;
            let line_len = line.len() + 1; // +1 for newline

            if let Ok(entry) = PafEntry::parse_from_line(&line) {
                let chain_entry = ChainEntryInfo {
                    offset,
                    length: line_len,
                    chain_pos: entry.chain_pos,
                };
                chains
                    .entry(entry.chain_id)
                    .or_insert_with(Vec::new)
                    .push(chain_entry);
            }

            offset += line_len as u64;
        }

        Ok(())
    }

    /// Load all entries for a specific chain
    pub fn load_chain(&self, chain_id: u64) -> Result<Vec<PafEntry>, Box<dyn Error>> {
        let entries = match self.chains.get(&chain_id) {
            Some(entries) => entries,
            None => return Ok(Vec::new()),
        };

        debug!("Loading chain {} with {} entries", chain_id, entries.len());

        if self.is_compressed {
            self.load_chain_compressed(chain_id, entries)
        } else {
            self.load_chain_uncompressed(chain_id, entries)
        }
    }

    /// Load chain entries from a compressed file
    fn load_chain_compressed(
        &self,
        chain_id: u64,
        entries: &[ChainEntryInfo],
    ) -> Result<Vec<PafEntry>, Box<dyn Error>> {
        let gzi_path = format!("{}.gzi", self.paf_file);
        let gzi_index = bgzf::gzi::read(gzi_path)?;

        let file = File::open(&self.paf_file)?;
        let mut reader = bgzf::Reader::new(file);
        let mut paf_entries = Vec::with_capacity(entries.len());

        for entry_info in entries {
            let mut buffer = vec![0; entry_info.length];

            reader.seek_by_uncompressed_position(&gzi_index, entry_info.offset)?;
            reader.read_exact(&mut buffer)?;

            let line = std::str::from_utf8(&buffer[..buffer.len() - 1])?; // Remove newline
            if let Ok(paf_entry) = PafEntry::parse_from_line(line) {
                // Verify chain_id matches
                if paf_entry.chain_id == chain_id {
                    paf_entries.push(paf_entry);
                } else {
                    return Err(format!(
                        "Chain ID mismatch: expected {}, got {}",
                        chain_id, paf_entry.chain_id
                    )
                    .into());
                }
            }
        }

        Ok(paf_entries)
    }

    /// Load chain entries from an uncompressed file
    fn load_chain_uncompressed(
        &self,
        chain_id: u64,
        entries: &[ChainEntryInfo],
    ) -> Result<Vec<PafEntry>, Box<dyn Error>> {
        let file = File::open(&self.paf_file)?;
        let mut reader = BufReader::new(file);
        let mut paf_entries = Vec::with_capacity(entries.len());

        for entry_info in entries {
            let mut buffer = vec![0; entry_info.length];

            reader.seek(SeekFrom::Start(entry_info.offset))?;
            reader.read_exact(&mut buffer)?;

            let line = std::str::from_utf8(&buffer[..buffer.len() - 1])?; // Remove newline
            if let Ok(paf_entry) = PafEntry::parse_from_line(line) {
                // Verify chain_id matches
                if paf_entry.chain_id == chain_id {
                    paf_entries.push(paf_entry);
                } else {
                    return Err(format!(
                        "Chain ID mismatch: expected {}, got {}",
                        chain_id, paf_entry.chain_id
                    )
                    .into());
                }
            }
        }

        Ok(paf_entries)
    }

    /// Save the index to a file
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn Error>> {
        let file = File::create(path)?;
        let writer = io::BufWriter::new(file);
        bincode::serialize_into(writer, self)?;
        Ok(())
    }

    /// Load an index from a file
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let file = File::open(path)?;
        let reader = io::BufReader::new(file);
        let index: ChainIndex = bincode::deserialize_from(reader)?;
        Ok(index)
    }

    /// Get a list of all chain IDs in the index
    pub fn get_chain_ids(&self) -> Vec<u64> {
        self.chains.keys().cloned().collect()
    }

    /// Get the number of chains in the index
    pub fn num_chains(&self) -> usize {
        self.chains.len()
    }
}
