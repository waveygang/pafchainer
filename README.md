# pafchainer

A tool for processing PAF (Pairwise Alignment Format) files and merging alignment chains using the Wavefront Alignment Algorithm (WFA).

## Overview

`pafchainer` merges consecutive alignments in a chain by connecting them with wavefront alignments. It erodes the boundaries of adjacent alignments and fills in the gaps with new alignments, creating a single, contiguous alignment path. This is particularly useful for genomic alignment workflows that require complete end-to-end alignments.

## Features

- Processes PAF files with chain information
- Uses the Wavefront Alignment Algorithm for fast and accurate gap alignments
- Indexes chains for efficient processing
- Handles both compressed (BGZF) and uncompressed PAF files
- Configurable boundary erosion size
- Outputs in PAF or SAM format
- Multi-threaded processing

## Installation

```bash
cargo install pafchainer
```

Or build from source:

```bash
git clone https://github.com/waveygang/pafchainer
cd pafchainer
cargo build --release
```

## Usage

```bash
pafchainer --paf input.paf --query query.fa --target target.fa [OPTIONS]
```

### Required Arguments:

- `--paf, -p`: Input PAF file with chain information
- `--query, -q`: Query FASTA file
- `--target, -t`: Target FASTA file

### Optional Arguments:

- `--erosion-size, -e`: Size of boundary erosion in base pairs (default: 100)
- `--output, -o`: Output file (default: stdout)
- `--sam`: Output in SAM format instead of PAF
- `--threads, -t`: Number of threads to use (default: 4)
- `--verbose, -v`: Verbosity level (0=error, 1=info, 2=debug)

## Examples

Basic usage:
```bash
pafchainer -p alignments.paf -q query.fa -t reference.fa -o connected.paf
```

Output in SAM format:
```bash
pafchainer -p alignments.paf -q query.fa -t reference.fa --sam -o aligned.sam
```

## How It Works

1. Builds or loads a chain index for efficient access to PAF entries
2. For each chain, erodes boundaries between adjacent entries
3. Performs wavefront alignment to connect the eroded regions
4. Merges CIGAR strings to create a single continuous alignment
5. Outputs the merged alignment in PAF or SAM format

## License

[MIT License](LICENSE)
