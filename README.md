# pafchainer

A tool to process PAF (Pairwise Alignment Format) files and connect alignment chains using WFA (Wavefront Alignment Algorithm).

## Overview

`pafchainer` merges consecutive alignments in a chain by connecting them with wavefront alignments. It erodes the boundaries of adjacent alignments and fills in the gaps with new alignments, resulting in a single, contiguous alignment path. This tool is useful for creating complete end-to-end alignments from chain fragments.

## Features

- Processes PAF files containing chain information
- Uses WFA for fast and accurate gap alignments
- Configurable boundary erosion size
- Outputs in PAF or SAM format
- Multi-threaded processing for improved performance

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

## License

[MIT License](LICENSE)
