# bam_mergetool

This project is designed to merge multiple BAM files into a BED format file, summing up the signal values across corresponding regions. The tool automates the merging process and outputs a BED file, saving you from managing the space required to handle multiple BAM files.

## Getting Started

To use this tool, follow these steps:

### 1. Add Dependency

In your Rust `Cargo.toml` file, add the following dependency to pull the code directly from GitHub:

```toml
[dependencies]
bam_mergetool = { git = "https://github.com/cheng-cheng-x/bam_mergetool", branch = "master" }
```   
## 2. Use in Your Code

In your main function, include the following usage:

```rust
use bam_mergetool::bam_merge;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    bam_merge(
        input_name,
        Some(&file_list_path),
        &bam_dir,
        &work_dir,
        &chromosome_file,
        &chromosomes
    )?;
    Ok(())
}
```

### Parameters

- **input_name**: The name for the analysis, e.g., `k562` for the K562 cell line.
- **file_list_path**: A text file with download URLs for the BAM files, one per line. If you already have the BAM files downloaded, set this parameter to `None`.
- **bam_dir**: The directory where BAM files should be downloaded.
- **work_dir**: The directory for intermediate files and the final results.
- **chromosome_file**: A file with chromosome sizes. You can use `hg38.chrom.sizes` for GRCh38 as a reference.
- **chromosomes**: The chromosomes you want to filter.

## 3. Tool Requirements

Ensure that your environment includes `samtools` and `bedtools`. The program will automatically locate these tools. If you need to specify custom paths for these tools, you can modify their locations in the `bam_merge` function.
