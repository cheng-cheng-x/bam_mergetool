// 更新 htslib 的导入路径
use rust_htslib::bam::{self, Reader, Read};

// 如果不需要 HashMap，可以删除这行
use std::collections::HashMap;

// 保留 File 的导入
use std::fs::File;

// 如果不需要 Write，可以删除 self 和 Write
use std::io::{self, Write};

// 如果不需要 Path，可以删除这行
use std::path::Path;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bam_file = "path/to/your.bam"; // 替换为实际的 BAM 文件路径
    let output_file = "output.bedgraph"; // 输出的 BEDGraph 文件路径

    // 打开 BAM 文件
    let mut reader = Reader::from_path(bam_file)?;
    
    // 存储覆盖度数据，key 为染色体，value 为对应的覆盖度向量
    let mut coverage_data: HashMap<String, Vec<u32>> = HashMap::new();

    // 遍历 BAM 文件中的每一条 alignment
    // 遍历 BAM 文件中的每一条 alignment
for record in reader.records() {
    let record = record?;
    
        let chrom = String::from_utf8_lossy(record.reference_name()).to_string();
        let start = record.start() as usize;
        let end = record.end() as usize;

        // 只处理 chr1 到 chrY 的染色体
        if chrom.starts_with("chr") {
            if let Ok(chrom_num) = chrom[3..].parse::<u32>() {
                if (1..=22).contains(&chrom_num) || chrom == "chrX" || chrom == "chrY" {
                    let chrom_length = coverage_data.entry(chrom.clone()).or_insert(vec![0; 0]);
                    if chrom_length.len() < end {
                        chrom_length.resize(end, 0); // 确保覆盖度向量足够长
                    }
                    add_coverage(chrom_length, start, end);
                }
            }
        
    }
}

    // 输出覆盖度信息
    let mut output = File::create(output_file)?;
    for (chrom, coverage) in coverage_data {
        for (pos, count) in coverage.iter().enumerate() {
            if *count > 0 {
                writeln!(output, "{}\t{}\t{}\t{}", chrom, pos, pos + 1, count)?;
            }
        }
    }

    Ok(())
}

fn add_coverage(coverage: &mut Vec<u32>, start: usize, end: usize) {
    for i in start..end {
        if i < coverage.len() {
            coverage[i] += 1;
        }
    }
}
