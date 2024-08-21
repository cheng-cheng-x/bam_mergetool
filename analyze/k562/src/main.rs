use std::process::{Command, Stdio};
use std::fs::{File, metadata};  
use std::io::{self, Write,BufReader, BufRead};
use std::path::Path;
use std::fs;
use std::collections::HashMap;
use std::thread;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use anyhow::{Context, Result};
use reqwest::Client;
use tokio::fs as tokio_fs; 
use tokio::io::AsyncWriteExt;
use std::env;


async fn download_file(url: &str, destination: &str) -> Result<(), anyhow::Error> {
    let client = Client::new();
    let mut response = client.get(url).send().await?;

    let mut file = tokio_fs::File::create(destination).await?;
    
    while let Some(chunk) = response.chunk().await? {
        file.write_all(&chunk).await?;
    }
    
    Ok(())
}

async fn download_files_from_list(file_list_path: &str, bam_dir: &str) -> Result<(), anyhow::Error> {
    let file_content = tokio_fs::read_to_string(file_list_path).await?;
    let urls: Vec<String> = file_content.lines().map(String::from).collect();
    
    let mut tasks = vec![];
    for url in urls {
        let file_name = url.split('/').last().unwrap();
        let destination = format!("{}/{}", bam_dir, file_name);
        
        tasks.push(tokio::spawn(async move {
            if !tokio_fs::metadata(&destination).await.is_ok() {
                if let Err(e) = download_file(&url, &destination).await {
                    eprintln!("Error downloading {}: {:?}", url, e);
                }
            } else {
                println!("File already exists: {}", destination);
            }
        }));
    }
    
    for task in tasks {
        task.await.unwrap();
    }
    
    Ok(())
}
//遍历目录下的每个文件
fn get_bam_files(bam_dir: &str) -> Vec<String> {
    let mut bam_files = Vec::new();

    if let Ok(entries) = fs::read_dir(bam_dir) {
        for entry in entries {
            if let Ok(entry) = entry {
                let path = entry.path();
                if path.is_file() && path.extension().and_then(|ext| ext.to_str()) == Some("bam") {
                    if let Some(path_str) = path.to_str() {
                        bam_files.push(path_str.to_string());
                    }
                }
            }
        }
    }

    bam_files
}

//创建bam文件的索引

fn index_bam(bam_file: &str) -> io::Result<()> {
    let bai_file = format!("{}.bai", bam_file);

    if !Path::new(&bai_file).exists() {
        println!("Indexing BAM file: {}", bam_file);
        
        // 直接使用指定的 samtools 路径
        let samtools_path = "/data/haocheng/miniconda3/envs/rust/bin/samtools";

        // 检查 samtools 是否存在
        if !Path::new(samtools_path).exists() {
            eprintln!("samtools not found at the specified path.");
            return Err(io::Error::new(io::ErrorKind::NotFound, "samtools not found"));
        }

        let status = Command::new(samtools_path)
            .arg("index")
            .arg(bam_file)
            .status()?;

        if status.success() {
            println!("Index created: {}", bai_file);
        } else {
            eprintln!("Error in creating BAM index.");
        }
    } else {
        println!("Index already exists: {}", bai_file);
    }
    Ok(())
}

fn get_coverage(bam_file: &str, chromosomes: &[&str], temp_output_file: &str) -> io::Result<()> {
    // 执行 bedtools genomecov 命令以计算覆盖度，结果保存到临时文件
    let status = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "source /data/haocheng/miniconda3/bin/activate rust && bedtools genomecov -ibam {} -bg > {}",
            bam_file,
            temp_output_file
        ))
        .status()?;

    if status.success() {
        println!("Coverage data written to {}", temp_output_file);
    } else {
        eprintln!("Error in executing bedtools.");
    }

    Ok(())
}

fn filter_chromosome_data(temp_output_file: &str, chromosomes: &[&str], output_file: &str) -> io::Result<()> {
    let mut temp_file = File::open(temp_output_file)?;
    let reader = io::BufReader::new(temp_file);

    let mut output_file = File::create(output_file)?;
    let mut writer = io::BufWriter::new(output_file);

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        if chromosomes.iter().any(|&chr| parts[0] == chr) {
            if !parts[0].contains('_') {
                writeln!(writer, "{}", line)?;
            }
        }
    }

    Ok(())
}



//判断文件是否有值
fn is_file_empty(file_path: &str) -> io::Result<bool> {
    if Path::new(file_path).exists() {
        let metadata = std::fs::metadata(file_path)?;
        Ok(metadata.len() == 0)
    } else {
        Ok(true) // 如果文件不存在，也视为 "空"
    }
}

//构建结构体
struct ChromosomeRegion {
    starts: [Vec<i32>; 2],
    ends: [Vec<i32>; 2],
}

//读取标准染色体长度
fn read_chromosome_size(file_path: &str, chromosome: &str) -> io::Result<i32> {
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);
    
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 2 && parts[0] == chromosome {
            // 返回对应染色体的长度
            return parts[1].parse::<i32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e));
        }
    }
    
    // 如果没有找到指定的染色体
    Err(io::Error::new(io::ErrorKind::NotFound, format!("Chromosome {} not found", chromosome)))
}

fn create_chromosome_regions(chromosome_length: usize) -> Vec<ChromosomeRegion> {
    let mut regions = Vec::new();

    // 根据染色体长度创建相应数量的 ChromosomeRegion 对象
    for _ in 0..chromosome_length {
        // 创建具有两个空 Vec<i32> 的 ChromosomeRegion 对象
        let region = ChromosomeRegion { 
            starts: [Vec::new(), Vec::new()], 
            ends: [Vec::new(), Vec::new()] 
        };
        regions.push(region);
    }

    regions
}

//写入数据
fn build_chromosome_data(file_path: &str) -> io::Result<HashMap<String, Vec<(usize, usize, i32)>>> {
    let mut chromosome_data = HashMap::new();

    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        
        let chromosome = parts[0].to_string();
        let start = parts[1].parse::<usize>().unwrap();
        let end = parts[2].parse::<usize>().unwrap();
        let value = parts[3].parse::<i32>().unwrap();

            chromosome_data
                .entry(chromosome)
                .or_insert_with(Vec::new)
                .push((start, end, value));
        
    }

    Ok(chromosome_data)
}
fn assign_values_to_regions(
    chromosome_data: &HashMap<String, Vec<(usize, usize, i32)>>,
    regions: &mut Vec<ChromosomeRegion>,
    chromosome: &str,
    file_index: usize
) -> io::Result<()> {
    if let Some(data) = chromosome_data.get(chromosome) {
        for &(start, end, value) in data {
            let start_index = start - 1; // 转换为索引
            let end_index = end - 1; // 转换为索引

            regions[start_index].starts[file_index].push(value);
            regions[end_index].ends[file_index].push(value);
        }
    } else {
        eprintln!("Chromosome {} not found in data.", chromosome);
    }

    Ok(())
}



//寻找有值对象的数量
fn count_non_empty_regions(regions: &[ChromosomeRegion]) -> usize {
    regions.iter().filter(|region| 
        !region.starts[0].is_empty() || !region.starts[1].is_empty() || 
        !region.ends[0].is_empty() || !region.ends[1].is_empty()
    ).count()
}

//合并
fn generate_bedgraph(
    regions: &Vec<ChromosomeRegion>,
    chromosome: &str,
    output_bedgraph_path: &str,
    include_initial_line: bool, // 布尔参数
) -> io::Result<()> {
    let mut bedgraph_file = File::create(output_bedgraph_path)?;

    let mut prev_end = 1; // 初始化起始位置
    let mut current_value = 0;
    let mut lines: Vec<String> = Vec::new(); // 用于存储所有行

    // 根据 include_initial_line 的值决定是否初始化第一行
    if include_initial_line {
        lines.push(format!("{}\t{}\t{}\t{}", chromosome, 1, 1, 0));
    }

    for (i, region) in regions.iter().enumerate() {
        // 只要有一个 starts 或 ends 非空，即可开始处理
        if !region.starts[0].is_empty() || !region.starts[1].is_empty() || !region.ends[0].is_empty() || !region.ends[1].is_empty() {
            let end = i + 1; // 当前索引 +1 作为 end 的值

            // 计算 value，忽略空 Vec
            let starts_0_sum: i32 = region.starts[0].iter().sum();
            let starts_1_sum: i32 = region.starts[1].iter().sum();
            let ends_0_sum: i32 = region.ends[0].iter().sum();
            let ends_1_sum: i32 = region.ends[1].iter().sum();

            // 存储当前行
            lines.push(format!("{}\t{}\t{}\t{}", chromosome, prev_end, end, current_value));
            current_value += starts_0_sum + starts_1_sum - ends_0_sum - ends_1_sum;

            prev_end = end; // 更新 prev_end，为下一次循环做准备
        }
    }

    // 写入去除第一行后的内容
    for line in lines.iter().skip(1){
        writeln!(bedgraph_file, "{}", line)?;
    }

    Ok(())
}

type ChromosomeData = HashMap<String, Vec<(usize, usize, i32)>>;
fn process_chromosomes(
    chromosomes: &[&str],
    chromosome_file: &str,
    chromosome_data1: &ChromosomeData,
    chromosome_data2: &ChromosomeData,
    work_dir: &str,
    output_file1: &str,
    output_file2: &str,
    num_threads: usize, // 添加线程数作为参数
    some_flag: bool,
) -> Result<()> {
    // 设置 Rayon 使用的线程数量
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .context("Failed to build thread pool")?;

    // 使用 Rayon 处理染色体数据
    pool.install(|| {
        chromosomes.par_iter().try_for_each(|&chromosome| -> Result<()> {
            // 读取染色体长度
            let chromosome_length = read_chromosome_size(chromosome_file, chromosome)
                .context("Failed to read chromosome size")?;

            // 创建 ChromosomeRegion 对象的组合
            let mut regions = create_chromosome_regions(chromosome_length as usize);

            println!("Created {} chromosome regions for {}", regions.len(), chromosome);

            // 给 regions 赋值
            assign_values_to_regions(&chromosome_data1, &mut regions, chromosome, 0)
                .context("Failed to assign values from first file")?;
            println!("Assigned values from {} to regions", output_file1);
            
            assign_values_to_regions(&chromosome_data2, &mut regions, chromosome, 1)
                .context("Failed to assign values from second file")?;
            println!("Assigned values from {} to regions", output_file2);

            println!("Assigned values to {} chromosome regions for {}", regions.len(), chromosome);

            // 统计非空的 regions 数量
            let non_empty_count = count_non_empty_regions(&regions);
            println!("Number of non-empty regions for {}: {}", chromosome, non_empty_count);

            // 生成 bedgraph 文件
            let output_bedgraph_path = format!("{}/output_{}.bedgraph", work_dir, chromosome);
            generate_bedgraph(&regions, chromosome, &output_bedgraph_path, some_flag)
                .context("Failed to generate bedgraph file")?;

            Ok(())
        })
    })
}


fn merge_bedgraph_files(work_dir: &str, chromosomes: &[&str]) -> io::Result<HashMap<String, Vec<(usize, usize, i32)>>> {
    let mut chromosome_data = HashMap::new();

    for chromosome in chromosomes {
        let file_path = format!("{}/output_{}.bedgraph", work_dir, chromosome);
        
        // 打开文件并检查其是否存在
        if Path::new(&file_path).exists() {
            let file = File::open(&file_path)?;
            let reader = io::BufReader::new(file);

            // 读取文件内容并合并到 chromosome_data
            for line in reader.lines() {
                let line = line?;
                let parts: Vec<&str> = line.split('\t').collect();

                let chrom = parts[0].to_string();
                let start = parts[1].parse::<usize>().unwrap();
                let end = parts[2].parse::<usize>().unwrap();
                let value = parts[3].parse::<i32>().unwrap();

                chromosome_data
                    .entry(chrom)
                    .or_insert_with(Vec::new)
                    .push((start, end, value));
            }
        } else {
            println!("Warning: File {} does not exist", file_path);
        }
    }

    Ok(chromosome_data)
}


fn main() -> Result<()> {
    
    let rt = tokio::runtime::Runtime::new()?;
    
    rt.block_on(async {
        let file_list_path = "/data/haocheng/data/ACAT/txt/k562.txt";
        let bam_dir = "/data/haocheng/data/bam/k562";
        
        download_files_from_list(file_list_path, bam_dir).await
    })?;
    // 定义 BAM 目录和文件路径

    let bam_dir = "/data/haocheng/data/bam/k562";

    let work_dir = format!("{}/work", bam_dir);

    // 检查或创建 work 目录
    if !Path::new(&work_dir).exists() {
        std::fs::create_dir(&work_dir)?;
    }
    let mut bam_files = get_bam_files(bam_dir);
    
    // 创建 BAM 文件的索引
for bam_file in &bam_files {
    index_bam(bam_file)?;
}
    let first_bam = bam_files.remove(0);// 移除第一个元素
    let second_bam = bam_files.remove(0);
    



    let output_file1 = format!("{}/coverage1.bed", work_dir);
    let output_file2 = format!("{}/coverage2.bed", work_dir);

    let chromosomes = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ];

// 创建临时文件
let temp_output_file1 = format!("{}/temp_coverage1.bed", work_dir);
let temp_output_file2 = format!("{}/temp_coverage2.bed", work_dir);

// 克隆 
let first_bam_clone = first_bam.clone();
let second_bam_clone = second_bam.clone();
let temp_output_file1_clone = temp_output_file1.clone();
let temp_output_file2_clone = temp_output_file2.clone();
let output_file1_clone = output_file1.clone();
let output_file2_clone = output_file2.clone();
let chromosomes_clone = chromosomes.to_vec();
let chromosomes_clone1 = chromosomes.to_vec();
let chromosomes_clone2 = chromosomes.to_vec();



    if is_file_empty(&output_file1_clone)? {
        get_coverage(&first_bam_clone, &chromosomes_clone1, &temp_output_file1_clone)?;
        filter_chromosome_data(&temp_output_file1_clone, &chromosomes_clone1, &output_file1_clone)?;
    } else {
        println!("{} is not empty. Skipping coverage calculation.", output_file1_clone);
    }
    let chromosome_data1 = build_chromosome_data(&output_file1_clone)?;
    println!("Processed {} with data from {}", chromosome_data1.len(), output_file1_clone);

    if is_file_empty(&output_file2_clone)? {
        get_coverage(&second_bam_clone, &chromosomes_clone2, &temp_output_file2_clone)?;
        filter_chromosome_data(&temp_output_file2_clone, &chromosomes_clone2, &output_file2_clone)?;
    } else {
        println!("{} is not empty. Skipping coverage calculation.", output_file2_clone);
    }
    let chromosome_data2 = build_chromosome_data(&output_file2_clone)?;
    println!("Processed {} with data from {}", chromosome_data2.len(), output_file2_clone);



let chromosome_file = "/data/haocheng/data/hg38.chrom.sizes"; // 这里应替换为想要的染色体长度文件路径

let num_threads = 6; // 或根据需要设置线程数

// 调用封装好的函数
process_chromosomes(
    &chromosomes,
    chromosome_file,
    &chromosome_data1,
    &chromosome_data2,
    &work_dir,
    &output_file1,
    &output_file2,
    num_threads,
    true
)?;
// 删除文件
fs::remove_file(&output_file1)?;
fs::remove_file(&output_file2)?;
 

while !bam_files.is_empty() {
    let first_bam = bam_files.remove(0); // 移除第一个元素
    let output_file1 = format!("{}/coverage1.bed", work_dir);
    
    // 创建临时文件
    let temp_output_file1 = format!("{}/temp_coverage1.bed", work_dir);

    if is_file_empty(&output_file1).unwrap() {
        get_coverage(&first_bam, &chromosomes, &temp_output_file1).unwrap();
        filter_chromosome_data(&temp_output_file1, &chromosomes, &output_file1).unwrap();
    } else {
        println!("{} is not empty. Skipping coverage calculation.", output_file1);
    }

    let chromosome_data1 = build_chromosome_data(&output_file1).unwrap();
    println!("Processed {} with data from {}", chromosome_data1.len(), output_file1);

    // 合并所有 bedgraph 文件的数据
    let chromosome_data2 = merge_bedgraph_files(&work_dir, &chromosomes)?;
    println!("Processed {} with data", chromosome_data2.len());

    // 调用封装好的函数
    process_chromosomes(
        &chromosomes,
        chromosome_file,
        &chromosome_data1,
        &chromosome_data2,
        &work_dir,
        &output_file1,
        &output_file2,
        num_threads ,
        false
    )?;
    fs::remove_file(&output_file1)?;
}





    Ok(())

}
