use std::process::{Command, Stdio};
use std::fs::{File, metadata};  
use std::io::{self, Write,BufReader, BufRead};
use std::path::Path;
use std::fs;
use std::collections::{HashMap,BTreeMap};
use std::thread;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use anyhow::{Context, Result};
use reqwest::Client;
use tokio::fs as tokio_fs; 
use tokio::io::AsyncWriteExt;
use std::env;
use std::error::Error;
use std::sync::{Arc, Mutex, RwLock,mpsc};




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
    starts: Vec<i32>,
    ends: Vec<i32>,
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
    let mut regions = Vec::with_capacity(chromosome_length);

    // 根据染色体长度创建相应数量的 ChromosomeRegion 对象
    for _ in 0..chromosome_length {
        // 创建具有两个空 Vec<i32> 的 ChromosomeRegion 对象
        let region = ChromosomeRegion {
            starts: Vec::new(),
            ends: Vec::new(),
        };
        regions.push(region);
    }

    regions
}

fn process_data_chunk(
    chunk: Vec<String>,
    chromosome_data: Arc<Mutex<HashMap<String, Vec<(usize, usize, i32)>>>>,
) {
    for line in chunk {
        let parts: Vec<&str> = line.split('\t').collect();
        let chromosome = parts[0].to_string();
        let start = parts[1].parse::<usize>().unwrap();
        let end = parts[2].parse::<usize>().unwrap();
        let value = parts[3].parse::<i32>().unwrap();

        chromosome_data
            .lock()
            .unwrap()
            .entry(chromosome)
            .or_insert_with(Vec::new)
            .push((start, end, value));
    }
}
fn process_file_in_chunks(
    file_path: &str,
    chunk_size: usize,
) -> Result<HashMap<String, Vec<(usize, usize, i32)>>> {
    let file = File::open(file_path)?;
    let mut reader = BufReader::new(file);
    let chromosome_data = Arc::new(Mutex::new(HashMap::new()));

    let pool = ThreadPoolBuilder::new().build().unwrap();

    loop {
        let mut buffer = Vec::with_capacity(chunk_size);
        let mut line = String::new();
        for _ in 0..chunk_size {
            if reader.read_line(&mut line)? == 0 {
                break; // End of file
            }
            buffer.push(line.trim().to_string());
            line.clear();
        }

        if buffer.is_empty() {
            break; // No more data
        }

        let data_clone = Arc::clone(&chromosome_data);
        pool.spawn(move || {
            process_data_chunk(buffer, data_clone);
        });
    }

    pool.install(|| {});

    Ok(Arc::try_unwrap(chromosome_data)
        .unwrap()
        .into_inner()
        .unwrap())
}

/* 
fn build_chromosome_data(file_path: &str) -> io::Result<HashMap<String, Vec<(usize, usize, i32)>>> {
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);

    // 创建一个线程安全的HashMap，使用Arc和RwLock包装
    let chromosome_data = Arc::new(RwLock::new(HashMap::new()));

    // 使用par_bridge将Iterator转换为ParallelIterator
    reader.lines().par_bridge().for_each(|line| {
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split('\t').collect();

            let chromosome = parts[0].to_string();
            let start = parts[1].parse::<usize>().unwrap();
            let end = parts[2].parse::<usize>().unwrap();
            let value = parts[3].parse::<i32>().unwrap();

            // 写入操作需要加锁
            let mut data = chromosome_data.write().unwrap();
            data.entry(chromosome)
                .or_insert_with(Vec::new)
                .push((start, end, value));
        }
    });

    // 将Arc解包并返回结果
    let final_data = Arc::try_unwrap(chromosome_data).unwrap().into_inner().unwrap();
    Ok(final_data)
}

*/
fn build_chromosome_data(
    file_path: &str,
    data: Arc<RwLock<HashMap<String, HashMap<usize, (i32, i32)>>>>
) -> io::Result<()> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    // 使用Rayon的默认线程池
    reader.lines()
        .par_bridge()
        .for_each(|line| {
            if let Ok(line) = line {
                let parts: Vec<&str> = line.split('\t').collect();
                let chromosome = parts[0].to_string();
                let start_pos = parts[1].parse::<usize>().unwrap();
                let end_pos = parts[2].parse::<usize>().unwrap();
                let value = parts[3].parse::<i32>().unwrap();

                // 使用写锁保护对data的可变访问
                let mut data = data.write().unwrap();
                let chromosome_entry = data.entry(chromosome.clone()).or_insert_with(HashMap::new);

                // 更新start_pos
                let (current_value_starts, current_value_ends) = chromosome_entry
                    .get(&start_pos)
                    .map_or((value, 0), |&(prev_value_starts, prev_value_ends)| (prev_value_starts + value, prev_value_ends));

                chromosome_entry.insert(start_pos, (current_value_starts, current_value_ends));

                // 更新end_pos
                let (current_value_starts, current_value_ends) = chromosome_entry
                    .get(&end_pos)
                    .map_or((0, value), |&(prev_value_starts, prev_value_ends)| (prev_value_starts, prev_value_ends + value));

                chromosome_entry.insert(end_pos, (current_value_starts, current_value_ends));
            }
        });

    Ok(())
}

fn write_chromosome_data_to_file(
    file_path: &str,
    data: Arc<RwLock<HashMap<String, HashMap<usize, (i32, i32)>>>>
) -> io::Result<()> {
    let data = data.read().unwrap();
    let mut file = File::create(file_path)?;

    for (chromosome, inner_map) in data.iter() {
        for (pos, (value_starts, value_ends)) in inner_map.iter() {
            writeln!(file, "{}\t{}\t{}\t{}", chromosome, pos, value_starts, value_ends)?;
        }
    }

    Ok(())
}

fn read_chromosome_data_from_file(
    file_path: &str,
) -> io::Result<Arc<RwLock<HashMap<String, HashMap<usize, (i32, i32)>>>>> {
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);
    let mut data: HashMap<String, HashMap<usize, (i32, i32)>> = HashMap::new();

    for line in reader.lines() {
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() == 4 {
                let chromosome = parts[0].to_string();
                let pos = parts[1].parse::<usize>().unwrap();
                let value_starts = parts[2].parse::<i32>().unwrap();
                let value_ends = parts[3].parse::<i32>().unwrap();

                let inner_map = data.entry(chromosome.clone()).or_insert_with(HashMap::new);
                inner_map.insert(pos, (value_starts, value_ends));
            }
        }
    }

    Ok(Arc::new(RwLock::new(data)))
}
//写入数据
/* 
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

    */
    fn assign_values_to_regions(
        chromosome_data: &HashMap<String, HashMap<usize, (i32, i32)>>,
        regions: &mut Vec<ChromosomeRegion>,
        chromosome: &str,
        file_index: usize
    ) -> io::Result<()> {
        if let Some(data) = chromosome_data.get(chromosome) {
            for (&position, &(value_starts, value_ends)) in data {
                let start_index = position - 1; // 转换为索引
                let end_index = position - 1; // 转换为索引
    
                // 将值分配到指定的区域
                regions[start_index].starts.push(value_starts);
                regions[end_index].ends.push(value_ends);
            }
        } else {
            eprintln!("Chromosome {} not found in data.", chromosome);
        }
    
        Ok(())
    }


//寻找有值对象的数量
fn count_non_empty_regions(regions: &[ChromosomeRegion]) -> usize {
    regions.iter().filter(|region| 
        !region.starts.is_empty() || !region.ends.is_empty()
    ).count()
}


//合并
fn generate_bedgraph(
    regions: &Vec<ChromosomeRegion>,
    chromosome: &str,
    output_bedgraph_path: &str,
    include_initial_line: bool,
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
        if !region.starts.is_empty() || !region.ends.is_empty() {
            let end = i + 1; // 当前索引 + 1 作为 end 的值

            // 计算 starts 和 ends 的总和
            let starts_sum: i32 = region.starts.iter().sum();
            let ends_sum: i32 = region.ends.iter().sum();

            // 存储当前行
            lines.push(format!("{}\t{}\t{}\t{}", chromosome, prev_end, end, current_value));
            current_value += starts_sum - ends_sum;

            prev_end = end; // 更新 prev_end，为下一次循环做准备
        }
    }

    // 写入去除第一行后的内容
    for line in lines.iter().skip(if include_initial_line { 1 } else { 0 }) {
        writeln!(bedgraph_file, "{}", line)?;
    }

    Ok(())
}

/* 
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

*/





type ChromosomeData = HashMap<String, HashMap<usize, (i32, i32)>>;

fn process_chromosomes(
    chromosomes: &[&str],
    chromosome_file: &str,
    chromosome_data: &ChromosomeData,
    work_dir: &str,
    some_flag: bool,
) -> Result<()> {
    // 遍历染色体数据，顺序处理
    for &chromosome in chromosomes {
        // 读取染色体长度
        let chromosome_length = read_chromosome_size(chromosome_file, chromosome)
            .context("Failed to read chromosome size")?;

        // 创建 ChromosomeRegion 对象的组合
        let mut regions = create_chromosome_regions(chromosome_length as usize);

        println!("Created {} chromosome regions for {}", regions.len(), chromosome);

        // 给 regions 赋值
        assign_values_to_regions(&chromosome_data, &mut regions, chromosome, 0)
            .context("Failed to assign values from chromosome data")?;
        println!("Assigned values to regions for {}", chromosome);

        // 统计非空的 regions 数量
        let non_empty_count = count_non_empty_regions(&regions);
        println!("Number of non-empty regions for {}: {}", chromosome, non_empty_count);

        // 生成 bedgraph 文件
        let output_bedgraph_path = format!("{}/output_{}.bedgraph", work_dir, chromosome);
        generate_bedgraph(&regions, chromosome, &output_bedgraph_path, some_flag)
            .context("Failed to generate bedgraph file")?;
    }

    Ok(())
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

fn cleanup_files(bam_dir: &str) -> io::Result<()> {
    if let Ok(entries) = fs::read_dir(bam_dir) {
        for entry in entries {
            if let Ok(entry) = entry {
                let path = entry.path();
                if path.is_file() && (path.extension().and_then(|ext| ext.to_str()) == Some("bam") || path.extension().and_then(|ext| ext.to_str()) == Some("bai")) {
                    fs::remove_file(path)?;
                }
            }
        }
    }
    Ok(())
}

pub fn process_files(input_name: &str) -> Result<(), Box<dyn std::error::Error>> {
    let rt = tokio::runtime::Runtime::new()?;
    
    rt.block_on(async {
        let file_list_path = format!("/data/haocheng/data/ACAT/txt/{}.txt", input_name);
        let bam_dir = format!("/data/haocheng/data/bam/{}", input_name);
        
        download_files_from_list(&file_list_path, &bam_dir).await
    })?;
    
    // 定义 BAM 目录和文件路径
    let bam_dir = format!("/data/haocheng/data/bam/{}", input_name);

    let work_dir = format!("{}/work", bam_dir);

    // 检查或创建 work 目录
    if !Path::new(&work_dir).exists() {
        std::fs::create_dir(&work_dir)?;
    }
    let mut bam_files = get_bam_files(&bam_dir);
    
    // 创建 BAM 文件的索引
    for bam_file in &bam_files {
        index_bam(bam_file)?;
    }
    let first_bam = bam_files.remove(0); // 移除第一个元素

    let output_file1 = format!("{}/coverage1.bed", work_dir);

    let chromosomes = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ];

    // 创建临时文件
    let temp_output_file1 = format!("{}/temp_coverage1.bed", work_dir);




    if is_file_empty(&output_file1).unwrap_or(false) {
        get_coverage(&first_bam, &chromosomes, &temp_output_file1).unwrap();
        filter_chromosome_data(&temp_output_file1, &chromosomes, &output_file1).unwrap();
    } else {
        println!("{} is not empty. Skipping coverage calculation.", output_file1);
    }
    println!("bedtools done");

    let num_threads = 4;

// 初始化数据结构
let chromosome_data = Arc::new(RwLock::new(HashMap::new()));

// 处理数据文件
build_chromosome_data(&output_file1, Arc::clone(&chromosome_data)).unwrap();

while !bam_files.is_empty() {
    let first_bam = bam_files.remove(0); 
    let output_file1 = format!("{}/coverage1.bed", work_dir);
    let temp_output_file1 = format!("{}/temp_coverage1.bed", work_dir);

    if is_file_empty(&output_file1).unwrap_or(false) {
        get_coverage(&first_bam, &chromosomes, &temp_output_file1).unwrap();
        filter_chromosome_data(&temp_output_file1, &chromosomes, &output_file1).unwrap();
    } else {
        println!("{} is not empty. Skipping coverage calculation.", output_file1);
    }
    println!("bedtools done");
    build_chromosome_data(&output_file1, Arc::clone(&chromosome_data)).unwrap();
}












let chromosome_file = "/data/haocheng/data/hg38.chrom.sizes";
let chromosome_data_ref = chromosome_data.read().unwrap();






    // 调用封装好的函数
    process_chromosomes(
        &chromosomes,
        chromosome_file,
        &chromosome_data_ref,
        &work_dir,
        //num_threads,
        true
    )?;



println!("All processing complete!");
// 将数据写入文件
//let chromosome_data_path = format!("{}/chromosome_data.txt", work_dir);
//write_chromosome_data_to_file(&chromosome_data_path, Arc::clone(&chromosome_data))?;

// 从文件中读取数据
//let chromosome_data = read_chromosome_data_from_file(&chromosome_data_path)?;

 //let chromosome_file = "/data/haocheng/data/hg38.chrom.sizes";





    let bam_dir = format!("/data/haocheng/data/bam/{}", input_name );
    cleanup_files(&bam_dir )?;
    println!("Cleanup completed successfully.");

    Ok(())
}




fn main() -> Result<()> {
   // let param = "GM_test"; // 这里替换为实际参数
   // process_files(param);
  
    let param = "GM_2"; // 这里替换为实际参数
    process_files(param);
    let param = "GM_3"; // 这里替换为实际参数
    process_files(param);

    let param = "HG_1"; // 这里替换为实际参数
    process_files(param);
    let param = "HG_2"; // 这里替换为实际参数
    process_files(param);
    let param = "HG_3"; // 这里替换为实际参数
    process_files(param);



    Ok(())

}
