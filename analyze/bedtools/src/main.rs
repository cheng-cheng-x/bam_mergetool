use std::process::{Command, Stdio};
use std::fs::{File, metadata};  
use std::io::{self, Write,BufReader, BufRead};
use std::path::Path;
use std::fs;

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
        let status = Command::new("samtools")
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


fn get_coverage(bam_file: &str, chromosomes: &[&str], output_file: &str) -> io::Result<()> {
    // 打开输出文件，用于写入结果
    let output = File::create(output_file)?;

    // 执行 bedtools genomecov 命令以计算覆盖度
    let status = Command::new("sh")
    .arg("-c")
    .arg(format!(
        "source /data/haocheng/miniconda3/bin/activate rust && bedtools genomecov -ibam {} -bg | grep -E '^({})' > {}",
        bam_file,
        chromosomes.join("|"),
        output_file
    ))
    .status()?;


    if status.success() {
        println!("Coverage data for selected chromosomes written to {}", output_file);
    } else {
        eprintln!("Error in executing bedtools or grep.");
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
fn assign_values_to_regions(file_path: &str, regions: &mut Vec<ChromosomeRegion>, chromosome: &str, file_index: usize) -> io::Result<()> {
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);

    let mut found_chromosome = false; // 标记是否找到匹配的染色体

    for line in reader.lines() {
        let line = line.unwrap();
        let parts: Vec<&str> = line.split('\t').collect(); // 使用 '\t' 作为分隔符

        // 检查当前行的染色体
        if parts[0] == chromosome {
            found_chromosome = true; // 找到匹配的染色体

            let start = parts[1].parse::<usize>().unwrap();
            let end = parts[2].parse::<usize>().unwrap();

            let start_index = start - 1; // 转换为索引
            let end_index = end - 1; // 转换为索引

            let value = parts[3].parse::<i32>().unwrap();

            regions[start_index].starts[file_index].push(value);
            regions[end_index].ends[file_index].push(value);
        } else if found_chromosome {
            // 一旦已经找到匹配的染色体，再遇到其他染色体就跳出循环
            break;
        }
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
    output_bedgraph_path: &str
) -> io::Result<()> {
    let mut bedgraph_file = File::create(output_bedgraph_path)?;

    let mut prev_end = 1; // 初始化起始位置
    let mut current_value = 0;
    let mut lines: Vec<String> = Vec::new(); // 用于存储所有行
    // 初始化第一行
    //lines.push(format!("{}\t{}\t{}\t{}", chromosome, 1, 1, 0));

    for (i, region) in regions.iter().enumerate() {
        // 只要有一个 starts 或 ends 非空，即可开始处理
        if !region.starts[0].is_empty() || !region.starts[1].is_empty() || !region.ends[0].is_empty() || !region.ends[1].is_empty() {
            let end = i+1; // 当前索引 +1 作为 end 的值

            // 计算 value，忽略空 Vec
            let starts_0_sum: i32 = region.starts[0].iter().sum();
            let starts_1_sum: i32 = region.starts[1].iter().sum();
            let ends_0_sum: i32 = region.ends[0].iter().sum();
            let ends_1_sum: i32 = region.ends[1].iter().sum();

           // current_value += starts_0_sum + starts_1_sum - ends_0_sum - ends_1_sum;

            // 存储当前行
            lines.push(format!("{}\t{}\t{}\t{}", chromosome, prev_end, end, current_value));
            current_value += starts_0_sum + starts_1_sum - ends_0_sum - ends_1_sum;

            prev_end = end; // 更新 prev_end，为下一次循环做准备
        }
    }

    // 写入去除第一行后的内容
    for line in lines.iter().skip(1) {
        writeln!(bedgraph_file, "{}", line)?;
    }

    Ok(())
}


fn main() -> io::Result<()> {
    // 定义 BAM 目录和文件路径
    let bam_dir = "/data/haocheng/data/bam/k562";
    let work_dir = format!("{}/work", bam_dir);

    // 检查或创建 work 目录
    if !Path::new(&work_dir).exists() {
        std::fs::create_dir(&work_dir)?;
    }
    let bam_files = get_bam_files(bam_dir);
    let first_bam = &bam_files[0];
    let second_bam = &bam_files[1];
    let output_file1 = format!("{}/coverage1.bed", work_dir);
    let output_file2 = format!("{}/coverage2.bed", work_dir);

    let chromosomes = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ];

    // 处理第一个 BAM 文件
    if is_file_empty(&output_file1)? {
        get_coverage(&first_bam, &chromosomes, &output_file1)?;
    } else {
        println!("{} is not empty. Skipping coverage calculation.", output_file1);
    }

    // 处理第二个 BAM 文件
    if is_file_empty(&output_file2)? {
        get_coverage(&second_bam, &chromosomes, &output_file2)?;
    } else {
        println!("{} is not empty. Skipping coverage calculation.", output_file2);
    }

   let chromosome_file = "/data/haocheng/data/hg38.chrom.sizes"; // 这里应替换为想要的染色体长度文件路径
    
   let chromosome = "chr1"; // 替换想要查询的染色体名称
 
 // 读取染色体长度
   let chromosome_length = read_chromosome_size(chromosome_file, chromosome)?;
 
 // 创建 ChromosomeRegion 对象的组合
   let  mut regions = create_chromosome_regions(chromosome_length as usize);
 
 // 现在 regions 中包含 chromosome_length 个 ChromosomeRegion 对象
    println!("Created {} chromosome regions for {}", regions.len(), chromosome);
   

    assign_values_to_regions(&output_file1, &mut regions, chromosome, 0)?;
    assign_values_to_regions(&output_file2, &mut regions, chromosome, 1)?;
 

    println!("Assigned values to {} chromosome regions for {}", regions.len(), chromosome);


   
 // 统计非空的 regions 数量
    let non_empty_count = count_non_empty_regions(&regions);
    println!("Number of non-empty regions: {}", non_empty_count);
    

    let output_bedgraph_path = format!("{}/output_chr1.bedgraph", work_dir);
    generate_bedgraph(&regions, chromosome, &output_bedgraph_path)?;





    Ok(())

}
