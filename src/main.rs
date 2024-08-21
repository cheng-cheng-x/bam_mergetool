use bam_merge_tools::bam_merge;

fn main() -> Result<(), Box<dyn std::error::Error>> {

     // 定义输入参数
     let input_name = "GM_3";

     // 定义路径前缀
     let acat_txt_dir = "/data/haocheng/data/ACAT/txt";
     let bam_base_dir = "/data/haocheng/data/bam";
 
     // 基于 input_name 生成其他路径和参数
     let file_list_path = format!("{}/{}.txt", acat_txt_dir, input_name);
     let bam_dir = format!("{}/{}", bam_base_dir, input_name);
     let work_dir = format!("{}/work", bam_dir);
     let chromosome_file ="/data/haocheng/data/hg38.chrom.sizes";
 
     // 定义染色体列表
     let chromosomes = [
         "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
         "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
         "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
     ];
 
     // 调用 process_files 函数
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
