
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(S4Vectors)
library(BiocGenerics)
library(stats4)
library(BiocParallel)
library(parallel)
library(rstudioapi)
restart_and_run <- function(script_path) {
    temp_file <- tempfile()
    writeLines(script_path, temp_file)
    rstudioapi::restartSession()
    .Last <- function() {
        load("temp_data.RData")
        script_path <- readLines(temp_file)
        source(script_path)
    }
}
# 使用所有可用核心减去1个核心
num_cores <- 4
cl <- makeCluster(num_cores, type = "FORK")

# 获取所有已导出的文件列表
files <- list.files("/data/haocheng/data/result/GM/", pattern = "20GM_result_.*\\.bigWig", full.names = TRUE)

# 初始化一个空的 GRanges 对象
result <- GRanges()

# 并行读取数据并组合成一个结果对象
import_result <- parLapply(cl, seq_along(files), function(i) {
  file <- files[i]
  import(file)
})

# 合并所有导入的结果
result <- do.call(c, import_result)

# 停止集群
stopCluster(cl)

# 清除不需要的数据
rm(import_result)
gc()
seqinfo <- Seqinfo(genome="hg38")
# 初始化计数器
count <- 21
batch_size <- 10
bw2_files <- list.files("/data/haocheng/data/ACAT/GM/", pattern = "*.bigWig", full.names = TRUE)
# 移除已分配给bw1和bw2的文件路径
bw2_files <- bw2_files[-c(1:22)]

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
for (bw2_file_index in seq_along(bw2_files)) {
  bw2_file <- bw2_files[bw2_file_index]
  print(paste("处理第",count+1,"文件:", bw2_file))
  
  # 读BigWig文件
  bw2 <- import(bw2_file)
 
  print("读取bw2数据成功")
  # 提取score列
  scores <- score(bw2)
  
  # 归一化到0-100范围
  normalized_scores <- (scores - min(scores)) / (max(scores) - min(scores)) * 1000
  
  # 将归一化后的scores赋值回bw对象
  score(bw2) <- normalized_scores
  print("归一化成功")
  # 清除不需要的对象以释放内存
  rm(scores, normalized_scores)
  gc()
  bw2_chr <- subset(bw2, seqnames(bw2) %in% chromosomes)
seqlevels(bw2_chr) <- seqlevels(seqinfo)
seqinfo(bw2_chr) <- seqinfo

  print("筛选bw_2")
  rm(bw2)
  gc()
  bw1_chr <-result
  seqlevels(bw1_chr) <- seqlevels(seqinfo)
  seqinfo(bw1_chr) <- seqinfo
    print("筛选bw_1")
    rm(result)
    print("删除result")
  
# 假设 bw1_chr 和 bw2_chr 已经是 GRangesList 对象，如果不是，请转换为 GRangesList 对象
bw1_chr <- GRangesList(bw1_chr)
bw2_chr <- GRangesList(bw2_chr)

# 建立并行集群
cl <- makeCluster(4, type = "FORK")
clusterExport(cl, c("bw1_chr", "bw2_chr"))
clusterEvalQ(cl, {
  library(GenomicRanges)
  library(BiocParallel)
})

print("开始计算重叠区域")
# 使用 parLapply 并行执行 findOverlaps
overlaps <- parLapply(cl, seq_along(bw1_chr), function(i) {
  findOverlaps(bw1_chr[[i]], bw2_chr[[i]])
})
print("重叠区域计算完成")
clusterExport(cl, c("overlaps"))
print("开始计算重叠区间的具体范围")
# 使用 parLapplyLB 并行计算重叠区间的具体范围并创建GRanges对象
overlap_granges <- parLapplyLB(cl, seq_along(bw1_chr), function(i) {
  # 提取每对 bw1_chr[i] 和 bw2_chr[i] 的 ranges
  ranges1 <- ranges(bw1_chr[[i]])
  ranges2 <- ranges(bw2_chr[[i]])
  
  # 提取对应的 overlaps[[i]] 的 Hits 对象
  hits <- overlaps[[i]]
  
  # 计算交集
  overlap <- pintersect(ranges1[queryHits(hits)], ranges2[subjectHits(hits)])
  
  # 根据公式计算分数
  count <- length(hits)
  n <- count + 2
  m <- count + 1
  combined_scores <- (score(bw1_chr[[i]])[queryHits(hits)]) * (m / n) + (score(bw2_chr[[i]])[subjectHits(hits)]) / n
  
  # 构建 GRanges 对象
  GRanges(
    seqnames = seqnames(bw1_chr[[i]])[queryHits(hits)],
    ranges = overlap,
    strand = strand(bw1_chr[[i]])[queryHits(hits)],
    score = combined_scores
  )
})

print("重叠区间计算完成")

# 清除不需要的对象以释放内存

gc()

print("开始处理非重叠区域")
# 非重叠区间处理
non_overlap_granges1 <- parLapplyLB(cl, seq_along(bw1_chr), function(i) {
  hits <- overlaps[[i]]
  non_overlaps <- setdiff(seq_along(bw1_chr[[i]]), queryHits(hits))
  
  GRanges(
    seqnames = seqnames(bw1_chr[[i]])[non_overlaps],
    ranges = ranges(bw1_chr[[i]])[non_overlaps],
    strand = strand(bw1_chr[[i]])[non_overlaps],
    score = score(bw1_chr[[i]])[non_overlaps]
  )
})

non_overlap_granges2 <- parLapplyLB(cl, seq_along(bw2_chr), function(i) {
  hits <- overlaps[[i]]
  non_overlaps <- setdiff(seq_along(bw2_chr[[i]]), subjectHits(hits))
  
  GRanges(
    seqnames = seqnames(bw2_chr[[i]])[non_overlaps],
    ranges = ranges(bw2_chr[[i]])[non_overlaps],
    strand = strand(bw2_chr[[i]])[non_overlaps],
    score = score(bw2_chr[[i]])[non_overlaps]
  )
})
print("非重叠区间处理完成")

# 清除不需要的对象以释放内存
rm(bw1_chr, bw2_chr, overlaps)
gc()

stopCluster(cl)  # 停止并行集群
gc()

# 整合所有的 overlap_granges 和 non_overlap_granges 到一个 GRanges 对象中
result_overlap <- do.call(c, overlap_granges)
result_non_overlap1 <- do.call(c, non_overlap_granges1)
result_non_overlap2 <- do.call(c, non_overlap_granges2)

# 清除不需要的对象以释放内存
rm(overlap_granges, non_overlap_granges1, non_overlap_granges2)
gc()

result <- c(result_overlap, result_non_overlap1, result_non_overlap2)

# 清除不需要的对象以释放内存
rm(result_overlap, result_non_overlap1, result_non_overlap2)
print(paste("清除不需要的对象完成", count))
gc()
# 设置 seqinfo
seqinfo(result) <- seqinfo
print("处理完成")

# 输出 result
print(result)
  # 更新计数器
  count <- count + 1
  print(result)
  
 
# 从bw2_files中排除当前处理的文件，更新列表
bw2_files <- setdiff(bw2_files, bw2_file)
}
num_cores <- 4  # 使用所有可用核心减去1个核心
cl <- makeCluster(num_cores, type = "FORK")
# 导出必要的对象到集群
clusterExport(cl, varlist = c("result"))

parLapply(cl, unique(seqnames(result)), function(chr) {
  export(result[seqnames(result) == chr], paste0("/data/haocheng/data/result/GM/",count,"GM_result_", chr, ".bigWig"))}
  restart_and_run("/data/haocheng/project/GM2.r")