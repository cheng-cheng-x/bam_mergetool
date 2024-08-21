library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(S4Vectors)
library(BiocGenerics)
library(stats4)
library(BiocParallel)
library(foreach)
library(doParallel)

# 设置并行参数
num_cores <- 4  # 使用较少的核心数
registerDoParallel(num_cores)

# 读取两个BigWig文件
print("读取BigWig文件")
bw1 <- import("/data/haocheng/data/ACAT/k562/ENCFF798GCW.bigWig")
bw2 <- import("/data/haocheng/data/ACAT/k562/ENCFF920FUW.bigWig")

# 分别归一化
print("归一化BigWig文件")
normalize_scores <- function(bw) {
  scores <- score(bw)
  normalized_scores <- (scores - min(scores)) / (max(scores) - min(scores)) * 1000
  score(bw) <- normalized_scores
  return(bw)
}

bw1 <- normalize_scores(bw1)
bw2 <- normalize_scores(bw2)

# 保留指定的染色体
print("筛选染色体")
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
bw1 <- subset(bw1, seqnames(bw1) %in% chromosomes)
bw2 <- subset(bw2, seqnames(bw2) %in% chromosomes)

# 分批查找重叠区间
print("查找重叠区间")
find_overlaps_parallel <- function(bw1, bw2, chunk_size = 1e6) {  # 增大块大小，减少同时处理的块数
  num_chunks <- ceiling(length(bw1) / chunk_size)
  overlaps_list <- foreach(i = 1:num_chunks, .combine = c, .packages = c("GenomicRanges", "IRanges", "GenomeInfoDb", "S4Vectors")) %dopar% {
    start_index <- (i - 1) * chunk_size + 1
    end_index <- min(i * chunk_size, length(bw1))
    findOverlaps(bw1[start_index:end_index], bw2)
  }
  return(overlaps_list)
}

overlaps <- find_overlaps_parallel(bw1, bw2)

# 清理临时对象
rm(bw1, bw2)
gc()

# 分批处理重叠区间并计算分数
print("计算重叠区间分数")
split_overlaps <- split(overlaps, ceiling(seq_along(overlaps) / 50000))  # 分块大小改为50000
combined_scores_list <- foreach(overlap_chunk = split_overlaps, .combine = c, .packages = c("S4Vectors")) %dopar% {
  (score(bw1)[queryHits(overlap_chunk)] + score(bw2)[subjectHits(overlap_chunk)]) / 2
}

# 清理临时对象
rm(overlaps)
gc()

# 合并所有计算结果
print("合并计算结果")
combined_scores <- do.call(c, combined_scores_list)

# 清理临时对象
rm(combined_scores_list)
gc()

# 提取重叠区间并创建结果对象
print("创建结果对象")
result <- pintersect(bw1[queryHits(overlaps)], bw2[subjectHits(overlaps)])
score(result) <- combined_scores

# 清理临时对象
rm(combined_scores)
gc()

# 提取不重叠的部分
print("提取不重叠部分")
non_overlap_bw1 <- bw1[setdiff(seq_along(bw1), queryHits(overlaps))]
non_overlap_bw2 <- bw2[setdiff(seq_along(bw2), subjectHits(overlaps))]

# 清理临时对象
rm(queryHits, subjectHits)
gc()

# 合并结果
print("合并结果")
result <- c(result, non_overlap_bw1, non_overlap_bw2)

# 清理临时对象
rm(non_overlap_bw1, non_overlap_bw2)
gc()

# 导出结果
print("导出结果")
export(result, paste0("/data/haocheng/data/result/k562/k562_result.bigWig"))

# 清理内存
rm(result)
gc()

print("处理完成")
