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
library(doParallel)
# 初始化计数器
count <- 0
batch_size <- 6
# 定义GRCh38的seqinfo
seqinfo <- Seqinfo(genome="hg38")
# 定义BigWig文件路径
bw2_files <- list.files("/data/haocheng/data/ACAT/GM12878/", pattern = "*.bigWig", full.names = TRUE)

bw1 <- import(bw2_files[1])
bw2 <- import(bw2_files[2])

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
bw1_chr <- subset(bw1, seqnames(bw1) %in% chromosomes)
bw2_chr <- subset(bw2, seqnames(bw2) %in% chromosomes)
print("筛选出指定染色体的区域")
# 确保 seqlevels 一致
seqlevels(bw1_chr) <- seqlevels(seqinfo)
seqlevels(bw2_chr) <- seqlevels(seqinfo)

# 设置 seqinfo
seqinfo(bw1_chr) <- seqinfo
seqinfo(bw2_chr) <- seqinfo

rm(bw1, bw2)
gc()

# 分别归一化
# 提取score列
scores1 <- score(bw1_chr)
scores2 <- score(bw2_chr)

# 归一化到0-1000范围
normalized_scores1 <- (scores1 - min(scores1)) / (max(scores1) - min(scores1)) * 1000
normalized_scores2 <- (scores2 - min(scores2)) / (max(scores2) - min(scores2)) * 1000

# 将归一化后的scores赋值回bw对象
score(bw1_chr) <- normalized_scores1
score(bw2_chr) <- normalized_scores2

# 清除不需要的对象以释放内存
rm(scores1, scores2, normalized_scores1, normalized_scores2)
gc()

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
gc()
# 设置 seqinfo
seqinfo(result) <- seqinfo
print("处理完成")

# 输出 result
print(result)
bw2_files <- bw2_files[-c(1:2)]
 count <- count + 1

export(result, paste0("/data/haocheng/data/result/GM12878/GM12878_result.bigWig"))

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
  print("筛选bw_2")
  rm(bw2)
  gc()
  bw1_chr <-result
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
gc()
# 设置 seqinfo
seqinfo(result) <- seqinfo
print("处理完成")

# 输出 result
print(result)
count <- count + 1
# 从bw2_files中排除当前处理的文件，更新列表
  bw2_files <- setdiff(bw2_files, bw2_file)
 

# 导出为bigWig文件
export.bw(result, "/data/haocheng/data/result/A549/A549_result.bigWig")


}