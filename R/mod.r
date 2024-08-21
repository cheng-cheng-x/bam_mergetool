library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(S4Vectors)
library(BiocGenerics)
library(stats4)
library(BiocParallel)

# 读取两个BigWig文件
bw1 <- import("/data/haocheng/data/ACAT/k562/ENCFF798GCW.bigWig")
bw2 <- import("/data/haocheng/data/ACAT/k562/ENCFF920FUW.bigWig")

#分别归一化
# 提取score列
scores <- score(bw1)

# 归一化到0-100范围
normalized_scores <- (scores - min(scores)) / (max(scores) - min(scores)) * 1000

# 将归一化后的scores赋值回bw对象
score(bw1) <- normalized_scores

# 提取score列
scores <- score(bw2)

# 归一化到0-100范围
normalized_scores <- (scores - min(scores)) / (max(scores) - min(scores)) * 1000

# 将归一化后的scores赋值回bw对象
score(bw2) <- normalized_scores
# 清除不需要的对象以释放内存
rm(scores, normalized_scores)
gc()  #

chromosomes <- c("chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
bw1_chr <- subset(bw1, seqnames(bw1) %in% chromosomes)
bw2_chr <- subset(bw2, seqnames(bw2) %in% chromosomes)
rm(bw1,bw2)
gc()

# 找到两个BigWig文件中重叠的区间
overlaps <- findOverlaps(bw1_chr, bw2_chr)
print("找到重叠区间完成")

# 获取重叠区间的具体范围
overlap_ranges <- pintersect(ranges(bw1_chr)[queryHits(overlaps)], ranges(bw2_chr)[subjectHits(overlaps)])
print("获取重叠区间的具体范围完成")
print(overlap_ranges)
#print(overlaps)
# 接下来使用new_ranges创建GRanges对象

gc()
result <- GRanges(seqnames = seqnames(bw1_chr)[queryHits(overlaps)], ranges = overlap_ranges)
print(result)
# 提取重叠区间的分数并计算和
bw1_scores <- score(bw1_chr)[queryHits(overlaps)]
bw2_scores <- score(bw2_chr)[subjectHits(overlaps)]
rm(bw1_scores,bw2_scores,result_scores)
combined_scores <- (score(bw1_chr)[queryHits(overlaps)] + score(bw2_chr)[subjectHits(overlaps)]) / 2
print("提取重叠区间的分数并计算和完成")
#rm(overlaps)
gc()

# 将分数赋值给result对象
score(result) <- combined_scores
print("将分数赋值给result对象完成")
print(result)
gc()
# 提取bw1中不重叠的部分
non_overlap_indices_bw1 <- setdiff(seq_along(bw1_chr), queryHits(overlaps))
non_overlap_bw1 <- bw1_chr[non_overlap_indices_bw1]
print("提取bw1中不重叠的部分完成")

# 提取bw2中不重叠的部分
non_overlap_indices_bw2 <- setdiff(seq_along(bw2_chr), subjectHits(overlaps))
non_overlap_bw2 <- bw2_chr[non_overlap_indices_bw2]
print("提取bw2中不重叠的部分完成")

# 将不重叠的部分添加到result中
result <- c(result, non_overlap_bw1, non_overlap_bw2)
print("将不重叠的部分添加到result中完成")

# 清除不需要的对象以释放内存
rm(bw1_chr, bw2_chr, overlaps, overlap_ranges, breakpoints, new_ranges, bw1_scores, bw2_scores, combined_scores, result_scores, non_overlap_indices_bw1, non_overlap_bw1, non_overlap_indices_bw2, non_overlap_bw2)
gc()
print("清除不需要的对象完成")
print(result)

# 筛选特定序列名称 chr1 的区间
#result_chr <- subset(result, seqnames(result) == "chr1")
for (chr in chromosomes) {
  result_chr <- subset(result, seqnames(result) == chr)
  export(result_chr, paste0("/data/haocheng/result/HG/2_HG_new", chr, "_result.bigWig"))
  print(paste("导出染色体", chr, "的数据完成"))
  # 删除当前染色体的结果以释放内存
  rm(result_chr)
  result <- result[seqnames(result) != "chr1"]

  # 强制进行垃圾回收
  gc()
  print(paste("染色体", chr, "的数据已删除并释放内存"))
}