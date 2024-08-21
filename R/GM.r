
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



# 定义BigWig文件路径
bw2_files <- list.files("/data/haocheng/data/ACAT/GM/", pattern = "*.bigWig", full.names = TRUE)

bw1 <- parLapply(cl, list(bw2_files[1]), import)[[1]]
bw2 <- parLapply(cl, list(bw2_files[2]), import)[[1]]


chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
bw1_chr <- subset(bw1, seqnames(bw1) %in% chromosomes)
bw2_chr <- subset(bw2, seqnames(bw2) %in% chromosomes)
print("筛选出指定染色体的区域")
rm(bw1,bw2)
gc()
#分别归一化
# 提取score列
scores <- score(bw1_chr)

# 归一化到0-100范围
normalized_scores <- (scores - min(scores)) / (max(scores) - min(scores)) * 1000

# 将归一化后的scores赋值回bw对象
score(bw1) <- normalized_scores

# 提取score列
scores <- score(bw2_chr)

# 归一化到0-100范围
normalized_scores <- (scores - min(scores)) / (max(scores) - min(scores)) * 1000

# 将归一化后的scores赋值回bw对象
score(bw2) <- normalized_scores
# 清除不需要的对象以释放内存
rm(scores, normalized_scores)
gc()  #

# 找到两个BigWig文件中重叠的区间
overlaps <- findOverlaps(bw1_chr, bw2_chr)
print("找到重叠区间完成")


bw1_chr <- uncompress(bw1_chr)
bw2_chr <- uncompress(bw2_chr)

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

# 移除已分配给bw1和bw2的文件路径
bw2_files <- bw2_files[-c(1, 2)]


num_cores <- 4  # 使用所有可用核心减去1个核心
cl <- makeCluster(num_cores, type = "FORK")
# 导出必要的对象到集群
clusterExport(cl, varlist = c("result"))

parLapply(cl, unique(seqnames(result)), function(chr) {
  export(result[seqnames(result) == chr], paste0("/data/haocheng/data/result/GM/GM_result_", chr, ".bigWig"))
})

print("导出染色体的数据完成")

stopCluster(cl)
  #export(result, paste0("/data/haocheng/data/result/GM/GM_result.bigWig"))
  print(paste("导出染色体的数据完成"))
  










#rm(bw1,bw2)
gc()



# 定义染色体列表
chromosomes <- c( "chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

# 初始化计数器
count <- 0
batch_size <- 6


for (bw2_file_index in seq_along(bw2_files)) {
  bw2_file <- bw2_files[bw2_file_index]
  print(paste("处理文件:", bw2_file))
  
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
  n<-count+3
  m<-count+2
  # 提取重叠区间的分数并计算和
  combined_scores <- (score(bw1_chr)[queryHits(overlaps)]) * (m / n)+ (score(bw2_chr)[subjectHits(overlaps)])/ n
  print("提取重叠区间的分数并计算和完成")
  #rm(overlaps)
  gc()
  
  # 将分数赋值给result对象
  score(result) <- combined_scores
  print("将分数赋值给result对象完成")
  print(result)
  rm(combined_scores)
  gc()
  # 提取bw1中不重叠的部分
  non_overlap_bw1 <- bw1_chr[setdiff(seq_along(bw1_chr), queryHits(overlaps))]
  print("提取bw1中不重叠的部分完成")
  
  # 提取bw2中不重叠的部分
  non_overlap_bw2 <- bw2_chr[setdiff(seq_along(bw2_chr), subjectHits(overlaps))]
  print("提取bw2中不重叠的部分完成")
  
  rm(bw1_chr, bw2_chr, overlaps, overlap_ranges)
  gc()
  # 将不重叠的部分添加到result中
  result <- c(result, non_overlap_bw1, non_overlap_bw2)
  print("将不重叠的部分添加到result中完成")
  
  # 清除不需要的对象以释放内存
  rm(bw1_chr, bw2_chr, overlaps, overlap_ranges, breakpoints, new_ranges, bw1_scores, bw2_scores, combined_scores, result_scores, non_overlap_indices_bw1, non_overlap_bw1, non_overlap_indices_bw2, non_overlap_bw2)
  gc()
  # 更新计数器
  count <- count + 1
  print(paste("清除不需要的对象完成", count))
  print(result)
  

  
  # 每处理6个文件保存一次结果
  if (count %% batch_size == 0) {
    for (chr in chromosomes) {
      result_chr <- subset(result, seqnames(result) == chr)
      export(result_chr, paste0("C:/Users/Administrator/haocheng/acat/result/HG/2_HG_new", chr, "_result.bigWig"))
      print(paste("导出染色体", chr, "的数据完成"))
      # 删除当前染色体的结果以释放内存
      rm(result_chr)
      # 强制进行垃圾回收
      gc()
      print(paste("染色体", chr, "的数据已删除并释放内存"))
    }
  }
  # 从bw2_files中排除当前处理的文件，更新列表
  bw2_files <- setdiff(bw2_files, bw2_file)

}
 

# 所有染色体处理完成后，进行数据导出
for (chr in chromosomes) {
  result_chr <- subset(result, seqnames(result) == chr)
  export(result_chr, paste0("C:/Users/Administrator/haocheng/acat/result/HG/2_HG_new", chr, "_result.bigWig"))
  print(paste("导出染色体", chr, "的数据完成"))
  # 删除当前染色体的结果以释放内存
  rm(result_chr)
  result <- result[seqnames(result) != "chr1"]
  
  # 强制进行垃圾回收
  gc()
  print(paste("染色体", chr, "的数据已删除并释放内存"))
}




num_cores <- 4  # 使用所有可用核心减去1个核心
cl <- makeCluster(num_cores, type = "FORK")
# 导出必要的对象到集群
clusterExport(cl, varlist = c("result"))

parLapply(cl, unique(seqnames(result)), function(chr) {
  export(result[seqnames(result) == chr], paste0("/data/haocheng/data/result/GM/GM_result_", chr, ".bigWig"))
})

print("导出染色体的数据完成")

stopCluster(cl)
