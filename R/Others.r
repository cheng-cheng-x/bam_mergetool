
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
library(scales)
# 初始化计数器
count <- 0
batch_size <- 6



# 定义BigWig文件路径
bw2_files <- list.files("/data/haocheng/data/ACAT/Others/", pattern = "*.bigWig", full.names = TRUE)

bw1 <- import(bw2_files[1])
bw2 <- import(bw2_files[2])

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
bw1_chr <- subset(bw1, seqnames(bw1) %in% chromosomes)
bw2_chr <- subset(bw2, seqnames(bw2) %in% chromosomes)

print("筛选出指定染色体的区域")
rm(bw1,bw2)
gc()
#分别归一化
# 提取score列
scores <- score(bw1_chr)

# 归一化到0-1000范围
normalized_scores <- rescale(scores, to = c(0, 1000))
# 将归一化后的scores赋值回bw对象
score(bw1_chr) <- normalized_scores

# 提取score列
scores <- score(bw2_chr)

# 归一化到0-1000范围
normalized_scores <- rescale(scores, to = c(0, 1000))
# 将归一化后的scores赋值回bw对象
score(bw2_chr) <- normalized_scores
# 清除不需要的对象以释放内存
rm(scores, normalized_scores)
gc()  #
# 设置并行参数
# 找到两个BigWig文件中重叠的区间
overlaps <- findOverlaps(bw1_chr, bw2_chr)
print("找到重叠区间完成")

# 获取重叠区间的具体范围
query_hits <- queryHits(overlaps)
subject_hits <- subjectHits(overlaps)

# 预先获取ranges和seqnames
bw1_ranges <- ranges(bw1_chr)
bw2_ranges <- ranges(bw2_chr)
bw1_seqnames <- seqnames(bw1_chr)
bw2_seqnames <- seqnames(bw2_chr)

# 计算重叠区间
overlap_ranges <- pintersect(bw1_ranges[query_hits], bw2_ranges[subject_hits])
print("获取重叠区间的具体范围完成")

# 创建GRanges对象
result <- GRanges(seqnames = bw1_seqnames[query_hits], ranges = overlap_ranges)

# 提取重叠区间的分数
bw1_scores <- score(bw1_chr)[query_hits]
bw2_scores <- score(bw2_chr)[subject_hits]

# 计算加权后的分数
n <- count + 2
m <- count + 1
combined_scores <- (bw1_scores * (m / n)) + (bw2_scores / n)
print("提取重叠区间的分数并计算加权和完成")

# 将分数赋值给result对象
score(result) <- combined_scores
print("将分数赋值给result对象完成")

# 清除不需要的对象以释放内存
rm(bw1_scores, bw2_scores, overlaps)
gc()  # 触发垃圾回收
print("清理中间结果完成")

# 提取bw1中不重叠的部分
non_overlap_indices_bw1 <- setdiff(seq_along(bw1_chr), query_hits)
non_overlap_bw1 <- bw1_chr[non_overlap_indices_bw1]
print("提取bw1中不重叠的部分完成")

# 提取bw2中不重叠的部分
non_overlap_indices_bw2 <- setdiff(seq_along(bw2_chr), subject_hits)
non_overlap_bw2 <- bw2_chr[non_overlap_indices_bw2]
print("提取bw2中不重叠的部分完成")
rm(non_overlap_indices_bw1,non_overlap_indices_bw2)
gc()
# 将不重叠的部分添加到result中
result <- c(result, non_overlap_bw1, non_overlap_bw2)
print("将不重叠的部分添加到result中完成")

# 清除不需要的对象以释放内存
rm(non_overlap_bw1, non_overlap_bw2, query_hits, subject_hits)
gc()  # 触发垃圾回收
print("清理最终结果完成")
print(result)
# 移除已分配给bw1和bw2的文件路径
bw2_files <- bw2_files[-c(1, 2)]
#export(result, paste0("/data/haocheng/data/result/Others/Others_result.bigWig"))
count<-count+1

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
for (bw2_file_index in seq_along(bw2_files)) {
  bw2_file <- bw2_files[bw2_file_index]
  print(paste("处理第",count+1,"文件:", bw2_file))
  
  # 读BigWig文件
  bw2 <- import(bw2_file)
 
  print("读取bw2数据成功")
    bw2_chr <- subset(bw2, seqnames(bw2) %in% chromosomes)
    rm(bw2)
    gc()
  # 提取score列
  scores <- score(bw2_chr)
  
 # 归一化到0-1000范围
normalized_scores <- rescale(scores, to = c(0, 1000)) 
  
  # 将归一化后的scores赋值回bw对象
  score(bw2_chr) <- normalized_scores
  print("归一化成功")
  # 清除不需要的对象以释放内存
  rm(scores, normalized_scores)
  gc()



  print("筛选bw_2")
  
  gc()
  bw1_chr <-result

    print("筛选bw_1")
    rm(result)
    print("删除result")
  
  # 找到两个BigWig文件中重叠的区间
overlaps <- findOverlaps(bw1_chr, bw2_chr)
print("找到重叠区间完成")

# 获取重叠区间的具体范围
query_hits <- queryHits(overlaps)
subject_hits <- subjectHits(overlaps)

# 预先获取ranges和seqnames
bw1_ranges <- ranges(bw1_chr)
bw2_ranges <- ranges(bw2_chr)
bw1_seqnames <- seqnames(bw1_chr)
bw2_seqnames <- seqnames(bw2_chr)

# 计算重叠区间
overlap_ranges <- pintersect(bw1_ranges[query_hits], bw2_ranges[subject_hits])
print("获取重叠区间的具体范围完成")

# 创建GRanges对象
result <- GRanges(seqnames = bw1_seqnames[query_hits], ranges = overlap_ranges)

# 提取重叠区间的分数
bw1_scores <- score(bw1_chr)[query_hits]
bw2_scores <- score(bw2_chr)[subject_hits]

# 计算加权后的分数
n <- count + 2
m <- count + 1
combined_scores <- (bw1_scores * (m / n)) + (bw2_scores / n)
print("提取重叠区间的分数并计算加权和完成")

# 将分数赋值给result对象
score(result) <- combined_scores
print("将分数赋值给result对象完成")

# 清除不需要的对象以释放内存
rm(bw1_scores, bw2_scores, overlaps)
gc()  # 触发垃圾回收
print("清理中间结果完成")

# 提取bw1中不重叠的部分
non_overlap_indices_bw1 <- setdiff(seq_along(bw1_chr), query_hits)
non_overlap_bw1 <- bw1_chr[non_overlap_indices_bw1]
print("提取bw1中不重叠的部分完成")

# 提取bw2中不重叠的部分
non_overlap_indices_bw2 <- setdiff(seq_along(bw2_chr), subject_hits)
non_overlap_bw2 <- bw2_chr[non_overlap_indices_bw2]
print("提取bw2中不重叠的部分完成")
rm(non_overlap_indices_bw1,non_overlap_indices_bw2)
gc()
# 将不重叠的部分添加到result中
result <- c(result, non_overlap_bw1, non_overlap_bw2)
print("将不重叠的部分添加到result中完成")

# 清除不需要的对象以释放内存
rm(non_overlap_bw1, non_overlap_bw2, query_hits, subject_hits)
gc()  # 触发垃圾回收
print("清理最终结果完成")
print(result)

 # 更新计数器
  count <- count + 1
  
  # 从bw2_files中排除当前处理的文件，更新列表
  bw2_files <- setdiff(bw2_files, bw2_file)

}
print(bw2_files)
export(result, paste0("/data/haocheng/data/result/Others/Others_result.bigWig"))
