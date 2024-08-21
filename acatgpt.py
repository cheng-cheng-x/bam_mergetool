import numpy as np
import pyBigWig
import subprocess
import gc
import os
from concurrent.futures import ThreadPoolExecutor

# 读取两个BigWig文件
bw1_path = "/data/haocheng/data/ACAT/k562/ENCFF798GCW.bigWig"
bw2_path = "/data/haocheng/data/ACAT/k562/ENCFF920FUW.bigWig"

# 定义需要的染色体
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
               "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

# 读取染色体的分数并返回
def read_scores(bw, chrom):
    length = bw.chroms()[chrom]  # 获取染色体长度
    scores = bw.values(chrom, 0, length)  # 提取整个染色体的分数
    return scores

# 归一化分数
def normalize_scores(bw):
    all_scores = []

    # 使用 ThreadPoolExecutor 并行读取分数
    with ThreadPoolExecutor(max_workers=100) as executor:  # 设置适当的线程数
        futures = {executor.submit(read_scores, bw, chrom): chrom for chrom in chromosomes if chrom in bw.chroms()}
        for future in futures:
            chrom = futures[future]
            scores = future.result()
            if scores is not None:
                all_scores.extend(scores)  # 收集所有分数

    # 计算整体的最小值和最大值
    min_score = np.min(all_scores)
    max_score = np.max(all_scores)
    range_score = max_score - min_score

    normalized_scores = {}
    for chrom in chromosomes:
        if chrom in bw.chroms():  # 检查染色体是否在 BigWig 文件中
            print(f"正在归一化染色体: {chrom}")
            scores = read_scores(bw, chrom)  # 提取整个染色体的分数
            if scores is not None:
                if range_score == 0:  # 避免除以零
                    normalized_scores[chrom] = np.zeros_like(scores)
                else:
                    normalized_scores[chrom] = (scores - min_score) / range_score * 1000

    return normalized_scores

# 打开 BigWig 文件
bw1 = pyBigWig.open(bw1_path)
bw2 = pyBigWig.open(bw2_path)

# 归一化两个 BigWig 文件
print("开始归一化第一个 BigWig 文件")
normalized_bw1_scores = normalize_scores(bw1)
print("第一个 BigWig 文件归一化完成")

print("开始归一化第二个 BigWig 文件")
normalized_bw2_scores = normalize_scores(bw2)
print("第二个 BigWig 文件归一化完成")

# 清除不需要的对象以释放内存
del bw1, bw2
gc.collect()

# 使用wiggletools进行叠加并导出结果为一个BigWig文件
output_file = "/data/haocheng/data/result/k562/combined_result.bigWig"

# 创建临时文件路径
temp_bw1_path = "/data/haocheng/data/result/k562/temp_bw1.bigWig"
temp_bw2_path = "/data/haocheng/data/result/k562/temp_bw2.bigWig"

# 删除已存在的临时文件
for temp_path in [temp_bw1_path, temp_bw2_path]:
    if os.path.exists(temp_path):
        os.remove(temp_path)

# 创建临时 BigWig 文件并写入归一化后的结果
with pyBigWig.open(temp_bw1_path, "w") as temp_bw1:
    for chrom, scores in normalized_bw1_scores.items():
        temp_bw1.addHeader([(chrom, len(scores))])
        temp_bw1.addEntries(chrom, np.arange(len(scores)), values=scores)

with pyBigWig.open(temp_bw2_path, "w") as temp_bw2:
    for chrom, scores in normalized_bw2_scores.items():
        temp_bw2.addHeader([(chrom, len(scores))])
        temp_bw2.addEntries(chrom, np.arange(len(scores)), values=scores)

# 使用wiggletools进行叠加
command = f"wiggletools add {temp_bw1_path} {temp_bw2_path} > {output_file}"
print(f"正在执行命令: {command}")
subprocess.run(command, shell=True, check=True)

# 删除临时文件
os.remove(temp_bw1_path)
os.remove(temp_bw2_path)

print("所有数据已合并并导出为一个BigWig文件完成")
