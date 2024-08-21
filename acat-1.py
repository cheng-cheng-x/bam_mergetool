import pyBigWig
'''
# 打开 BigWig 文件
bw = pyBigWig.open('/data/haocheng/data/ACAT/k562/ENCFF798GCW.bigWig')
/data/haocheng/data/ACAT/k562/ENCFF920FUW.bigWig
# 获取特定染色体的前10个区间的数据
chrom = 'chr1'  # 替换为你感兴趣的染色体
intervals = bw.intervals(chrom)

# 打印前10个区间的数据
for i, entry in enumerate(intervals):
    if i >= 10:  # 只打印前10个区间
        break
    start, end, value = entry
    print(f"Chromosome: {chrom}, Start: {start}, End: {end}, Value: {value}")

# 关闭文件
bw.close()
'''
import pyBigWig
import subprocess

def normalize_bigwig(input_bw_path, output_bw_path, chrom_sizes_path, max_value=1000):
    # 打开输入 BigWig 文件
    bw = pyBigWig.open(input_bw_path)
    
    # 读取染色体大小信息
    chrom_sizes = {}
    with open(chrom_sizes_path) as f:
        for line in f:
            chrom, size = line.strip().split()
            chrom_sizes[chrom] = int(size)
    
    # 读取所有区间数据
    all_intervals = []
    for chrom in bw.chroms():
        intervals = bw.intervals(chrom)
        if intervals:
            all_intervals.extend([(chrom, interval[0], interval[1], interval[2]) for interval in intervals])
    
    # 计算整体最大值
    values = [interval[3] for interval in all_intervals]
    max_val = max(values)
    scale_factor = max_value / max_val
    
    # 创建输出 BigWig 文件
    bw_out = pyBigWig.open(output_bw_path, "w")
    bw_out.addHeader(list(chrom_sizes.items()))
    
    # 归一化数据并写入输出 BigWig 文件
    all_intervals.sort(key=lambda x: (x[0], x[1]))  # 按染色体和起始位置排序
    for chrom, start, end, value in all_intervals:
        bw_out.addEntries([chrom], [start], ends=[end], values=[value * scale_factor])
    
    bw.close()
    bw_out.close()
    
    return max_val

# 定义文件路径
save_path = "/data/haocheng/data/hg38.chrom.sizes"
print(f"使用已下载的文件: {save_path}")

file1 = "/data/haocheng/data/ACAT/k562/ENCFF798GCW.bigWig"
file2 = "/data/haocheng/data/ACAT/k562/ENCFF920FUW.bigWig"
output_bw = "/data/haocheng/data/result/k562/merged.bw"

# 归一化 BigWig 文件
file1_bw_normalized = "/data/haocheng/data/result/k562/ENCFF798GCW_normalized.bigWig"
file2_bw_normalized = "/data/haocheng/data/result/k562/ENCFF920FUW_normalized.bigWig"

max_val_file1 = normalize_bigwig(file1, file1_bw_normalized, save_path)
max_val_file2 = normalize_bigwig(file2, file2_bw_normalized, save_path)

print(f"文件 {file1} 的最大值: {max_val_file1}")
print(f"文件 {file2} 的最大值: {max_val_file2}")

# 运行 bigWigMerge 命令
subprocess.run(["bigWigMerge", file1_bw_normalized, file2_bw_normalized, output_bw])

# 打开合并后的 BigWig 文件
bw_merged = pyBigWig.open(output_bw)

# 读取所有区间数据
all_intervals_merged = []
for chrom in bw_merged.chroms():
    intervals = bw_merged.intervals(chrom)
    if intervals:
        all_intervals_merged.extend(intervals)

# 计算合并后的最大值
values_merged = [interval[2] for interval in all_intervals_merged]
max_val_merged = max(values_merged)

print(f"合并后的 BigWig 文件的最大值: {max_val_merged}")

bw_merged.close()


file1 = "/data/haocheng/data/ACAT/k562/ENCFF798GCW.bigWig"
file2 = "/data/haocheng/data/ACAT/k562/ENCFF920FUW.bigWig"
output_file = "/data/haocheng/data/result/k562/merged.bw"