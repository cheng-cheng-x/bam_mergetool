import os
import subprocess
import pandas as pd

bam_dir = "/data/haocheng/data/bam/GM"  # BAM 文件目录
bam_files = os.listdir(bam_dir)  # 列出所有 BAM 文件
chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]  # 选择的染色体

dataframes = []  # 存储 DataFrame 的列表

for bam_file in bam_files[:2]:  # 只取前两个 BAM 文件
    bam_path = os.path.join(bam_dir, bam_file)
    depth_data = []
    
    for chrom in chromosomes:
        # 使用 sambamba 计算深度，指定单个染色体
        command = ["sambamba", "depth", "base", bam_path, "-L", chrom]
        print(f"Running command: {' '.join(command)}")  # 打印调试信息
        
        try:
            depth_output = subprocess.check_output(command, universal_newlines=True)
            print(f"Output for {bam_file} - {chrom} received")  # 打印调试信息
            
            # 将深度输出存储为 DataFrame
            for line in depth_output.strip().split('\n')[1:]:  # 跳过表头
                fields = line.split('\t')
                depth_data.append((fields[0], fields[1], int(fields[3])))  # (CHROM, POS, DEPTH)
        
        except subprocess.CalledProcessError as e:
            print(f"Error processing {bam_file} - {chrom}: {e}")
    
    df = pd.DataFrame(depth_data, columns=['CHROM', 'POS', 'DEPTH'])
    dataframes.append(df)

# 确保至少有两个 DataFrame 可以进行处理
if len(dataframes) < 2:
    raise ValueError("Not enough BAM files processed to perform merging and normalization")

# 计算最大深度
max1 = dataframes[0]['DEPTH'].max()
max2 = dataframes[1]['DEPTH'].max()
norm_max = max(max1, max2)

# 计算 n 和 m
n = len(dataframes[0])
m = len(dataframes[1])

# 归一化并合并数据
dataframes[0]['DEPTH_NORM'] = (dataframes[0]['DEPTH'] / max1) * norm_max * (n / m)
dataframes[1]['DEPTH_NORM'] = (dataframes[1]['DEPTH'] / max2) * norm_max * (1 / m)

# 合并 DataFrame
merged_df = pd.merge(dataframes[0], dataframes[1], on=['CHROM', 'POS'], suffixes=('_file1', '_file2'), how='outer')
merged_df['COMBINED_DEPTH'] = merged_df['DEPTH_NORM_file1'].fillna(0) + merged_df['DEPTH_NORM_file2'].fillna(0)

# 打印合并后的 DataFrame
print(merged_df[['CHROM', 'POS', 'COMBINED_DEPTH']])
