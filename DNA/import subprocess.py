import os
import subprocess
import pandas as pd

# BAM 文件目录
bam_dir = "/data/haocheng/data/bam/GM"
bam_files = os.listdir(bam_dir)

# 设定要分析的染色体
chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# 用于存储 DataFrame 的列表
dataframes = []

# 检查 BAM 文件中包含的染色体
print("Checking chromosomes in BAM files...")
chromosomes_in_bam = set()
for bam_file in bam_files[:2]:  # 只取前两个 BAM 文件
    bam_path = os.path.join(bam_dir, bam_file)
    output = subprocess.check_output(
        ["samtools", "idxstats", bam_path],
        universal_newlines=True
    )
    for line in output.strip().split('\n'):
        chrom, length, count, _ = line.split()
        chromosomes_in_bam.add(chrom)

# 过滤染色体
filtered_chromosomes = [chrom for chrom in chromosomes if chrom in chromosomes_in_bam]
print(f"Filtered chromosomes: {filtered_chromosomes}")

# 计算深度并存储结果
print("Calculating depth for each BAM file...")
for bam_file in bam_files[:2]:  # 只处理前两个 BAM 文件
    bam_path = os.path.join(bam_dir, bam_file)
    for chrom in filtered_chromosomes:
        command = ["samtools", "depth", "-r", chrom, bam_path]
        print(f"Running command for {bam_file} on {chrom}...")
        try:
            depth_output = subprocess.check_output(
                command,
                universal_newlines=True,
                stderr=subprocess.STDOUT  # 捕获标准错误输出
            )
        except subprocess.CalledProcessError as e:
            print(f"Error processing {bam_file} on {chrom}: {e.output}")
            continue
        
        # 将深度输出存储为 DataFrame
        depth_data = []
        for line in depth_output.strip().split('\n'):
            fields = line.split('\t')
            depth_data.append((fields[0], fields[1], int(fields[2])))  # (CHROM, POS, DEPTH)
        
        df = pd.DataFrame(depth_data, columns=['CHROM', 'POS', 'DEPTH'])
        dataframes.append(df)
        print(f"Completed depth calculation for {bam_file} on {chrom}")

# 合并所有 DataFrame
if dataframes:
    final_df = pd.concat(dataframes, ignore_index=True)
    print("All depth calculations completed. Final DataFrame created.")
    # 你可以在这里进行进一步的数据处理或分析
    print(final_df)
else:
    print("No dataframes to concatenate.")
