import os
import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

# 定义 BAM 文件目录
bam_dir = "/data/haocheng/data/bam/GM/"
bam_files = sorted([f for f in os.listdir(bam_dir) if f.endswith('.bam')])  # 确保文件按字母顺序排序
first_bam = os.path.join(bam_dir, bam_files[0])  # sample1.bam
second_bam = os.path.join(bam_dir, bam_files[1])  # sample2.bam

# 定义要分析的染色体列表
chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# 使用 sambamba 计算深度并将结果保存到 DataFrame
def get_depth_for_chrom(bam_file, chrom):
    depth_file = f"/dev/shm/depth_{chrom}.txt"  # 为每个染色体使用不同的临时文件
    filtered_bam_file = f"/tmp/filtered_{chrom}.bam"  # 临时 BAM 文件
    # 过滤特定染色体并计算深度
    subprocess.run(["sambamba", "view", "-f", "bam", bam_file, chrom, "-o", filtered_bam_file])
    subprocess.run(["sambamba", "depth", "base", filtered_bam_file, "-o", depth_file])
    # 读取深度结果
    depth_df = pd.read_csv(depth_file, sep='\t', header=0, names=['chrom', 'position', 'ref', 'depth'])
    # 删除临时文件
    os.remove(filtered_bam_file)
    os.remove(depth_file)
    return depth_df


def get_depth(bam_file, num_threads=24):  # 添加 num_threads 参数
    with ThreadPoolExecutor(max_workers=num_threads) as executor:  # 设置最大线程数
        # 使用多线程处理染色体
        results = list(executor.map(lambda chrom: get_depth_for_chrom(bam_file, chrom), chromosomes))
    # 合并所有染色体的深度数据
    return pd.concat(results)

depth1 = get_depth(first_bam)
depth2 = get_depth(second_bam)

# 计算两个文件中的最大深度
max_depth = max(depth1['depth'].max(), depth2['depth'].max())

# 归一化并计算加权深度
merged_depth = pd.merge(
    depth1[['chrom', 'position', 'depth']],
    depth2[['chrom', 'position', 'depth']],
    on=['chrom', 'position'],
    how='outer',
    suffixes=('_1', '_2')
).fillna(0)
n=1
m=2
# 直接计算加权深度
merged_depth['weighted_depth_1'] = (merged_depth['depth_1'] / max_depth) * n
merged_depth['weighted_depth_2'] = (merged_depth['depth_2'] / max_depth) * m

# 计算总深度
merged_depth['total_depth'] = merged_depth['weighted_depth_1'] + merged_depth['weighted_depth_2']



# 输出到文件，只包含总深度
output_file = "/data/haocheng/data/bam/GM/merged_normalized_depth.txt"
merged_depth[['chrom', 'position', 'total_depth']].to_csv(output_file, sep='\t', header=True, index=False)
