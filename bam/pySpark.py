from pyspark.sql import SparkSession
from pyspark.sql import functions as F
import subprocess
import os

# 创建 SparkSession
spark = SparkSession.builder \
    .appName("Normalization") \
    .getOrCreate()

# 函数：计算 BAM 文件的覆盖度
def get_coverage(bam_file):
    cmd = ["samtools", "depth", bam_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout

# 获取目录下的所有 BAM 文件
bam_dir = "/data/haocheng/data/bam/GM/"
bam_files = [f for f in os.listdir(bam_dir) if f.endswith('.bam')]
bam_files = sorted(bam_files)  # 确保文件按字母顺序排序

# 初始化 m 和 n
m = 2
n = 1
bam_file1 = os.path.join(bam_dir, bam_files[0])
bam_file2 = os.path.join(bam_dir, bam_files[1])

coverage1 = get_coverage(bam_file1)
coverage2 = get_coverage(bam_file2)

# 将覆盖度数据转换为 DataFrame
def create_df(coverage_data):
    data = [line.split() for line in coverage_data.strip().split("\n")]
    return spark.createDataFrame(data, schema=["chrom", "position", "coverage"]).withColumn("coverage", F.col("coverage").cast("int"))

coverage1_df = create_df(coverage1)
coverage2_df = create_df(coverage2)

# 合并两个 DataFrame
combined_df = coverage1_df.join(coverage2_df, on=["chrom", "position"], how="outer").fillna(0)
# 计算最大覆盖度
max_cov1 = combined_df.agg(F.max("coverage")).first()[0]
max_cov2 = combined_df.agg(F.max("coverage")).first()[0]
# 归一化
norm_max = max(max_cov1, max_cov2)
m = 2  # 根据你的需求设置
n = 1  # 根据你的需求设置
normalized_df = combined_df.withColumn(
    "normalized_coverage",
    (F.col("coverage1") / max_cov1 * norm_max * (n/m)) +
    (F.col("coverage2") / max_cov2 * norm_max * (1/m))
)
# 收集结果到本地
result = normalized_df.select("chrom", "position", "normalized_coverage").rdd.map(lambda row: f"{row.chrom}\t{row.position}\t{row.normalized_coverage}").collect()
# 输出结果
normalized_df.select("chrom", "position", "normalized_coverage") \
    .write.csv("path/to/output/normalized_coverage.csv", sep="\t", header=True)

# 关闭 SparkSession
spark.stop()



while bam_files:
    next_bam = os.path.join(bam_dir, bam_files[0])
    print(f"现在正在处理 BAM 文件: {next_bam}")

    # 计算新的覆盖度
    coverage2 = get_coverage(next_bam)

    # 归一化并处理结果
    print("归一化覆盖度...")
    result = normalize_with_spark(result, coverage2, 1000, m, n)

    # 从数组中移除已使用的 BAM 文件
    bam_files = bam_files[1:]  # 删除第一个元素
    m += 1
    n += 1

# 将最终结果写入文件
output_file = "/data/haocheng/data/bed/GM/GM_normalized_coverage.bedgraph"

print(f"所有 BAM 文件处理完成，结果将输出到 {output_file}")

if result:
    with open(output_file, 'w') as f:
        f.write("\n".join(result))
    print("结果已成功写入文件。")
else:
    print("结果为空，未写入文件。")