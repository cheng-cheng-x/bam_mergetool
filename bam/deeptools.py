import subprocess
import pysam
import shutil
import os

# 设置BAM文件目录
bam_dir = "/data/haocheng/data/bam/GM/"
bam_files = sorted([f for f in os.listdir(bam_dir) if f.endswith('.bam')])

# 获取第一个和第二个BAM文件的路径
bam_file1 = os.path.join(bam_dir, bam_files[0])
bam_file2 = os.path.join(bam_dir, bam_files[1])

# 定义输出BigWig文件的路径
bw_file1 = bam_file1.replace('.bam', '.bw')
bw_file2 = bam_file2.replace('.bam', '.bw')
merged_bw_file = os.path.join(bam_dir, 'merged_output.bw')
temp_file = '/data/haocheng/data/bam/GM/temp.bw'

def get_chromosome_lengths(bam_file):
    # 使用 pysam 获取 BAM 文件头信息
    samfile = pysam.AlignmentFile(bam_file, "rb")
    chrom_lengths = {chrom: samfile.get_reference_length(chrom) for chrom in samfile.references if chrom.startswith('chr') and chrom != 'chrM'}
    samfile.close()
    return chrom_lengths

def convert_bam_to_bw(bam_file, bw_file):
    chrom_lengths = get_chromosome_lengths(bam_file)

    # 使用 chr1 到 chrY 的长度进行转换
    regions = [f"{chrom}:1:{length}" for chrom, length in chrom_lengths.items() if chrom in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]]

    for region in regions:
        command = [
            "bamCoverage",
            "-b", bam_file,
            "-o", bw_file,
            "--binSize", "1",
            "--region", region
        ]
        print(f"正在转换区域: {region}...")
        subprocess.run(command, check=True)
        print(f"{region} 转换完成。")

# 转换两个BAM文件
print(f"Converting {bam_file1} to {bw_file1}...")
convert_bam_to_bw(bam_file1, bw_file1)
print(f"Conversion complete: {bw_file1}")

print(f"Converting {bam_file2} to {bw_file2}...")
convert_bam_to_bw(bam_file2, bw_file2)
print(f"Conversion complete: {bw_file2}")

# 使用bigwigMerge叠加BigWig文件
def merge_bigwig_files(bw_file1, bw_file2, merged_bw_file):
    command = [
        'bigwigMerge',
        bw_file1,
        bw_file2,
        merged_bw_file
    ]
    subprocess.run(command, check=True)
    shutil.copy(merged_bw_file, temp_file)

# 叠加BigWig文件
print(f"Merging {bw_file1} and {bw_file2} into {merged_bw_file}...")
merge_bigwig_files(bw_file1, bw_file2, merged_bw_file)
print(f"Merged BigWig file saved as '{merged_bw_file}'.")

# 从数组中移除已使用的 BAM 文件
bam_files = bam_files[2:]

while bam_files:
    next_bam = os.path.join(bam_dir, bam_files[0])
    next_bw_file = next_bam.replace('.bam', '.bw')  # 新的 BigWig 文件名
    print(f"Converting {next_bam} to {next_bw_file}...")
    convert_bam_to_bw(next_bam, next_bw_file)
    print(f"Conversion complete: {next_bw_file}")

    print(f"Merging {next_bw_file} with {temp_file} into {merged_bw_file}...")
    merge_bigwig_files(next_bw_file, temp_file, merged_bw_file)
    print(f"Merged BigWig file saved as '{merged_bw_file}'.")

    bam_files = bam_files[1:]  # 删除第一个元素
