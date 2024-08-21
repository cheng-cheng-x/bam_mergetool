import os
import subprocess
# 定义要执行的命令
command = "aria2c -s 8 -x 16 -i /data/haocheng/data/ACAT/txt/Others.txt -d /data/haocheng/data/bam/Others"
# 使用 subprocess 运行命令
try:
    subprocess.run(command, shell=True, check=True)
    print("Download completed successfully.")
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")

# 设置BAM文件目录
bam_dir = "/data/haocheng/data/bam/Others/"
bam_files = sorted([f for f in os.listdir(bam_dir) if f.endswith('.bam')])
# 获取第一个和第二个BAM文件的路径
bam_file1 = os.path.join(bam_dir, bam_files[0])
bam_file2 = os.path.join(bam_dir, bam_files[1])
sample1_peak_file = f"{bam_dir}/sample1_peaks.narrowPeak"
sample2_peak_file = f"{bam_dir}/sample2_peaks.narrowPeak"
sample1_signal_file = f"{bam_dir}/sample1_signal.bed"
sample2_signal_file = f"{bam_dir}/sample2_signal.bed"
combined_signal_file = f"{bam_dir}/combined_signals.bed"
output_file = f"{bam_dir}/updated_combined_signals.bed"
# 检查每个 BAM 文件是否存在索引，如果没有，则添加索引
for bam_file in bam_files:
    bam_file_path = os.path.join(bam_dir, bam_file)
    index_file_path = bam_file_path + ".bai"  # 索引文件的路径

    # 检查索引文件是否存在
    if not os.path.exists(index_file_path):
        print(f"Index file for {bam_file} not found. Creating index...")
        # 使用 samtools 创建索引
        subprocess.run(f"samtools index {bam_file_path}", shell=True, check=True)
        print(f"Index file created for {bam_file}.")
    else:
        print(f"Index file for {bam_file} already exists.")
# 运行MACS2以获得峰值
macs2_command1 = f"macs2 callpeak -t {bam_file1} -f BAM -g 3.1e9 -n sample1 --outdir {bam_dir} --keep-dup all"
macs2_command2 = f"macs2 callpeak -t {bam_file2} -f BAM -g 3.1e9 -n sample2 --outdir {bam_dir} --keep-dup all"
subprocess.run(macs2_command1, shell=True, check=True)
subprocess.run(macs2_command2, shell=True, check=True)
# 提取必要列并保存到新的文件中
awk_command1 = f"awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$10}}' {sample1_peak_file} > {sample1_signal_file}"
awk_command2 = f"awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$10}}' {sample2_peak_file} > {sample2_signal_file}"
# 运行 awk 命令
subprocess.run(awk_command1, shell=True, check=True)
subprocess.run(awk_command2, shell=True, check=True)
# 对输入文件进行排序
sort_sample1_command = f"bedtools sort -i {sample1_signal_file} > {sample1_signal_file}.sorted"
sort_sample2_command = f"bedtools sort -i {sample2_signal_file} > {sample2_signal_file}.sorted"
# 运行排序命令
subprocess.run(sort_sample1_command, shell=True, check=True)
subprocess.run(sort_sample2_command, shell=True, check=True)
# 使用 bedtools unionbedg 进行信号叠加，并将空白地方设置为0
unionbedg_command = f"bedtools unionbedg -i {sample1_signal_file}.sorted {sample2_signal_file}.sorted -filler 0 > {combined_signal_file}"
subprocess.run(unionbedg_command, shell=True, check=True)
print(f"Combined signals written to {combined_signal_file}")
output_file = f"{bam_dir}/updated_combined_signals.bed"
# awk 命令
awk_command = f"awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"($4 + $5)}}' {combined_signal_file} > {sample1_signal_file}"
# 调用 awk 命令
result =subprocess.run(awk_command, shell=True, check=True)
print(f"Processed data written to {sample1_signal_file}")
# 从数组中移除已使用的 BAM 文件
bam_files = bam_files[2:]  # 删除前两个元素
# 删除中间文件，只保留最终结果文件
# 处理剩余的 BAM 文件
while bam_files:
    next_bam = os.path.join(bam_dir, bam_files[0])  # 获取下一个 BAM 文件的路径
    print(f"现在正在处理数据: {next_bam}")
    sample2_peak_file = f"{bam_dir}/sample2_peaks.narrowPeak"
    sample2_signal_file = f"{bam_dir}/sample2_signal.bed"
    combined_signal_file = f"{bam_dir}/combined_signals.bed"
    output_file = f"{bam_dir}/updated_combined_signals.bed"
    print(f"Running MACS2 on {next_bam}...")
    macs2_command2 = f"macs2 callpeak -t {next_bam} -f BAM -g 3.1e9 -n sample2 --outdir {bam_dir} --keep-dup all"
    subprocess.run(macs2_command2, shell=True, check=True)
    print(f"Extracting necessary columns from {sample2_peak_file} to {sample2_signal_file}...")
    awk_command2 = f"awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$10}}' {sample2_peak_file} > {sample2_signal_file}"
    subprocess.run(awk_command2, shell=True, check=True)
    # 对输入文件进行排序
    sort_sample1_command = f"bedtools sort -i {sample1_signal_file} > {sample1_signal_file}.sorted"
    sort_sample2_command = f"bedtools sort -i {sample2_signal_file} > {sample2_signal_file}.sorted"
    # 运行排序命令
    subprocess.run(sort_sample1_command, shell=True, check=True)
    subprocess.run(sort_sample2_command, shell=True, check=True)
    # 使用 bedtools unionbedg 进行信号叠加，并将空白地方设置为0
    unionbedg_command = f"bedtools unionbedg -i {sample1_signal_file}.sorted {sample2_signal_file}.sorted -filler 0 > {combined_signal_file}"
    subprocess.run(unionbedg_command, shell=True, check=True)
    print(f"Combined signals written to {combined_signal_file}")
    print(f"Processing combined signals in {combined_signal_file} to {sample1_signal_file}...")
    awk_command = f"awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"($4 + $5)}}' {combined_signal_file} > {sample1_signal_file}"
    result = subprocess.run(awk_command, shell=True, check=True)
    print(f"Processed data written to {sample1_signal_file}")
    # 从数组中移除已使用的 BAM 文件
    bam_files = bam_files[1:]  # 删除第一个元素