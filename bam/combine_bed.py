import os
import subprocess
# 设置BAM文件目录
bed_dir = "/data/haocheng/data/bam/result/"
bed_work="/data/haocheng/data/bam/result/work"
bed_files = sorted([f for f in os.listdir(bed_dir) if f.endswith('.bed')])
# 获取第一个和第二个BAM文件的路径
bed_file1 = os.path.join(bed_dir, bed_files[0])
bed_file2 = os.path.join(bed_dir, bed_files[1])
# 定义排序后的输出文件路径，使用固定名称
sorted_bed_file1 = os.path.join(bed_work, "sorted_file1.bed")
sorted_bed_file2 = os.path.join(bed_work, "sorted_file2.bed")
# 使用bedtools sort命令对第一个BED文件排序，并将结果重定向到输出文件
with open(sorted_bed_file1, 'w') as out1:
    subprocess.run(["bedtools", "sort", "-i", bed_file1], stdout=out1)
# 使用bedtools sort命令对第二个BED文件排序，并将结果重定向到输出文件
with open(sorted_bed_file2, 'w') as out2:
    subprocess.run(["bedtools", "sort", "-i", bed_file2], stdout=out2)
print(f"Sorted files saved as {sorted_bed_file1} and {sorted_bed_file2}")
# 定义unionbedg的输出文件路径
union_bed_file = os.path.join(bed_work, "unionbedg_output.bed")
# 使用bedtools unionbedg命令对两个排序后的BED文件进行叠加，并将结果重定向到输出文件
with open(union_bed_file, 'w') as out:
    subprocess.run([
        "bedtools", "unionbedg",
        "-i", sorted_bed_file1, sorted_bed_file2,
        # "-header",  # 如果需要保留头部信息可以添加这个选项
    ], stdout=out)
print(f"Unionbedg output saved as {union_bed_file}")
# 定义最终处理后的输出文件路径
final_output_file = os.path.join(bed_work, "final_output.bed")
# 使用awk进行操作，将第1-3列保留，并将第4列和第5列相加，输出到新的文件中
awk_command = f"awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"($4 + $5)}}' {union_bed_file} > {final_output_file}"
subprocess.run(awk_command, shell=True)
print(f"Final output saved as {final_output_file}")

bed_files = bed_files[2:]  
while bed_files:
    next_bed = os.path.join(bed_dir, bed_files[0])  # 获取下一个 BAM 文件的路径
    print(f"现在正在处理数据: {next_bed}")
    with open(sorted_bed_file1, 'w') as out1:
     subprocess.run(["bedtools", "sort", "-i", final_output_file], stdout=out1)
    with open(sorted_bed_file2, 'w') as out2:
     subprocess.run(["bedtools", "sort", "-i", next_bed], stdout=out2)
    print(f"Sorted files saved as {sorted_bed_file1} and {sorted_bed_file2}")
    # 定义unionbedg的输出文件路径
    union_bed_file = os.path.join(bed_work, "unionbedg_output.bed")
    # 使用bedtools unionbedg命令对两个排序后的BED文件进行叠加，并将结果重定向到输出文件
    with open(union_bed_file, 'w') as out:
        subprocess.run([
            "bedtools", "unionbedg",
            "-i", sorted_bed_file1, sorted_bed_file2,
            # "-header",  # 如果需要保留头部信息可以添加这个选项
        ], stdout=out)
    print(f"Unionbedg output saved as {union_bed_file}")
    # 定义最终处理后的输出文件路径
    final_output_file = os.path.join(bed_work, "final_output.bed")
    # 使用awk进行操作，将第1-3列保留，并将第4列和第5列相加，输出到新的文件中
    awk_command = f"awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"($4 + $5)}}' {union_bed_file} > {final_output_file}"
    subprocess.run(awk_command, shell=True)
    print(f"Final output saved as {final_output_file}")
    bam_files = bam_files[1:]  # 删除第一个元素
