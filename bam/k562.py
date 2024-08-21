import subprocess
import os
import tempfile
# 获取目录下的所有 BAM 文件
bam_dir = "/data/haocheng/data/bam/k562/"
bam_files = [f for f in os.listdir(bam_dir) if f.endswith('.bam')]
bam_files = sorted(bam_files)  # 确保文件按字母顺序排序

# 初始化 m 和 n
m = 2
n = 1

# 存储归一化后的结果prin
result = ""

# 初始覆盖度计算，获取前两个 BAM 文件
def get_coverage(bam_file):
    cmd = ["bedtools", "genomecov", "-ibam", bam_file, "-bg"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout

first_bam = os.path.join(bam_dir, bam_files[0])
second_bam = os.path.join(bam_dir, bam_files[1])

print(f"初次计算覆盖度，使用 BAM 文件: {first_bam} 和 {second_bam}")

coverage1 = get_coverage(first_bam)
coverage2 = get_coverage(second_bam)

# 将覆盖度数据写入临时文件
with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed') as temp_file1:
    temp_file1.write(coverage1)
    temp_file1_path = temp_file1.name

with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed') as temp_file2:
    temp_file2.write(coverage2)
    temp_file2_path = temp_file2.name

# 合并覆盖度
cmd = ["bedtools", "unionbedg", "-i", temp_file1_path, temp_file2_path]
merged_coverage_result = subprocess.run(cmd, capture_output=True, text=True)
if merged_coverage_result.returncode != 0:
    print(f"Error running bedtools unionbedg: {merged_coverage_result.stderr}")
    merged_coverage = ""
else:
    merged_coverage = merged_coverage_result.stdout

# 删除临时文件
os.remove(temp_file1_path)
os.remove(temp_file2_path)
#print(merged_coverage)


# 使用 awk 进行归一化
def normalize_with_awk(coverage_data, norm_max, m, n):
    awk_script = f"""
    BEGIN {{
        OFS = "\\t"  # 设置输出字段分隔符为制表符
        norm_max = {norm_max}
        m = {m}
        n = {n}
        max_cov4 = 0
        max_cov5 = 0
    }}
    {{
        if ($4 > max_cov4) max_cov4 = $4
        if ($5 > max_cov5) max_cov5 = $5
        data[NR] = $0
    }}
    END {{
        norm_max = (max_cov4 > max_cov5) ? max_cov4 : max_cov5
        if (norm_max == max_cov4) {{
            for (i = 1; i <= NR; i++) {{
                split(data[i], fields, OFS)
                norm4 = fields[4]
                norm5 = (fields[5] / max_cov5) * norm_max
                norm_cov = (norm4 * (n/m)) + (norm5 * (1/m))
                print fields[1] "\t" fields[2] "\t" fields[3] "\t" norm_cov
            }}
        }} else {{
            for (i = 1; i <= NR; i++) {{
                split(data[i], fields, OFS)
                norm4 = (fields[4] / max_cov4) * norm_max
                norm5 = fields[5]
                norm_cov = (norm4 * (n/m)) + (norm5 * (1/m))
                print fields[1] "\t" fields[2] "\t" fields[3] "\t" norm_cov
            }}
        }}
    }}
    """
    process = subprocess.Popen(
        ['awk', awk_script],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    result = process.communicate(input=coverage_data)[0]
    return result


print("计算初次覆盖度归一化...")
result = normalize_with_awk(merged_coverage, 1000, m, n)
# 从数组中移除已使用的 BAM 文件
bam_files = bam_files[2:]  # 删除前两个元素

# 自增 m
m += 1
n += 1


# 将最终结果写入文件
output_file = "/data/haocheng/data/bed/k562/k562_normalized_coverage.bedgraph"
# 打印输出文件的路径
print(f"所有 BAM 文件处理完成，结果将输出到 {output_file}")
# 检查 result 是否为空
if result:
    # 使用 with 语句打开文件，以 'w' 模式写入
    with open(output_file, 'w') as f:
        f.write(result)
    print("结果已成功写入文件。")
else:
    print("结果为空，未写入文件。")




#result = result.replace(" ", "\t")
# 循环处理剩余的 BAM 文件
while bam_files:
    result = result.replace(" ", "\t")
    next_bam = os.path.join(bam_dir,bam_files[0])
    print(f"现在正在处理第 {m} 个数据: {next_bam}")
    # 使用当前的归一化覆盖度数据
    coverage1 = result
    # 计算新的覆盖度并存储prin
    print(f"计算 BAM 文件: {next_bam} 的覆盖度...")
    coverage2 = get_coverage(next_bam)
    # 将覆盖度数据写入临时文件
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed') as temp_file1:
     temp_file1.write(coverage1)
     temp_file1_path = temp_file1.name
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed') as temp_file2:
     temp_file2.write(coverage2)
     temp_file2_path = temp_file2.name
    # 合并覆盖度
    merged_coverage =""
    merged_coverage_result=None
    cmd = ["bedtools", "unionbedg", "-i", temp_file1_path, temp_file2_path]
    merged_coverage_result = subprocess.run(cmd, capture_output=True, text=True)
    merged_coverage = merged_coverage_result.stdout
    # 删除临时文件
    os.remove(temp_file1_path)
    os.remove(temp_file2_path)
    #print(merged_coverage)
    # 归一化并处理结果，直接累加到 result
    print("归一化覆盖度...")
    result = normalize_with_awk(merged_coverage, 1000, m, n)
    # 从数组中移除已使用的 BAM 文件
    bam_files = bam_files[1:]  # 删除第一个元素
    # 检查 bam_files 的状态
    print(f"剩余 BAM 文件列表: {bam_files}")
    # 自增 m
    m += 1
    n += 1
    # 打印处理进度
    # 将最终结果写入文件
    output_file = "/data/haocheng/data/bed/k562/k562_normalized_coverage.bedgraph"
   # 打印输出文件的路径
    #print(f"所有 BAM 文件处理完成，结果将输出到 {output_file}")
    # 检查 result 是否为空
    if result:
        # 使用 with 语句打开文件，以 'w' 模式写入
        with open(output_file, 'w') as f:
           f.write(result)
        print("结果已成功写入文件。")
    else:
        print("结果为空，未写入文件。")
    print(f"已处理 {m} 个 BAM 文件，剩余 {len(bam_files)} 个文件待处理...")

# 假设 get_coverage 和 normalize_with_awk 函数已定义
# 将最终结果写入文件
output_file = "/data/haocheng/data/bed/k562/k562_normalized_coverage.bedgraph"

# 打印输出文件的路径
print(f"所有 BAM 文件处理完成，结果将输出到 {output_file}")

# 检查 result 是否为空
if result:
    # 使用 with 语句打开文件，以 'w' 模式写入
    with open(output_file, 'w') as f:
        f.write(result)
    print("结果已成功写入文件。")
else:
    print("结果为空，未写入文件。")


'''
max_cov1 = max(float(line.split()[3]) for line in coverage1.strip().split('\n'))
max_cov2 = max(float(line.split()[3]) for line in coverage2.strip().split('\n'))
print(f"最大覆盖度值：coverage1: {max_cov1}, coverage2: {max_cov2}")
'''
'''
with open(output_file, 'r') as f:
     result = f.read()  # 将文件内容读取到 result
     '''