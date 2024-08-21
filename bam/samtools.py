import os
import subprocess

# 获取目录下的所有 BAM 文件
bam_dir = "/data/haocheng/data/bam/GM/"
bam_files = [f for f in os.listdir(bam_dir) if f.endswith('.bam')]
bam_files = sorted(bam_files)  # 确保文件按字母顺序排序

# 初始化 m 和 n
m = 2
n = 1

# 存储归一化后的结果
result = ""

# 函数：计算 BAM 文件的覆盖度
def get_coverage(bam_file):
    cmd = f"samtools depth {bam_file} | awk '($1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/)'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout


# 使用 awk 进行归一化
def normalize_with_awk(coverage1, coverage2, norm_max, m, n):
    awk_script = f"""
    BEGIN {{
        OFS = "\\t"  # 设置输出字段分隔符为制表符
        norm_max = {norm_max}
        m = {m}
        n = {n}
        max_cov1 = 0
        max_cov2 = 0
    }}
    {{
        split($0, fields, OFS)
        key = fields[1] OFS fields[2] OFS fields[3]
        if (NR <= lines_in_coverage1) {{
            cov1[key] = fields[4]
            if (fields[4] > max_cov1) max_cov1 = fields[4]
        }} else {{
            cov2[key] = fields[4]
            if (fields[4] > max_cov2) max_cov2 = fields[4]
        }}
    }}
    END {{
        for (key in cov1) {{
            if (!(key in cov2)) cov2[key] = 0
        }}
        for (key in cov2) {{
            if (!(key in cov1)) cov1[key] = 0
        }}
        norm_max = (max_cov1 > max_cov2) ? max_cov1 : max_cov2
        for (key in cov1) {{
            norm1 = (cov1[key] / max_cov1) * norm_max
            norm2 = (cov2[key] / max_cov2) * norm_max
            norm_cov = (norm1 * (n/m)) + (norm2 * (1/m))
            total_cov = cov1[key] + cov2[key]
            print key, total_cov, norm_cov
        }}
    }}
    """
    # 合并两个覆盖度数据
    combined_coverage = coverage1 + "\n" + coverage2
    # 计算第一部分覆盖度数据的行数
    lines_in_coverage1 = coverage1.strip().count("\n") + 1
    awk_script = f"""lines_in_coverage1={lines_in_coverage1}; {awk_script}"""
    process = subprocess.Popen(
        ['awk', awk_script],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    result, error = process.communicate(input=combined_coverage)
    if process.returncode != 0:
        print(f"Error running awk: {error}")
    return result


def normalize_with_awk(coverage1, coverage2, norm_max, m, n):
    awk_script = f"""
    BEGIN {{
        OFS = "\\t"  # 设置输出字段分隔符为制表符
        norm_max = {norm_max}
        m = {m}
        n = {n}
        max_cov1 = 0
        max_cov2 = 0
    }}
    {{
        split($0, fields, OFS)
        key = fields[1] OFS fields[2] OFS fields[3]
        if (NR <= lines_in_coverage1) {{
            cov1[key] = fields[4]
            if (fields[4] > max_cov1) max_cov1 = fields[4]
        }} else {{
            cov2[key] = fields[4]
            if (fields[4] > max_cov2) max_cov2 = fields[4]
        }}
    }}
    END {{
        for (key in cov1) {{
            if (!(key in cov2)) cov2[key] = 0
        }}
        for (key in cov2) {{
            if (!(key in cov1)) cov1[key] = 0
        }}
        norm_max = (max_cov1 > max_cov2) ? max_cov1 : max_cov2
        for (key in cov1) {{
            norm1 = (cov1[key] / max_cov1) * norm_max
            norm2 = (cov2[key] / max_cov2) * norm_max
            norm_cov = (norm1 * (n/m)) + (norm2 * (1/m))
            total_cov = cov1[key] + cov2[key]
            print key, total_cov, norm_cov
        }}
    }}
    """
    # 计算第一部分覆盖度数据的行数
    lines_in_coverage1 = coverage1.strip().count("\n") + 1
    awk_script = f"""lines_in_coverage1={lines_in_coverage1}; {awk_script}"""
    # 调用 awk 脚本进行处理
    process = subprocess.Popen(
        ['awk', awk_script],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    result, error = process.communicate(input=coverage1 + "\n" + coverage2)
    if process.returncode != 0:
        print(f"Error running awk: {error}")
    return result
# 计算初始覆盖度
first_bam = os.path.join(bam_dir, bam_files[0])
second_bam = os.path.join(bam_dir, bam_files[1])

print(f"初次计算覆盖度，使用 BAM 文件: {first_bam} 和 {second_bam}")

coverage1 = get_coverage(first_bam)
coverage2 = get_coverage(second_bam)

# 归一化初始覆盖度
print("计算初次覆盖度归一化...")
result = normalize_with_awk(coverage1, coverage2, 1000, m, n)

# 循环处理剩余的 BAM 文件
bam_files = bam_files[2:]  # 删除前两个元素

while bam_files:
    next_bam = os.path.join(bam_dir, bam_files[0])
    print(f"现在正在处理 BAM 文件: {next_bam}")

    # 计算新的覆盖度
    coverage2 = get_coverage(next_bam)

    # 归一化并处理结果
    print("归一化覆盖度...")
    result = normalize_with_awk(result, coverage2, 1000, m, n)

    # 从数组中移除已使用的 BAM 文件
    bam_files = bam_files[1:]  # 删除第一个元素
    m += 1
    n += 1

# 将最终结果写入文件
output_file = "/data/haocheng/data/bed/GM/GM_normalized_coverage.bedgraph"

print(f"所有 BAM 文件处理完成，结果将输出到 {output_file}")

if result:
    with open(output_file, 'w') as f:
        f.write(result)
    print("结果已成功写入文件。")
else:
    print("结果为空，未写入文件。")
