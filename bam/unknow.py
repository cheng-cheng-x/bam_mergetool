import os
import subprocess
# 获取 BAM 文件目录
bam_dir = "/data/haocheng/data/bam/GM/"
bam_files = [f for f in os.listdir(bam_dir) if f.endswith('.bam')]
bam_files = sorted(bam_files)  # 确保文件按字母顺序排序
# 初始化 m 和 n 的值
m = 2  # 用于归一化的参数
n = 1  # 用于归一化的参数
bam_file1 = os.path.join(bam_dir, bam_files[0])  # 第一个 BAM 文件
bam_file2 = os.path.join(bam_dir, bam_files[1])  # 第二个 BAM 文件
# 函数：计算 BAM 文件的覆盖度
def get_coverage(bam_file):
    cmd = f"samtools depth {bam_file} | awk '($1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/)'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout
# 获取两个 BAM 文件的覆盖度
coverage1 = get_coverage(bam_file1)  # 获取第一个 BAM 文件的覆盖度
coverage2 = get_coverage(bam_file2)  # 获取第二个 BAM 文件的覆盖度

def normalize_with_awk(coverage1, coverage2, norm_max, m, n):
    awk_script = f"""
    BEGIN {{
        OFS = "\\t";  # 设置输出字段分隔符为制表符
        norm_max = {norm_max};
        m = {m};
        n = {n};
        max_cov1 = 0;
        max_cov2 = 0;
    }}
    {{
        split($0, fields, OFS);
        key = fields[1] OFS fields[2];  # 使用染色体和位点作为键
        
        if (NR <= lines_in_coverage1) {{
            cov1[key] = fields[3];  # 覆盖度是第三列
            if (fields[3] > max_cov1) max_cov1 = fields[3];
        }} else {{
            cov2[key] = fields[3];  # 覆盖度是第三列
            if (fields[3] > max_cov2) max_cov2 = fields[3];
        }}
    }}
    END {{
        for (key in cov1) {{
            if (!(key in cov2)) cov2[key] = 0;
        }}
        for (key in cov2) {{
            if (!(key in cov1)) cov1[key] = 0;
        }}
        norm_max = (max_cov1 > max_cov2) ? max_cov1 : max_cov2;
        for (key in cov1) {{
            norm1 = (cov1[key] / max_cov1) * norm_max;
            norm2 = (cov2[key] / max_cov2) * norm_max;
            norm_cov = (norm1 * (n/m)) + (norm2 * (1/m));
            total_cov = cov1[key] + cov2[key];
            split(key, parts, OFS);  # 将key分割为染色体和位点
            print parts[1], parts[2], norm_cov;  # 打印染色体、位点和最终值
        }}
    }}
    """
    # 计算第一部分覆盖度数据的行数
    lines_in_coverage1 = coverage1.strip().count("\n") + 1
    awk_script = f"""lines_in_coverage1={lines_in_coverage1}; {awk_script}"""
    # 合并覆盖度数据
    combined_coverage = coverage1 + "\n" + coverage2
    # 调用 awk 脚本进行处理
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


result = normalize_with_awk(coverage1, coverage2, norm_max=100, m=2, n=1)