print(f"Size of result object: {sys.getsizeof(result)} bytes")
print(f"Size of result object: {sys.getsizeof(merged_coverage)} bytes")
coverage1
print(f"Size of result object: {sys.getsizeof(coverage1)} bytes")



aria2c -d /data/haocheng/data/DNA -o Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz

gunzip /data/haocheng/data/DNA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

aria2c -d /data/haocheng/data/DNA -o Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip /data/haocheng/data/DNA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz









# 获取目录下的所有 BAM 文件
bam_files=("/data/haocheng/data/bam/k562/"*.bam)


# 初始化 m 和 n
m=2
n=1

# 存储归一化后的结果
result=""

# 初始覆盖度计算，获取前两个 BAM 文件
first_bam="${bam_files[0]}"
second_bam="${bam_files[1]}"

echo "初次计算覆盖度，使用 BAM 文件: $first_bam 和 $second_bam"

# 初次计算覆盖度并存储结果到数组
coverage1=$(bedtools genomecov -ibam "$first_bam" -bg)
coverage2=$(bedtools genomecov -ibam "$second_bam" -bg)

# 合并覆盖度
merged_coverage=$(echo "$coverage1" && echo "$coverage2" | bedtools unionbedg -i -)

# 归一化并处理结果
echo "计算初次覆盖度归一化..."
result=$(awk -v norm_max=1000 -v m="$m" -v n="$n" '
{
    if ($4 > max_cov4) max_cov4 = $4;  # 计算第四列最大值
    if ($5 > max_cov5) max_cov5 = $5;  # 计算第五列最大值
    data[NR] = $0;                      # 记录每一行
}
END {
    for (i = 1; i <= NR; i++) {
        split(data[i], fields, OFS);   # 分割存储的行
        norm4 = (fields[4] / max_cov4) * norm_max;  # 归一化第四列
        norm5 = (fields[5] / max_cov5) * norm_max;  # 归一化第五列
        print fields[1], fields[2], fields[3], (norm4 * (n/m)) + (norm5 * (1/m));  # 输出结果
    }
}' <<< "$merged_coverage")

# 清理临时变量
unset coverage1 coverage2 merged_coverage

# 从数组中移除已使用的 BAM 文件
bam_files=("${bam_files[@]:2}")  # 删除前两个元素
# 自增 m
m=$((m + 1))
n=$((n + 1))

# 循环处理剩余的 BAM 文件
while [ ${#bam_files[@]} -gt 0 ]; do
    # 读取数组中的下一个 BAM 文件
    next_bam="${bam_files[0]}"

    echo "现在正在处理第 $m 个数据: $next_bam"

    # 使用当前的归一化覆盖度数据
    coverage1="$result"
     result=""  # 清空 result

    # 计算新的覆盖度并存储
    echo "计算 BAM 文件: $next_bam 的覆盖度..."
    coverage2=$(bedtools genomecov -ibam "$next_bam" -bg)

    # 合并覆盖度数据
    merged_coverage=$(echo "$coverage1" && echo "$coverage2" | bedtools unionbedg -i -)

    # 归一化并处理结果，直接累加到 result
    echo "归一化覆盖度..."
    result=$(awk -v norm_max=1000 -v m="$m" -v n="$n" '
    {
        if ($4 > max_cov4) max_cov4 = $4;  # 计算第四列最大值
        if ($5 > max_cov5) max_cov5 = $5;  # 计算第五列最大值
        data[NR] = $0;                      # 记录每一行
    }
    END {
        for (i = 1; i <= NR; i++) {
            split(data[i], fields, OFS);   # 分割存储的行
            norm4 = (fields[4] / max_cov4) * norm_max;  # 归一化第四列
            norm5 = (fields[5] / max_cov5) * norm_max;  # 归一化第五列
            print fields[1], fields[2], fields[3], (norm4 * (n/m)) + (norm5 * (1/m));  # 输出结果
        }
    }' <<< "$merged_coverage")

    # 清理临时变量
    unset coverage2 merged_coverage

    # 从数组中移除已使用的 BAM 文件
    bam_files=("${bam_files[@]:1}")  # 删除第一个元素

    # 自增 m
    m=$((m + 1))
    n=$((n + 1))
    # 打印处理进度
    echo "已处理 $m 个 BAM 文件，剩余 ${#bam_files[@]} 个文件待处理..."
done

# 将最终结果写入文件
echo "所有 BAM 文件处理完成，结果输出到 /data/haocheng/data/bed/k562/k562_normalized_coverage.bedgraph"
echo "$result" > /data/haocheng/data/bed/k562/k562_normalized_coverage.bedgraph

echo "结果已成功写入文件."