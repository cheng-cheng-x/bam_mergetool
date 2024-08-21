aria2c -d /data/haocheng/data/DNA/-s 8 -x 16 http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
bigWigInfo /data/haocheng/k562.bigWig
bedGraphToBigWig /data/haocheng/data/bam/result/GM12878_no_chrEBV.bed hg38.chrom.sizes.1 /data/haocheng/data/bam/result/GM12878.bigwig
bedtools sort -i /data/haocheng/updated_combined_signals.bed > /data/haocheng/sorted_k562.bed
aria2c -s 8 -x 16 -i /data/haocheng/data/ACAT/txt/GM_3.txt -d /data/haocheng/data/bam/GM

grep -v "chrEBV" /data/haocheng/data/bam/GM12878/GM12878.bed > /data/haocheng/data/bam/result/GM12878_no_chrEBV.bed
awk '{if($4 > max) max=$4} END {print max}' /data/haocheng/data/bam/result/all_no_chrEBV_normalized10.bed


bedGraphToBigWig /data/haocheng/data/bam/result/all_no_chrEBV_normalized100.bed hg38.chrom.sizes.1 /data/haocheng/data/bam/result/All_100.bigwig

grep -v "chrEBV" /data/haocheng/data/bam/result/work/final_output.bed > /data/haocheng/data/bam/result/all_no_chrEBV.bed

aria2c -s 8 -x 16 -i /data/haocheng/data/ACAT/txt/k562_1.txt -d /data/haocheng/data/bam/k562












cd /data/haocheng/data/bed/k562
gunzip ENCFF842UZU.bed.gz
head ENCFF842UZU.bed
awk '{if($7>max) max=$5} END {print max}' ENCFF842UZU.bed
https://www.encodeproject.org/files/ENCFF778GYX/@@download/ENCFF778GYX.bam
bedtools genomecov -ibam /data/haocheng/data/bam/k562/ENCFF879XNY.bam -bg > coverage.bedgraph




# 计算覆盖度
bedtools genomecov -ibam /data/haocheng/data/bam/k562/ENCFF879XNY.bam -bg > coverage1.bedgraph
bedtools genomecov -ibam /data/haocheng/data/bam/k562/ENCFF203LSX.bam -bg > coverage2.bedgraph

# 查看合并后的覆盖度文件的前几行内容
head coverage1.bedgraph
# 查看合并后的覆盖度文件的前几行内容
head coverage2.bedgraph

# 合并覆盖度文件
bedtools unionbedg -i coverage1.bedgraph coverage2.bedgraph > merged_coverage.bedgraph
less merged_coverage.bedgraph

# 计算每列的最大覆盖度值
max_coverage_col4=$(awk 'BEGIN{max=0} {if($4>max) max=$4} END {print max}' merged_coverage.bedgraph)
max_coverage_col5=$(awk 'BEGIN{max=0} {if($5>max) max=$5} END {print max}' merged_coverage.bedgraph)

# 设置归一化的最大值
normalized_max=1000

# 归一化每列的覆盖度值并相加，输出到新的文件
awk -v max_cov4="$max_coverage_col4" -v max_cov5="$max_coverage_col5" -v norm_max="$normalized_max" 'BEGIN{OFS="\t"} {norm4=($4 / max_cov4) * norm_max; norm5=($5 / max_cov5) * norm_max; print $1, $2, $3, norm4 + norm5}' merged_coverage.bedgraph > /data/haocheng/data/bed/k562/normalized_coverage.bedgraph

# 清理临时文件
rm coverage1.bedgraph coverage2.bedgraph merged_coverage.bedgraph


# 计算每列的最大覆盖度值
max_coverage_col4=$(awk 'BEGIN{max=0} {if($4>max) max=$4} END {print max}' merged_coverage.bedgraph)
max_coverage_col5=$(awk 'BEGIN{max=0} {if($5>max) max=$5} END {print max}' merged_coverage.bedgraph)

# 设置归一化的最大值
normalized_max=1000

# 归一化每列的覆盖度值，后续乘法处理
awk -v m=2 -v n=1 -v max_cov4="$max_coverage_col4" -v max_cov5="$max_coverage_col5" -v norm_max="$normalized_max" 'BEGIN{OFS="\t"} {
    norm4=($4 / max_cov4) * norm_max;   # 归一化第四列
    norm5=($5 / max_cov5) * norm_max;   # 归一化第五列
    
    # 按新的逻辑进行相加
    print $1, $2, $3, (norm4 * (n/m)) + (norm5 * (1/m));
}' merged_coverage.bedgraph > /data/haocheng/data/bed/k562/normalized_coverage.bedgraph




# 查看合并后的覆盖度文件的前几行内容
head coverage1.bedgraph
# 查看合并后的覆盖度文件的前几行内容
head coverage2.bedgraph





# 计算覆盖度
bedtools genomecov -ibam /data/haocheng/data/bam/k562/ENCFF879XNY.bam -bg > coverage1.bedgraph
bedtools genomecov -ibam /data/haocheng/data/bam/k562/ENCFF203LSX.bam -bg > coverage2.bedgraph



# 合并覆盖度文件
bedtools unionbedg -i coverage1.bedgraph coverage2.bedgraph > merged_coverage.bedgraph
less merged_coverage.bedgraph

awk -v norm_max=1000 -v m=2 -v n=1 '
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
}' merged_coverage.bedgraph > /data/haocheng/data/bed/k562/normalized_coverage.bedgraph


#!/bin/bash

# 获取目录下的所有 BAM 文件
bam_files=("/data/haocheng/data/bam/k562/"*.bam)

# 初始化 m 和 n
m=2
n=1

# 初始覆盖度计算，获取前两个 BAM 文件
first_bam="${bam_files[0]}"
second_bam="${bam_files[1]}"

echo "初次计算覆盖度，使用 BAM 文件: $first_bam 和 $second_bam"

# 初次计算覆盖度
bedtools genomecov -ibam "$first_bam" -bg > coverage1.bedgraph
bedtools genomecov -ibam "$second_bam" -bg > coverage2.bedgraph

# 合并覆盖度文件
bedtools unionbedg -i coverage1.bedgraph coverage2.bedgraph > merged_coverage.bedgraph

# 归一化并输出结果
awk -v norm_max=1000 -v m="$m" -v n="$n" '
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
}' merged_coverage.bedgraph > /data/haocheng/data/bed/k562/normalized_coverage.bedgraph

echo "初次覆盖度计算和归一化完成，结果输出到 /data/haocheng/data/bed/k562/normalized_coverage.bedgraph"

# 从数组中移除已使用的 BAM 文件
bam_files=("${bam_files[@]:2}")  # 删除前两个元素

# 循环处理剩余的 BAM 文件
while [ ${#bam_files[@]} -gt 0 ]; do
    # 读取数组中的下一个 BAM 文件
    next_bam="${bam_files[0]}"

    echo "现在正在处理第 $m 个数据: $next_bam"

    # 将上一次的 normalized_coverage.bedgraph 作为 coverage1.bedgraph
    cp /data/haocheng/data/bed/k562/normalized_coverage.bedgraph coverage1.bedgraph

    # 计算新的覆盖度
    bedtools genomecov -ibam "$next_bam" -bg > coverage2.bedgraph

    # 合并覆盖度文件
    bedtools unionbedg -i coverage1.bedgraph coverage2.bedgraph > merged_coverage.bedgraph

    # 归一化并输出结果
    awk -v norm_max=1000 -v m="$m" -v n="$n" '
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
    }' merged_coverage.bedgraph > /data/haocheng/data/bed/k562/normalized_coverage.bedgraph

    echo "覆盖度计算和归一化完成，结果输出到 /data/haocheng/data/bed/k562/normalized_coverage.bedgraph"

    # 保存当前数组长度
    array_length=${#bam_files[@]}

    # 从数组中移除已使用的 BAM 文件
    bam_files=("${bam_files[@]:1}")  # 删除第一个元素

    # 自增 m
    m=$((m + 1))

    # 检查是否还有 BAM 文件剩余
    if [ $array_length -eq 0 ]; then
        break
    fi
done

echo "所有 BAM 文件处理完成。"






aria2c -s 16 -x 16 -i /data/haocheng/data/ACAT/txt/k562_2.txt -d /data/haocheng/data/bam/k562
aria2c -s 16 -x 16 -i /data/haocheng/data/ACAT/txt/GM12878.txt -d /data/haocheng/data/bam/GM12878



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

    # 从数组中移除已使用的 BAM 文件
    bam_files=("${bam_files[@]:1}")  # 删除第一个元素

    # 自增 m
   m=$((m + 1))
    n=$((n + 1))
    # 打印处理进度
    echo "已处理 $m 个 BAM 文件，剩余 ${#bam_files[@]} 个文件待处理..."
done

# 将最终结果写入文件
echo "所有 BAM 文件处理完成，结果输出到 /data/haocheng/data/bed/k562/normalized_coverage.bedgraph"
echo "$result" > /data/haocheng/data/bed/k562/k562_normalized_coverage.bedgraph

echo "结果已成功写入文件."



echo "计算覆盖度..."
coverage2=$(printf "%s\n" "${bam_files[@]}" | parallel --no-notice -j 4 "bedtools genomecov -ibam {} -bg")










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
