#!/usr/bin/env bash

# ==============================================================================
# Script Name: run_homer.sh
# Description: HOMER 自动化峰注释脚本 (支持自定义输出目录 + R语言自动清洗归类)
# ==============================================================================

# 1. 帮助文档与格式说明
if [ -z "$1" ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo -e "\033[1;36m==================================================\033[0m"
    echo -e "\033[1;36m           HOMER 自动化峰注释与清洗脚本           \033[0m"
    echo -e "\033[1;36m==================================================\033[0m"
    echo -e "\033[1;33m用法 (Usage):\033[0m"
    echo "  ./run_homer.sh <输入文件.bed> [参考基因组, 默认mm10] [输出目录, 默认与输入文件同目录]"
    echo ""
    echo -e "\033[1;33m【输入 BED 文件格式严格要求】\033[0m"
    echo "请确保文件为 制表符 (Tab) 分隔："
    echo "  · 第 1 列 (chr)   : 染色体名称   [必须，例如: chr1]"
    echo "  · 第 2 列 (start) : 起始坐标     [必须，例如: 10000]"
    echo "  · 第 3 列 (end)   : 终止坐标     [必须，例如: 10500]"
    echo "  · 第 4 列 (ID)    : Peak 名字    [可选，若无脚本会自动生成]"
    echo "  · 第 5 列 (Strand): 正负链 (+/-) [可选，若无/为空则全按正链处理]"
    echo ""
    echo -e "\033[1;33m示例 (Examples):\033[0m"
    echo "  ./run_homer.sh my_peaks.bed"
    echo "  ./run_homer.sh my_peaks.bed mm10 /home1/out/"
    echo -e "\033[1;36m==================================================\033[0m"
    exit 1
fi

# 2. 自动解析变量
input_bed=$1
genome=${2:-mm10} 
out_dir=${3:-$(dirname "$input_bed")} 

filename=$(basename -- "$input_bed")
base="${filename%.bed}"

# 确保输出目录存在
if [ ! -d "$out_dir" ]; then
    echo "📁 正在创建输出目录: $out_dir"
    mkdir -p "$out_dir"
fi

# 自动推导输出文件名
homer_bed="${out_dir}/${base}_homer_input.bed"
out_csv="${out_dir}/${base}_peak_anno.csv"
out_cleaned="${out_dir}/${base}_peak_anno_cleaned.csv"
out_log="${out_dir}/${base}_anno.log"

# 获取当前 run_homer.sh 脚本所在的目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# 自动拼接 bed2homer.sh 的路径
bed2homer_script="$SCRIPT_DIR/bed2homer.sh"

echo "=== 启动 HOMER 自动化注释流程 ==="
echo "输入文件  : $input_bed"
echo "参考基因组: $genome"
echo "输出目录  : $out_dir"

# 3. 运行第一步：格式转换
echo "[1/3] 正在清洗并转换为 HOMER 专用格式..."
sh "$bed2homer_script" "$input_bed" "$homer_bed"

if [ $? -ne 0 ]; then
    echo "❌ 格式转换失败，请对照 -h 说明检查 BED 文件格式！"
    exit 1
fi

# 4. 运行第二步：HOMER 注释
echo "[2/3] 正在运行 annotatePeaks.pl (可能需要几分钟)..."
annotatePeaks.pl "$homer_bed" "$genome" 1> "$out_csv" 2> "$out_log"

if [ $? -ne 0 ]; then
    echo "❌ HOMER 注释失败，请查看日志寻找报错原因: $out_log"
    exit 1
fi

# 5. 运行第三步：R 语言数据清洗与终极归类
echo "[3/3] 正在使用 R 进行注释结果的终极归类与清洗..."

# 生成临时 R 脚本
TMP_R_SCRIPT="${out_dir}/.temp_clean_${base}.R"

cat << 'EOF' > "$TMP_R_SCRIPT"
# 接收 Bash 传来的参数
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# 抑制包加载信息的输出，保持终端整洁
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

# 1. 读取原始注释结果
anno <- read.delim(input_file, sep = "\t", header = TRUE, check.names = FALSE, comment.char = "")

# 2. 清洗逻辑
anno_final <- anno %>%
  mutate(
    clean_anno = `Detailed Annotation` %>%
      str_replace_all("^[^|]+\\|\\s*([^|]+)\\s*\\|.*$", "\\1") %>%
      str_replace_all("\\s*\\(.*\\)", "") %>%
      str_remove("\\..*$") %>%
      str_remove("\\?$") %>%
      str_trim(),

    final_category = case_when(
      clean_anno %in% c("promoter-TSS", "Promoter") ~ "Promoter",
      clean_anno %in% c("DNA", "DNA Transposon") ~ "DNA Transposon",
      clean_anno %in% c("intron", "Intron") ~ "Intron",
      clean_anno %in% c("intergenic", "Intergenic") ~ "Intergenic",
      clean_anno %in% c("exon", "Exon") ~ "Exon",
      clean_anno == "TTS" ~ "TTS",
      clean_anno %in% c("3' UTR", "3UTR", "3'UTR") ~ "3'UTR",
      clean_anno %in% c("5' UTR", "5UTR", "5'UTR") ~ "5'UTR",
      clean_anno %in% c("CpG Island", "CpG") ~ "CpG",
      clean_anno %in% c("LINE", "SINE", "LTR") ~ clean_anno, # 简写合并相同项
      clean_anno %in% c("non-coding", "ncRNA", "miRNA") ~ "Non-coding",
      clean_anno %in% c("RNA", "snRNA", "scRNA", "rRNA", "srpRNA") ~ "RNA repeats",
      clean_anno %in% c("Simple_repeat", "Low_complexity", "Satellite", "RC") ~ "Other repeats",
      clean_anno %in% c("Other", "Unknown", "pseudo", "tRNA") ~ "Other",
      TRUE ~ "Other"
    )
  )

# 3. 导出最终清洗文件
write.table(anno_final, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
EOF

# 运行临时 R 脚本并传递文件路径参数
Rscript "$TMP_R_SCRIPT" "$out_csv" "$out_cleaned"

# 检查 R 脚本是否成功执行
if [ $? -eq 0 ]; then
    echo "✅ 流程全部完成！"
    echo "📄 原始注释结果: $out_csv"
    echo "🌟 清洗归类结果: $out_cleaned"
    echo "📝 运行日志记录: $out_log"
    # 删除临时 R 脚本
    rm "$TMP_R_SCRIPT"
else
    echo "❌ R 语言清洗步骤失败，请检查 R 环境是否安装了 dplyr 和 stringr 包。"
    rm "$TMP_R_SCRIPT"
    exit 1
fi