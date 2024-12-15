#!/bin/bash
# This script gathers all the output vcf file to a directory
# Author: Jason
# Date: 25/11/2024

# 定义结果目录和输出目录
RESULTS_DIR="/hdd1/home/f24_wyli/project/pipeline/analysis_results/results" 
OUTPUT_DIR="/hdd1/home/f24_wyli/project/pipeline/vcf_gathering"   
# 如果输出目录不存在，则创建
mkdir -p "${OUTPUT_DIR}"

# 遍历结果目录下的所有 SRR 子目录
for SRR_DIR in "${RESULTS_DIR}"/SRR*; do
    if [ -d "${SRR_DIR}" ]; then
        # 提取 SRR 编号
        SRR=$(basename "${SRR_DIR}")
        
        # 查找该目录下的 VCF 文件（假设扩展名为 .vcf）
        VCF_FILE=$(find "${SRR_DIR}/output/output.vcf" -maxdepth 1 -type f -name "*.vcf")
        
        if [ -f "${VCF_FILE}" ]; then
            # 将 VCF 文件复制到输出目录，重命名为 SRR 编号命名的文件
            cp "${VCF_FILE}" "${OUTPUT_DIR}/${SRR}.vcf"
            echo "已复制 ${VCF_FILE} 到 ${OUTPUT_DIR}/${SRR}.vcf"
        else
            echo "在 ${SRR_DIR} 中未找到 VCF 文件"
        fi
    fi
done

