import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# 文件路径列表，包含所有需要解析的geneBodyCoverage.r文件
folder = "/home/huzongyao1/Homework_2020/HW1/output/bam"
file_paths = [
    f"{folder}/res_WAligned.sortedByCoord.txt.geneBodyCoverage.r",
    f"{folder}/res_XAligned.sortedByCoord.txt.geneBodyCoverage.r",
    f"{folder}/res_YAligned.sortedByCoord.txt.geneBodyCoverage.r",
    f"{folder}/res_ZAligned.sortedByCoord.txt.geneBodyCoverage.r",
]

# 用于存储样本名和覆盖度数据的字典
data_dict = {}

# 从每个文件中提取样本名称和标准化基因覆盖度数值
for file_path in file_paths:
    with open(file_path, 'r') as file:
        lines = file.readlines()
        sample_name, data_str = lines[0].split(' <- c(')
        data_str = data_str.strip().strip(')')
        values = [float(x) for x in data_str.split(',')]
        data_dict[sample_name] = values

# 转换数据为DataFrame格式
df = pd.DataFrame(data_dict)

# 绘制覆盖度曲线图（2x2布局）
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle("Gene Body Coverage Across Samples", fontsize=16)

# 将每个样本绘制在子图中
# print(data_dict.items())
for i, (sample_name, values) in enumerate(data_dict.items()):
    ax = axes[i // 2, i % 2]  # 2行2列布局
    ax.plot(np.linspace(0, 100, len(values)), values, label=sample_name)
    ax.set_title(f"{sample_name}")
    ax.set_xlabel("Gene Body percentile (5'->3')")
    ax.set_ylabel("Normalized Read Count")
    ax.legend()

# 检查覆盖度最低的样本
min_coverage_sample = df.mean().idxmin()
print(f"Sample with the lowest average gene body coverage: {min_coverage_sample}")

plt.tight_layout(rect=[0, 0, 1, 0.95])  # 调整布局以适应标题
# 存入文件
plt.savefig("gene_body_coverage.png")
