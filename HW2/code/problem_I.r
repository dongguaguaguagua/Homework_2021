# 加载必要的库
library(org.Hs.eg.db)
library(tximport)
library(rtracklayer)
library(ggplot2)
library(sva)

### 1. 读取样本信息

# Directory where your Salmon output files are stored
salmon_dir <- "/home/huzongyao1/HW/HW2/raw_data1/HW2_sf"
# 获取所有样本的文件路径
sf_files <- list.files(salmon_dir, pattern = "*.sf", full.names = TRUE)
names(sf_files) <- c("N1", "N2", "N3", "N4", "N5", "T1", "T2", "T3", "T4", "T5")

# Load the GTF file and create a transcript-to-gene mapping
gtf <- import("/home/huzongyao1/HW/HW1/genome/Homo_sapiens.GRCh38.99.gtf")
# Filter for transcript and gene IDs
tx2gene <- as.data.frame(
    gtf[
        gtf$type == "transcript",
        c("transcript_id", "gene_id")
    ]
)
# Remove any duplicates
tx2gene <- unique(tx2gene[, c("transcript_id", "gene_id")])

### 2. 转换转录本水平为基因水平丰度

# 转换 gene_id 为基因符号 (SYMBOL)
# 使用 org.Hs.eg.db 进行转换
gene_symbol <- select(
    org.Hs.eg.db,
    keys = unique(tx2gene$gene_id),
    columns = "SYMBOL",
    keytype = "ENSEMBL"
)
# 将 gene_id 与基因符号合并
tx2gene <- merge(
    tx2gene,
    gene_symbol,
    by.x = "gene_id",
    by.y = "ENSEMBL",
    all.x = TRUE
)
tx2gene <- tx2gene[, c("transcript_id", "gene_id", "SYMBOL")]

### 3. 导入Salmon输出

txi <- tximport(
    sf_files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE,
    countsFromAbundance = "scaledTPM"
)

# Access gene-level TPM and log-transform
tpm_matrix <- txi$abundance
log2_tpm_matrix <- log2(tpm_matrix + 1)
gene_variances <- apply(log2_tpm_matrix, 1, var)
log2_tpm_filtered <- log2_tpm_matrix[gene_variances != 0, ]

### 4. 主成分分析：处理批次效应前

# 读取批次信息并进行PCA
batch_info <- read.csv("/home/huzongyao1/HW/HW2/raw_data1/batch.csv")
pca_result <- prcomp(t(log2_tpm_filtered), scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- rownames(pca_data)
pca_data <- merge(pca_data, batch_info, by.x = "Sample", by.y = "X")

### 5. 绘制PCA图

ggplot(pca_data, aes(x = PC1, y = PC2, shape = tissue, color = as.factor(batch))) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "PCA of Gene Expression Data", x = "PC1", y = "PC2") +
    theme_minimal() +
    scale_shape_manual(values = c(16, 17)) +
    scale_color_manual(values = c("blue", "red")) +
    theme(legend.title = element_blank())

### 6. 层次聚类

# 计算样本之间的距离矩阵
dist_matrix <- dist(t(log2_tpm_matrix))
# 进行层次聚类
hc <- hclust(dist_matrix)
# 绘制层次聚类树
plot(hc, labels = batch_info$X)

### 1. 选择那些通过方差过滤的基因

# 设定一个方差阈值（可以调整），只保留高变异的基因
variance_threshold <- quantile(gene_variances, 0.75) # 选择方差在前25%的基因
filtered_genes <- rownames(log2_tpm_matrix)[gene_variances > variance_threshold]
# 提取这些基因的TPM值矩阵
filtered_tpm_matrix <- log2_tpm_matrix[filtered_genes, ]

### 2. 使用COMBAT去除批次效应

# 从批次信息中提取批次
batch <- as.factor(batch_info$batch)
# 运行Combat算法，去除批次效应
combat_corrected <- ComBat(dat = filtered_tpm_matrix, batch = batch, par.prior = TRUE, prior.plots = FALSE)

### 3. 批次效应去除后的PCA

# 对经过Combat调整后的数据进行PCA
combat_pca_result <- prcomp(t(combat_corrected), scale. = TRUE)

# 提取PCA结果
combat_pca_data <- as.data.frame(combat_pca_result$x)
combat_pca_data$Sample <- rownames(combat_pca_data)
combat_pca_data$Tissue <- ifelse(grepl("^T", combat_pca_data$Sample), "Tumor", "Normal")
combat_pca_data$Batch <- batch_info$batch[match(combat_pca_data$Sample, batch_info$X)]

### 4. 可视化PCA结果

ggplot(combat_pca_data, aes(x = PC1, y = PC2, shape = Tissue, color = as.factor(Batch))) +
    geom_point(size = 3) +
    labs(title = "PCA After Batch Effect Removal", x = "PC1", y = "PC2") +
    theme_minimal()

### 5. 层次聚类分析

# 计算样本之间的距离矩阵
combat_dist_matrix <- dist(t(combat_corrected))

# 进行层次聚类
combat_hc <- hclust(combat_dist_matrix)

# 绘制层次聚类树
plot(combat_hc, labels = batch_info$Sample, main = "Hierarchical Clustering After Batch Effect Removal")
