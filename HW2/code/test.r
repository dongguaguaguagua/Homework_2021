# 加载必要的库
library(org.Hs.eg.db)
library(tximport)
library(rtracklayer)

# Directory where your Salmon output files are stored
salmon_dir <- "/home/huzongyao1/HW/HW2/raw_data1/HW2_sf"
# List all the Salmon quant.sf files
files <- list.files(salmon_dir, pattern = "N1.sf", full.names = TRUE)
# Create a transcript-to-gene mapping
# Load the GTF file
gtf_file <- "/home/huzongyao1/HW/HW1/genome/Homo_sapiens.GRCh38.99.gtf"
gtf <- import(gtf_file)
# Filter for transcript and gene IDs
tx2gene <- as.data.frame(gtf[gtf$type == "transcript", c("transcript_id", "gene_id")])
# Remove any duplicates
tx2gene <- unique(tx2gene[, c("transcript_id", "gene_id")])
# Import the Salmon output
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
# Access gene-level TPM
gene_tpm <- txi$abundance # TPM values
print(gene_tpm)
# 示例 Ensembl ID
gene_ids <- c("ENSG00000139618", "ENSG00000157764")

# 从 Ensembl ID 映射到基因符号
gene_symbols <- mapIds(org.Hs.eg.db,
    keys = gene_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
)

# 查看结果
print(gene_symbols)
