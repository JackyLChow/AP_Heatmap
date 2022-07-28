library(readr)
library(DESeq2)

### load raw counts ---
raw_counts <- read_csv("Documents/BFX_proj/_data/Chow_Genetics_2018_GSE115657/raw_counts.csv")
c_mtx <- as.matrix(raw_counts[, 3:ncol(raw_counts)])
rownames(c_mtx) <- raw_counts$Symbol

### load metadata ---
meta <- data.frame(read_csv("Documents/BFX_proj/_data/Chow_Genetics_2018_GSE115657/GSE115657_meta.csv"))
rownames(meta) <- meta$GSM_id
meta$strain <- sub("_.*", "", meta$sample_description)
write.csv(meta, "~/Documents/BFX_proj/AP_Heatmap/_input/meta.csv", row.names = F)

### normalize counts in DESeq2 ---
ds2_ <- DESeqDataSetFromMatrix(countData = c_mtx, colData = meta, design = ~ 1) # make DESeq2 object
ds2_ <- DESeq(ds2_) # run DESeq2

### extract normalized counts ---
n_count <- counts(ds2_, normalized = T)
n_count <- data.frame(gene = rownames(n_count), n_count)

write.csv(n_count, "~/Documents/BFX_proj/AP_Heatmap/_input/Chow_Genetics_normcounts.csv", row.names = F)
