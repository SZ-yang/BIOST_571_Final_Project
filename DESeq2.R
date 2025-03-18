library(DESeq2)
library(edgeR)
library(tidyr)

rnaseq_data_subject$counts[1:5, 1:5]
head(rnaseq_data_subject$design)

X <- as.data.frame(rnaseq_data_subject$design)
rownames(X) <- rep(1:nrow(X))

Y <- as.data.frame(rnaseq_data_subject$counts)
rownames(Y) <- rownames(X)

head(X)
head(t(Y)[1:5, 1:5])

condition <- factor(X$condition)
time <- factor(X$time)
design <- model.matrix(~ time*condition)

#####################
## DESeq
#####################

# Create DGEList object
d0 <- DGEList(counts = t(Y)) # func takes a df with rows as genes and columns as samples

# Analysis is performed on simulated data transformed using the vst func in DESeq2
Y_adjusted <- d0$counts + 1  # Add a pseudo-count to avoid zero values
dds <- DESeqDataSetFromMatrix(countData = Y_adjusted,
                              colData = X,
                              design = ~ time*condition)

# run DESeq
dds <- DESeq(dds)

# get results
de_res <- results(dds, tidy=TRUE)
rownames(de_res) <- de_res$row
head(de_res)

# padj and p-value
results_padj <- de_res %>%
  select(row, padj) %>%
  rename(Gene = row, deseq_padj = padj) %>%
  as.data.frame()

results_p <- de_res %>%
  select(row, pvalue) %>%
  rename(Gene = row, deseq_p = pvalue) %>%
  as.data.frame()

# Power
genes_true_de <- paste0("Gene", 1:500)
power_df <- results_p %>% filter(Gene %in% genes_true_de)

power_value <- sum(power_df$deseq_p < 0.05, na.rm = TRUE) / nrow(power_df)

# FDR
genes_non_de <- paste0("Gene", 501:nrow(results_p))
fdr_df <- results_p %>% filter(Gene %in% genes_non_de)

fdr_value <- sum(fdr_df$deseq_p < 0.05, na.rm = TRUE) / nrow(fdr_df)

results_df <- data.frame(
  Method = "DESeq2",
  Power = power_value,
  FDR = fdr_value
)

print(results_df)