library(edgeR)
library(DESeq2)
library(lmerSeq)
library(tidyr)

load("/Users/amywatt1/Documents/UW/BIOST_571_Final_Project/simulated_rnaseq_data.RData")
rnaseq_data_family$counts[1:5, 1:5]
head(rnaseq_data_family$design)


X <- as.data.frame(rnaseq_data_subject$design)
rownames(X) <- rep(1:nrow(X))

Y <- as.data.frame(rnaseq_data_subject$counts)
rownames(Y) <- rownames(X)

head(X)
head(t(Y)[1:5, 1:5])

group <- factor(X$condition)
design <- model.matrix(~group)

#####################
## edgeR
#####################

# Create DGEList object
d0 <- DGEList(counts = t(Y)) # func takes a df with rows as genes and columns as samples

# Filter genes
keep <- filterByExpr(d0, design)
d0 <- d0[keep, , keep.lib.sizes=FALSE]

# Normalization
d0 <- normLibSizes(d0)

# Estimate dispersion parameter
d1 <- estimateDisp(d0, design, robust=TRUE) # estimate common dispersion and tagwise dispersion in one run

# Differential Expression
fit <- glmFit(d1)
lrt <- glmLRT(fit)
edgeR_res <- topTags(lrt, adjust.method = "BH", sort.by = "none", n = Inf)

head(edgeR_res)

#####################
## limma voom 
#####################

# Create DGEList object
d0 <- DGEList(counts = t(Y))

# Filter genes
keep <- filterByExpr(d0, design)
d0 <- d0[keep, , keep.lib.sizes=FALSE]

# Normalization
d0 <- normLibSizes(d0)

# Voom transformation and calculation of variance weights
v <- voom(d0, design)

# Fitting linear models in limma
fit <- lmFit(v, design)
fit <- eBayes(fit) # Empirical Bayes smoothing of standard errors
limma_res <- topTable(fit, adjust.method = "BH", sort.by = "none", n = Inf)

head(limma_res)

#####################
## DESeq
#####################

# Create DGEList object
d0 <- DGEList(counts = t(Y)) # func takes a df with rows as genes and columns as samples

# Filter genes
keep <- filterByExpr(d0, design)
d0 <- d0[keep, , keep.lib.sizes=FALSE]

# Analysis is performed on simulated data transformed using the vst func in DESeq2
Y_adjusted <- d0$counts + 1  # Add a pseudo-count to avoid zero values
dds <- DESeqDataSetFromMatrix(countData = Y_adjusted,
                              colData = X,
                              design = ~ condition)

# run DESeq
dds <- DESeq(dds)

# get results
de_res <- results(dds, tidy=TRUE)
head(de_res)

#####################
## lmerSeq 
#####################

# Create DGEList object
d0 <- DGEList(counts = t(Y)) # func takes a df with rows as genes and columns as samples

# Filter genes
keep <- filterByExpr(d0, design)
d0 <- d0[keep, , keep.lib.sizes=FALSE]

# Analysis is performed on simulated data transformed using the vst func in DESeq2
Y_adjusted <- d0$counts + 1  # Add a pseudo-count to avoid zero values
dds <- DESeqDataSetFromMatrix(countData = Y_adjusted,
                              colData = X,
                              design = ~ condition)
dds <- DESeq(dds)
vsd.fixed <- varianceStabilizingTransformation(dds, blind=F)
vst_expr <- assay(vsd.fixed)


# Fit the Model
fit.lmerSeq <- lmerSeq.fit(form = ~ condition + (1|subject_id),
                           expr_mat = vst_expr,
                           sample_data = X,
                           REML = T)

##  Summarize the condition coefficients
lmer_res <- lmerSeq.summary(lmerSeq_results = fit.lmerSeq,
                                   coefficient = "conditioncontrol",
                                   p_adj_method = 'BH',
                                   ddf = 'Satterthwaite',
                                   sort_results = F, 
                                   include_singular = TRUE)

head(lmer_res$summary_table)


#########################
## Put the results together
#########################

results <- data.frame(
  edgeR = edgeR_res$table$FDR, 
  limma = limma_res$adj.P.Val,
  deseq = de_res$padj,
  lmerSeq = lmer_res$summary_table$p_val_adj
)

rownames(results) <- rownames(edgeR_res$table)


### calculate power for each method
# Step 1: Subset the rows gene1 to gene500
genes_to_select <- paste0("Gene", 1:500)
power_df <- results[genes_to_select, ] 

# Step 2: Apply the condition per column and count values <= 0.05
power_count_values_per_column <- apply(power_df, 2, function(x) sum(x <= 0.05))

# Output the result
power_count_values_per_column / nrow(power_df)

### calculate FDR for each method
# Step 1: Subset the rows gene1 to gene500
genes_to_select <- paste0("Gene", 501:nrow(results))
fdr_df <- results[genes_to_select, ] 

# Step 2: Apply the condition per column and count values <= 0.05
fdr_count_values_per_column <- apply(fdr_df, 2, function(x) sum(x <= 0.05, na.rm = TRUE))

# Output the result
fdr_count_values_per_column / nrow(fdr_df)

data.frame(
  power = power_count_values_per_column / nrow(power_df),
  fdr = fdr_count_values_per_column / nrow(fdr_df)
)
