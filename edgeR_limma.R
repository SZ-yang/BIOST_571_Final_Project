library(edgeR)

load("/Users/amywatt1/Documents/UW/BIOST_571_Final_Project/simulated_rnaseq_data.RData")
rnaseq_data_family$counts[1:5, 1:5]
head(rnaseq_data_family$design)



Y <- as.data.frame(rnaseq_data_family$counts)
X <- as.data.frame(rnaseq_data_family$design)

group <- factor(X$condition)
design <- model.matrix(~group)

#####################
## edgeR
#####################

# Create DGEList object
d0 <- DGEList(counts = t(Y)) # func takes a df with rows as genes and columns as samples

# Estimate dispersion parameter
d1 <- estimateDisp(d0, design, robust=TRUE) # estimate common dispersion and tagwise dispersion in one run

# Differential Expression
fit <- glmFit(d1)
lrt <- glmLRT(fit)
edgeR_res <- (lrt$table[1])[,1]

head(lrt$table)
head(edgeR_res)


#####################
## limma voom 
#####################

# Create DGEList object
d0 <- DGEList(counts = t(Y))

# Voom transformation and calculation of variance weights
v <- voom(d0, design)

# Fitting linear models in limma
fit <- lmFit(v, design)
fit <- eBayes(fit) # Empirical Bayes smoothing of standard errors
limma_res <- unname(fit$coefficients[,2])

head(fit$coefficients)
head(fit$p.value)
head(limma_res)
