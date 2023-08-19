## 26/07/2023
## Harvard Course
## https://hbctraining.github.io/DGE_workshop_salmon/lessons/01_DGE_setup_and_overview.html
## Gene-level differential expression analysis using DESeq2
## Loading libraries (Bioconductor, CRAN)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(AnnotationHub)
library(ensembldb)
library(GenomicFeatures) #M. Ono Omics
library(fgsea) # gene set enrichment analysis # M. Ono Omics
library(ReactomePA) # reactome pathway analysis # M. Ono Omics
library(VennDiagram) # venn diagram analysis # M. Ono Omics
library(Rsubread) # Jenna
library(gProfileR)
library(treemap)
library(SPIA)
library(apeglm)
library(annotables) # need to install
library(dplyr)
library(Biobase)
library(GenomeInfoDb)
library(IRanges)
library(S4Vectors)
library(GenomicRanges)
library(AnnotationHub)
library(ensembldb)

## Jenna - command line
# make feature count table after using STAR
# nice R 
# fc_PEO_RNA <- featureCounts(files=list.files(pattern='.bam$'),
                            # annot.ext='/data/jr2715/PEO1_4/RNA_seq/star/Homo_sapiens.GRCh38.104.gtf',
                            # isGTFAnnotationFile=TRUE, nthread=32, GTF.featureType="exon", 
                            # GTF.attrType="gene_id", useMetaFeatures=TRUE, isPairedEnd=T)
# define comparison 
# condition <- c("PEO1", "PEO1", "PEO1", "PEO4", "PEO4", "PEO4")
# run DESeq2 
# dds <- DESeq(dds)
# extract results 
# res <- results(dds)

# Harvard
## list all directories containing data
samples <- list.files(path = "./data", full.names = T, pattern="salmon$")

# obtain vector of all filenames including path
files <- file.path(samples, "quant.sf")

# all quant files have same name, useful to have names for each element
names(files) <- str_replace(samples, "./data/", "") %>% 
  str_replace(".salmon", "")

# M. ONO Omics
#files <- filepath ("." , c("DP_B_Het_R1", "DP_N_Het_R1"), "quant.sf")
#names (files) <- c("DP_B_Het_R1", "DP_N_Het_R1")

# load gene annotation for GrCh38
tx2gene <- read.delim("C:/Users/Acer/Desktop/DEanalysis_Harvard/data/
                      tx2gene_grch38_ens94.txt")

# M. ONO Omics
# txi. salmon <- tximport (files, type = "salmon", tx2gene = tx2gene, 
# ignoreTxVersion = TRUE)
# gets around Ensembl version numbers 

# view it
tx2gene %>% View()

# M. ONO Omics
# head(txi.salmon$counts)

# ?tximport parameters
# running it
txi <- tximport(files, type="salmon", 
                tx2gene=tx2gene[,c("tx_id", "ensgene")], 
                countsFromAbundance="lengthScaledTPM")

# view data
attributes(txi)

# look at counts
txi$counts %>% View()

# write counts to an object
data <- txi$counts %>% round() %>% data.frame()

# creating metadata
sampletype <- factor(c(rep("control",3), rep("MOV10_knockdown", 2), 
                       rep("MOV10_overexpression", 3)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

# M. ONO Omics - loading preprocessed count data
# counts <- read.table(file = "path/counts.csv", header = TRUE, sep = ",")
# metadata <- read.table (file = "path/metadata.csv", header = TRUE, sep = ",")

# inspecting by dim and summary
dim(meta)
dim(data)
summary(meta)
summary(data)

# distribution of counts barchart
ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

## the type of distribution of RNA-seq data 
## either Poisson or Negative Binomial depending on mean/variance relationship
#plotting mean vs variance
mean_counts <- apply(data[,6:8], 1, mean)        
# The second argument '1' of 'apply' function indicates the 
# function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

# variance across replicates tends to be greater than the mean (red line), 
# especially for genes with large mean expression levels.

# lowly expressed genes = scatter, “heteroscedasticity”,
# for a given expression level = high variation in variance.

# If mRNA proportions stayed exactly constant biological 
# replicates for a sample group = Poisson distribution (mean == variance)
# more replicates = more like Poisson

# this data = not fit Poisson distribution
# here = Negative Binomial -> mean < variance

# M. ONO Omics - did filter out lowest detectable reads
# logic <- rowSums ( counts ) > 10
# counts <- counts [logic ,]

# distinguish biological vs technical replicates (biological are more useful)
# cell lines = https://web.archive.org/web/20170807192514/http://www.labstats.net:80/articles/cell_culture_n.html

# increase in the number of replicates tends to 
# return more DE genes than increasing the sequencing depth

# DESeq2, EdgeR, NBPSeq, baySeq, EBSeq, vst, TSPM 
# Limma-Voom (less sensitive, for 20 per group replicates, big datasets)

## COUNT NORMALIZATION
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
# why - sequencing depth, gene length, RNA composition
# RPKM/FPKM (not recommended for between sample comparisons)

# DESeq - median of ratios
# 1) creates a pseudo-reference sample (row-wise geometric mean)
# 2) calculates ratio of each sample to the reference
# 3) calculates the normalization factor for each sample (size factor)
#  normalization_factor_sampleA <- median(c(1.28, 1.3, 1.39, 1.35, 0.59))
# 4) calculate the normalized count values using the normalization factor
# median of ratios = not ALL genes are differentially expressed; 
# the normalization factors account for sequencing depth and RNA composition t
# robust to imbalance in up-/down-regulation and large numbers of differentially expressed genes.

# match metadata with counts
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))

# creating a DESeq object (data set from matrix)
# dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

# view the thing
View(counts(dds))
dds
dim(dds)
class(dds)
colData(dds)

# median of ratios method - can be done in single stap by DESeq()
# in different steps:
# estimate the size factor
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
# retriving counts back
normalized_counts <- counts(dds, normalized=TRUE)
# save into a file later
write.table(normalized_counts, file="data/normalized_counts.txt", 
            sep="\t", quote=F, col.names=NA)
# it uses the raw counts and models the normalization 
# inside the Generalized Linear Model (GLM).

dds <- DESeq(dds)

# Genes with zero counts , with an extreme count outlier, with a low mean normalized counts
# DESeq filters them out automatically

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE) # slow with more than 20 samples
# with large datasets use vst()
# blind=TRUE argument results in a transformation unbiased to sample condition information.

# M. ONO Omics
# Obtain a gene expression matrix
# vsd <- vst(dds , blind = FALSE )

## PCA - by which factor to cluster the data
plotPCA(rld, intgroup="sampletype")
# By default the function uses the top 500 most variable genes. 
# You can change this by adding the ntop

# additional PCAs
# Input is a matrix of log transformed values
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))

# M. ONO Omics
# pca_plot <- plotPCA (vsd, intgroup = c("sample"))
# ggsave(filename = "output/pca_plot.pdf", plot = pca_plot, device = cairo_pdf, dpi = 300)
## prcomp
# Obtain a gene expression matrix
# exprs.data <- assay(vsd)
# Inspect exprs.data
# exprs.data [1:3, 1:3]
# pca_prcomp <- prcomp (t(exprs.data), scale = FALSE)
# Plot PCA result
# pdf ("output/PCA.pdf")
# par(mfrow =c(2,2))
# plot(pca_prcomp$x[ ,1:2], main =’Sample plot’, col = colData (vsd)$sample, pch = 19)
# text( pca_prcomp$x[ ,1:2], labels = rownames (pca_prcomp$x[ ,1:2]), pos = 4, cex = 0.8)
# abline (v=0, h=0, col=8)
# plot (pca_prcomp$rotation[ ,1:2], main = "Gene loading")
# abline (v=0, h=0, col=8)
# eig <- (pca_prcomp$sdev)^2; names(eig) <- 1:length(eig); eig <- 100*eig/sum(eig)
# barplot (eig, main = "Scree plot", xlab = "PC", ylab = "Percentage of explained variances")
# dev.off ()

## heatmap - the correlation of gene expression for all pairwise combinations of samples 
# Extract the rlog matrix from the object
rld_mat <- assay(rld)  # assay() is function from the "SummarizedExperiment" 
# package that was loaded when you loaded DESeq2
# Compute pairwise correlation values
rld_cor <- cor(rld_mat) # cor() is a base R function
head(rld_cor) # check the output of cor(), make note of the rownames and colnames
# Plot heatmap
pheatmap(rld_cor, annotation = meta)
# ?pheatmap
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)

## HYPOTHESIS TESTING
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Wald test is the default used for hypothesis testing when comparing two groups.
# Taking LFC, dividing it by its SE, resulting in a z-statistic
# compare z-statistic to standard normal distribution, compute p-value 
# (probability that a z-statistic at least as extreme as the observed value would be selected at random)
# If p-value is small = reject null hypothesis (i.e. the gene is differentially expressed).

# likelihood ratio test - compare more than 2 sample classes, genes changing in expression in any direction
# Wald - only one model per gene, LRT 0 two models per gene, fits are compared
# null hypothesis that the full model fits just as well as the reduced model
# Likelihood ratio test
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)

# each gene that has been tested will be associated with a p-value. It is this result which we use to 
# determine which genes are considered significantly differentially expressed. 
# However, we cannot use the p-value directly. Since each p-value is the result of a single test (single gene). 
# The more genes we test, the more we inflate the false positive rate. This is the multiple testing problem.

# correct p value for multiple testing
# Bonferroni - This is a very conservative approach with a high probability of false negatives
# FDR/Benjamini-Hochberg - control the expected FDR below a specified level given a list of independent p-values
# # Q-value / Storey method: The minimum FDR that can be attained when calling that feature significant. For example, 
# if gene X has a q-value of 0.013 it means that 1.3% of genes that show p-values at least as 
# small as gene X are false positives.

# By setting the FDR cutoff to < 0.05, we’re saying that the proportion of false positives we expect 
# amongst our differentially expressed genes is 5%. For example, if you call 500 genes as differentially 
# expressed with an FDR cutoff of 0.05, you expect 25 of them to be false positives.

# Check the size factors
sizeFactors(dds)

# Total number of raw counts per sample
colSums(counts(dds))

# Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

# results table
# Define contrasts, extract results table, and shrink the log2 fold changes
contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
res_tableOE <- results(dds, contrast=contrast_oe, alpha = 0.05)
# explore results
class(res_tableOE)
mcols(res_tableOE, use.names=T)
res_tableOE %>% data.frame() %>% View()
 
# The order of the names in the contrast determines the direction of fold change that is reported. 
# The name provided in the second element is the level that is used as baseline. So for example, 
# if we observe a log2 fold change of -2 this would mean the gene expression is lower in Mov10_oe 
# relative to the control. However, these estimates do not account for the large dispersion we observe
# with low read counts. To avoid this, the log2 fold changes calculated by the model need to be adjusted.

# generate more accurate log2 foldchange estimates, DESeq2 allows for the shrinkage of the LFC estimates toward zero
# Save the unshrunken results to compare
res_tableOE_unshrunken <- res_tableOE
# Apply fold change shrinkage
res_tableOE <- lfcShrink(dds, coef = "sampletype_MOV10_overexpression_vs_control", type = "apeglm", res=res_tableOE)
# browseVignettes("apeglm"), https://www.biostars.org/p/448959/
# Shrinking the log2 fold changes will not change the total number of genes that are identified 
# as significantly differentially expressed. The shrinkage of fold change is to help with downstream 
# assessment of results.
# use MA plot: mean of the normalized counts versus the log2 foldchanges for all genes tested
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
plotMA(res_tableOE, ylim=c(-2,2))

# LRT vs Wald test
# Likelihood ratio test
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
# Extract results
res_LRT <- results(dds_lrt)
# Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
# Subset to return genes with padj < 0.05
sigLRT_genes <- res_LRT_tb %>% 
  dplyr::filter(padj < padj.cutoff)
# Get number of significant genes
nrow(sigLRT_genes) # 7315
# Compare to numbers we had from Wald test
nrow(sigOE) # 4808
nrow(sigKD) # 2810

# DEA control vs knockdown
contrast_kd <-  c("sampletype", "MOV10_knockdown", "control")
res_tableKD <- results(dds, contrast=contrast_kd, alpha = 0.05)
res_tableKD <- lfcShrink(dds, coef = "sampletype_MOV10_knockdown_vs_control", type = "apeglm", res=res_tableKD)
# Summarize results
summary(res_tableOE, alpha = 0.05)
# Set thresholds
padj.cutoff <- 0.05
# filter
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
# subset only significant genes as per thresholds
sigOE <- res_tableOE_tb %>%
  dplyr::filter(padj < padj.cutoff)
sigOE
# use same p-adjusted adj.cutoff < 0.05), subset res_tableKD to report 
# the number of genes that are up- and down-regulated in Mov10_knockdown compared to control.
res_tableKD_tb <- res_tableKD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sigKD <- res_tableKD_tb %>%
  dplyr::filter(padj < padj.cutoff)
sigKD
# With large significant gene lists it can be hard to extract meaningful biological relevance. 
# To help increase stringency, one can also add a fold change threshold.
# ex. a new threshold lfc.cutoff and set it to 0.58 (remember that we are working 
# with log2 fold changes so this translates to an actual fold change of 1.5).


## M. ONO Omics - pairwise comparison
# Choose the reference sample group - use the exact name of a sample group
reference <- "MOV10_knockdown"
# Choose the target group
target <- "control"
res_ono <- results(dds, contrast = c("sampletype", target, reference))
# extract DEGS
result_ono <- as.data.frame(res_ono) # convert to a data frame object
# Remove NA
lg1 <- !is.na(result_ono$log2FoldChange)
lg2 <- !is.na(result_ono$padj)
result_ono <- result_ono[lg1&lg2,]
# Filter results by padj (adjusted p-value) and log2 fold change (absolute level)
lg3 <- result_ono$padj < 0.05
lg4 <- abs(result_ono$log2FoldChange) > 0
result_ono <- result_ono[lg3&lg4,]
# Order the result according to padj
ord_ono <- order(result_ono$padj, decreasing = FALSE)
result_ono <- result_ono[ord_ono,]
# Inspect
result_ono[1:3,]
summary(result_ono)

## ADVANCED VISUALIZATIONS
# metadata tibble
mov10_meta <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

# gene names, column with gene symbols
normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  left_join(tx2gene, by=c("gene" = "ensgene"))

# plotting one gene
# Find the Ensembl ID of MOV10
tx2gene[tx2gene$symbol == "MOV10", "ensgene"]
# Plot expression for single gene
d <- plotCounts(dds, gene="ENSG00000155363", intgroup="sampletype", returnData=TRUE) 
# Plot the MOV10 normalized counts, using the samplenames (rownames(d) as labels)
ggplot(d, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("MOV10") +
  theme(plot.title = element_text(hjust = 0.5))

# CONSTRUCTING A HEATMAP
# Extract normalized expression for significant genes from the OE and control samples (2:4 and 7:9)
norm_OEsig <- normalized_counts[,c(1:4,7:9)] %>% 
  dplyr::filter(gene %in% sigOE$gene) 
# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")
# Run pheatmap using the metadata data frame for the annotation
gene_map <- pheatmap(norm_OEsig[2:7], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", #scale="row", in which Z-scores are plotted,
         # rather than the actual normalized count value.
         fontsize_row = 10, 
         height = 20)
gene_map

grch38

# M. ONO Omics
# Calculate the variance of each gene expression
# exprs.data=assay(vsd)[[1]]
# var_genes=apply(exprs.data, 1, var)
# col.group=brewer.pal(10, "Set3")[as.factor(colData(vsd)$group)]
# Get the names for the top 100 most highly variable genes and subset data
# select_var=names(sort(var_genes, decreasing=TRUE))[1:100]
# hi_var=exprs.data[select_var,]
# dim(hi_var) # Get gene number and sample number as output
# pdf("output4/heatmap2.pdf")
# heatmap.2(as.matrix(hi_var), col=colorRampPalette(brewer.pal(10, "RdYlBu"))(256), trace="none",
# srtCol = 45, margins = c(12,9), cexCol = 1.2,
# density.info = "none", lhei = c(2, 8), keysize = 0.9,
# main="Top 100 Most Variable Genes", scale="row", ColSideColors=col.group)
# dev.off()

## VOLCANO PLOT
# Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_tableOE_tb <- res_tableOE_tb %>% 
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)
# Volcano plot
ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
# labelling genes
# Add all the gene symbols as a column from the grch38 table using bind_cols()
res_tableOE_tb <- bind_cols(res_tableOE_tb, symbol=tx2gene$symbol[match(res_tableOE_tb$gene, tx2gene$ensgene)])
# Create an empty column to indicate which genes to label
res_tableOE_tb <- res_tableOE_tb %>% mutate(genelabels = "")
# Sort by padj values 
res_tableOE_tb <- res_tableOE_tb %>% arrange(padj)
# Populate the genelabels column with contents of the gene symbols column for the first 10 rows, 
#i.e. the top 10 most significantly expressed genes
res_tableOE_tb$genelabels[1:10] <- as.character(res_tableOE_tb$symbol[1:10])
View(res_tableOE_tb)
ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_OE)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

# M. ONO Omics
# test = ifelse(result$padj < 0.01, TRUE, FALSE)
# upregulated = ifelse(test = result$log2FoldChange > 1, TRUE, FALSE)
# downregulated = ifelse(test = result$log2FoldChange < -1, TRUE, FALSE)
#ifelse is used to define the colour - check the logical operators
# state = ifelse(test == TRUE & upregulated == TRUE, 'Upregulated', 
#ifelse(test == TRUE & downregulated == TRUE, 'Downregulated','Not statistically significant'))
# genes = rownames(result)[test & (upregulated|downregulated)]
#Define which gene names are shown
# point_label = rownames(result)
# point_label[!(rownames(result) %in% genes)] = ""
# volcano_plot <- ggplot(data = result, mapping = aes(x = log2FoldChange, y = -log10(padj))) +
  # geom_point(mapping = aes(color = state, alpha = 0.9), size = 2) +
  # geom_text_repel(mapping = aes(label = point_label), max.overlaps=20) +
  # guides(alpha = "none") +
  # labs(x = 'log2(fold change)',
       # y = '-log10(adjusted p-value)') +
  # theme_minimal() +
  # ggtitle(paste0(target,' minus ', reference))
# ggsave(filename = paste0('output4/',target,'_minus_', reference, "_volcanoplot_SP_B_N.pdf"), 
# plot = volcano_plot, device = cairo_pdf, dpi = 300)

# fever lines of code with a package #NOT WORKING
DEGreport::degPlot(dds = dds, res = results, n = 20, xs = "type", group = "condition") # dds object is output from DESeq2
DEGreport::degVolcano(
  data.frame(results[,c("log2FoldChange","padj")]), # table - 2 columns
  plot_text = data.frame(res[1:10,c("log2FoldChange","padj","id")])) # table to add names
# Available in the newer version for R 3.4
DEGreport::degPlotWide(dds = dds, genes = row.names(results)[1:5], group = "condition")

## GENE CLUSTER ANALYSIS
# Subset results for faster cluster finding (for classroom demo purposes)
clustering_sig_genes <- sigLRT_genes %>%
  arrange(padj) %>%
  head(n=1000)
# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
# sets of genes that exhibit similar expression patterns across sample groups. 
# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
# What type of data structure is the `clusters` output?
class(clusters)
# Let's see what is stored in the `df` component
head(clusters$df)
# Extract the Group 1 genes
cluster_groups <- clusters$df
group1 <- clusters$df %>%
  filter(cluster == 1)
# The LRT test can be especially helpful when performing time course analyses. 
# We can use the LRT to explore whether there are any significant differences in 
# treatment effect between any of the timepoints.

# GENOMIC ANNOTATION
# Load libraries
# Connect to AnnotationHub
ah <- AnnotationHub()
ah
library(AnnotationHub) # no errors created here
sessionInfo()
# Access previous build of annotations
grch38
# Explore all species information available
unique(ah$species) %>% View()
# Explore the types of Data Objects available
unique(ah$rdataclass) %>% View()
# Explore the Data Providers
unique(ah$dataprovider) %>% View()
# Query AnnotationHub
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens
# Extract annotations of interest
human_ens <- human_ens[["AH75011"]]
# Extract gene-level information
genes(human_ens, return.type = "data.frame") %>% View()
# Extract transcript-level information
transcripts(human_ens, return.type = "data.frame") %>% View()
# Extract exon-level information
exons(human_ens, return.type = "data.frame") %>% View()
# Create a gene-level dataframe 
annotations_ahb <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>% 
  dplyr::filter(gene_id %in% res_tableOE_tb$gene)
# Wait a second, we don't have one-to-one mappings!
class(annotations_ahb$entrezid)
which(map(annotations_ahb$entrezid, length) > 1)
annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()
# OrgDb and not Ensdb
human_orgdb <- query(ah, c("Homo sapiens", "OrgDb"))
human_orgdb <- human_ens[["AH75742"]] # not s ubsettable
annotations_orgdb <- select(human_orgdb, res_tableOE_tb$gene, c("SYMBOL", "GENENAME", "ENTREZID"), "ENSEMBL")

which(is.na(annotations_ahb$gene_name)) %>% length()
which(duplicated(annotations_ahb$gene_name)) %>% length()
# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)
# How many rows does annotations_ahb have?
annotations_ahb %>% nrow()
# Return only the non-duplicated genes using indices
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]
# How many rows are we left with after removing?
annotations_ahb %>% nrow()
# Determine how many of the Entrez column entries are NA
which(is.na(annotations_ahb$entrezid)) %>%  length()

# to create a tx2gene file:
# Create a transcript dataframe
# txdb <- transcripts(human_ens, return.type = "data.frame") %>%
 # dplyr::select(tx_id, gene_id)
# txdb <- txdb[grep("ENST", txdb$tx_id),]

# Create a gene-level dataframe
# genedb <- genes(human_ens, return.type = "data.frame")  %>%
  # dplyr::select(gene_id, gene_name)

# Merge the two dataframes together
# annotations <- inner_join(txdb, genedb)


# PATHWAY AND FUNCTIONAL ANALYSIS
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html
# reactome, GSEA, Panther, DAVID, Diffbind for ChIP and CutNTag
# determine whether there is enrichment of known biological functions, interactions, or pathways
# identify genes’ involvement in novel pathways or networks by grouping genes together based on similar trends
# use global changes in gene expression by visualizing all genes being significantly up-
# or down-regulated in the context of external interaction data
# you should NOT use these tools to make conclusions about the pathways involved 
# in your experimental process. You will need to perform experimental validation of any suggested pathways.
# overrepresentation analysis, hypergeometric testing
# gene ontology, cluster profiler
# Cluster Profiler
# Merge the AnnotationHub dataframe with the results 
res_ids <- inner_join(res_tableOE_tb, annotations_ahb, by=c("gene"="gene_id"))    
# Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allOE_genes <- as.character(res_ids$gene)

## Extract significant results
sigOE <- dplyr::filter(res_ids, padj < 0.05)

sigOE_genes <- as.character(sigOE$gene)
# Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfiler_Mov10oe.csv")
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
# visualizing results
# Dotplot 
# dotplot shows the number of genes associated with the first 50 terms (size)
# and the p-adjusted values for these terms (color). This plot displays the top
# 50 GO terms by gene ratio (# genes related to GO term / total number of sig genes)
# not p-adjusted value.
dotplot(ego, showCategory=50)
# Save as PDF.... In the pop-up window, change:
# Orientation: to Landscape
# PDF size to 8 x 14 to give a figure of appropriate size for the text labels
# enrichment GO plot - not working
# shows the relationship between the top 50 most significantly
# enriched GO terms (padj.), by grouping similar terms together. 
# Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
# compare_cluster <- compareCluster(geneClusters = ego, 
                                       # fun = "enrichKEGG",
                                       # organism = "hsa",
                                       # pAdjustMethod = "BH",
                                       # universe = ego$entrezgene_id,
                                       # qvalueCutoff = 0.05)
# x2 <- pairwise_termsim(compare_cluster) 
# emapplot(x2)

# category netplot -  the relationships between the genes associated with the 
# top five most significant GO terms and the fold 
# changes of the significant genes associated with these terms
# To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene
# Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

# If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
OE_foldchanges <- ifelse(OE_foldchanges > 2, 2, OE_foldchanges)
OE_foldchanges <- ifelse(OE_foldchanges < -2, -2, OE_foldchanges)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

# Subsetting the ego results without overwriting original `ego` variable
ego2 <- ego
ego2@result <- ego@result[c(1,3,4,8,9),]
# Plotting terms of interest
cnetplot(ego2, 
         categorySize="pvalue", 
         foldChange=OE_foldchanges, 
         showCategory = 5, 
         vertex.label.font=6)

# GSEA
# gene-level statistics or log2 fold changes for all genes from the differential
# expression results, then look to see whether gene sets for particular biological
# pathways are enriched among the large positive or negative fold changes.
# The hypothesis of FCS methods is that although large changes in individual
# genes can have significant effects on pathways (and will be detected via ORA methods),
# weaker but coordinated changes in sets of functionally related genes 
# (i.e., pathways) can also have significant effects. 
# This type of analysis can be particularly helpful if the differential expression
# analysis only outputs a small list of significant DE genes.
#  Commonly used gene sets include those derived from KEGG pathways, 
# Gene Ontology terms, MSigDB, Reactome, or gene groups that share some
# other functional annotations, etc. Consistent perturbations over such gene
# sets frequently suggest mechanistic changes” [1].
# KEGG
# Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")
# Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]
# Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange
# Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrezid
# Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)
head(foldchanges)
# GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)
# Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result
# Write GSEA results to file
View(gseaKEGG_results)
write.csv(gseaKEGG_results, "results/gseaOE_kegg.csv", quote=F)
# We will all get different results for the GSEA because the permutations 
# performed use random reordering. If we would like to use the same permutations
# every time we run a function (i.e. we would like the same results every time
# we run the function), then we could use the set.seed(123456)
# Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa03040')
# integrate KEGG pathway data 
detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts
# Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
         pathway.id = "hsa03040",
         species = "hsa",
         limit = list(gene = 2, # value gives the max/min limit for foldchanges
                      cpd = 1))
# Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
  pathview(gene.data = foldchanges, pathway.id = gseaKEGG_results$ID[x], species = "hsa", 
           limit = list(gene = 2, cpd = 1))
}

purrr::map(1:length(gseaKEGG_results$ID), get_kegg_plots)
# the enrichment of BP Gene Ontology terms using gene set enrichment analysis:
# GSEA using gene sets associated with BP Gene Ontology terms
gseaGO <- gseGO(geneList = foldchanges, 
                OrgDb = org.Hs.eg.db, 
                ont = 'BP', 
                nPerm = 1000, 
                minGSSize = 20, 
                pvalueCutoff = 0.05,
                verbose = FALSE) 
gseaGO_results <- gseaGO@result
gseaplot(gseaGO, geneSetID = 'GO:0007423')
#  it is possible to supply your own gene set GMT file, such as a GMT for MSigDB using special clusterProfiler functions as shown below:
BiocManager::install("GSEABase")
library(GSEABase)
# Load in GMT file of gene sets (we downloaded from the Broad Institute for MSigDB)
# not working
# c2 <- read.gmt("/data/c2.cp.v6.0.entrez.gmt.txt")
# msig <- GSEA(foldchanges, TERM2GENE=c2, verbose=FALSE)
# msig_df <- data.frame(msig)

# Resources
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
# g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi
# DAVID - http://david.abcc.ncifcrf.gov/tools.jsp
# clusterProfiler - http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# GeneMANIA - http://www.genemania.org/
# GenePattern - http://www.broadinstitute.org/cancer/software/genepattern/ (need to register)
# WebGestalt - http://bioinfo.vanderbilt.edu/webgestalt/ (need to register)
# AmiGO - http://amigo.geneontology.org/amigo
# ReviGO (visualizing GO analysis, input is GO terms) - http://revigo.irb.hr/
# WGCNA - https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
# GSEA - http://software.broadinstitute.org/gsea/index.jsp
# SPIA - https://www.bioconductor.org/packages/release/bioc/html/SPIA.html
# GAGE/Pathview - http://www.bioconductor.org/packages/release/bioc/html/gage.html





# M. ONO Omics 
# Gene Set Enrichment Analysis (fgsea)
# Load pathway data
# load('ono_t_cells/pathways')
# the data (immunologic signatures) can be downloaded from WEHI  http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata
# https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#C7
# Remove genes without Entrez ID from the DEG result
# Note: fgsea will use the object pathway which defines pathways by entrez ID.
# Use the annotation data org.Mm.eg.db for retrieving Entrez ID and gene name and append Entrez ID (and gene name) to the DESeq2 output, result.
# result$entrez <- mapIds(org.Mm.eg.db,
                        # keys = rownames(result),
                        # column = 'ENTREZID',
                        # keytype = 'SYMBOL',
                        # multiVals = 'first')
# add Gene name
# result$gene_name <- mapIds(org.Mm.eg.db,
                           # keys = rownames(result),
                           # column = 'GENENAME',
                           # keytype = 'SYMBOL',
                           # multiVals = 'first')
# fgsea_data <- result[!is.na(result$entrez),]
# Extract rank data
# ranks <- fgsea_data$log2FoldChange
# names(ranks) <- fgsea_data$entrez
# Perform GSEA
# fgsea_res <- fgseaMultilevel(pathways,  stats = ranks)
# ord <- order(fgsea_res$padj)
# fgsea_res <- fgsea_res[ord,]
# Examine GSEA result
# fgsea_res[1:3,]
# summary(fgsea_res)
# Remove redundant characters from pathway names
# fgsea_res$pathway <- sub(fgsea_res$pathway, pattern= 'HALLMARK_', replacement='')#
# Examine the change
# fgsea_res[1:3,]
# Show pathways with Normalised enrichment score (NES) greater than 0.5
# fgsea_plot <- fgsea_res[abs(fgsea_res$NES) > 0.5,]
# ord <- order(fgsea_plot$NES)
# fgsea_plot <- fgsea_plot[ord,]
# fp <- ggplot(fgsea_plot, mapping = aes(x = reorder(pathway, NES), y = NES)) +
  # geom_col(mapping = aes(fill = padj > 0.05)) +
  # coord_flip() +
  # labs(x = 'Pathway', y = 'Normalized Enrichment Score') +
  # theme(legend.position = "none")
# ggsave(filename = paste0('output4/',target,'_minus_', reference, "_fgsea1.pdf"), plot = fp, device = cairo_pdf, dpi = 300)
# Enrichment plot for significant pathways############
# logic <- fgsea_plot$padj < 0.05
# pathway_tmp <- as.vector(as.data.frame(fgsea_plot)[logic,1])
# how many pathways have been identified as sigificant?
# length(pathway_tmp)
# The above will give zero if there is no significant pathway
# To obtain entrez IDs for the significant pathways
# i <- 1
# lg <- grepl(pattern=pathway_tmp[i], names(pathways))
# pe <- plotEnrichment(pathways[lg][[1]], stats = ranks) + labs(title='')
# ggsave(filename = paste0('output4/',target,'_minus_', reference, '_', names(pathways[lg]),".pdf"), plot = pe, device = cairo_pdf, dpi = 300)

#Hierachical clustering and heatmap by key DEGs
# library(gplots)
# Aim to visualise the expression of immunologically important genes in the DEG list from your Venn diagram analysis
# Import immunologically important genes:
# key_genes_1 = read.table(file='/project/data/ono/Practical/data/key_genes_to_display.txt', sep=',')
# key_genes_1=unique(as.character(key_genes_1[,1]))
# Intersect of your DEGs and the immunological genes
# key_genes_1 = unique(intersect(genes,key_genes_1))#genes are DEGs in the first comparison.
# exprs.data=assays(vsd)[[1]]
# col.group=brewer.pal(10, "Set3")[as.factor(colData(vsd)$group)]
# select_var=names(sort(var_genes, decreasing=TRUE))[1:100]
# mat=exprs.data[intersect(key_genes_1, rownames(exprs.data)),]
# dim(mat) 
# Heatmap with hierarchical clustering
# pdf("output4/heatmap6.pdf")
# heatmap.2(as.matrix(mat), col=colorRampPalette(brewer.pal(10, "RdYlBu"))(256), trace="none",
          # srtCol = 45, margins = c(12,9), cexCol = 1.2,
          # density.info = "none", lhei = c(2, 8), keysize = 0.9,
          # main="Key immunological DEGs", scale="row", ColSideColors=col.group)
# dev.off()

# Check expression of some DEG genes
# by Gene expression plots (boxplots)
# custom_colors1 = c("dodgerblue", "darkblue", "darkorange", "firebrick2", "darkred", "blueviolet", "chocolate", "cadetblue", "aquamarine", "grey", "brown", "yellow")
# g='Tnfrsf4'
# tmp <- plotCounts(dds, gene = g, intgroup = "group", returnData = TRUE)
# tmp$group <- metadata_SP$group
# tmp$bim <- metadata_SP$bim
# tmp$timer <- metadata_SP$timer
# tmp$gene <- g
# gene_plots <- ggplot(data = tmp, aes(x = group, y = count, fill = group)) +
  # geom_boxplot(
    # size = 0.25,
    # outlier.size = 0.75,
    # outlier.alpha = 0.75
  # ) +
  # scale_y_log10() +
  # scale_fill_manual(values = custom_colors1) +
  # labs(
    # x = "",
    # y = "Normalised count"
  # ) +
  # theme_bw() +
  # ggtitle(g)
# ggsave(filename = paste0("output4/gene_expr_boxplots_", g, ".pdf"), plot = gene_plots, device = cairo_pdf, dpi = 300)

# Venn diagram analysis
# library(VennDiagram)
# To obtain the list of DEGs upregulated in DP Blue+ (in comparison to DP Negative) in heterozygotes.
# reference <- "SP_N_HET"
# target <- "SP_B_HET"
# res <- results(dds, contrast = c("group", target, reference))
# result <- as.data.frame(res)#convert to a data frame object
# lg1 <- !is.na(result$log2FoldChange)
# lg2 <- !is.na(result$padj)
# result <- result[lg1&lg2,]
# lg3 <- result$padj < 0.05
# lg4 <- abs(result$log2FoldChange) > 0
# result <- result[lg3&lg4,]
# SP_B_N_HET = rownames(result)
# here we aim to obtain Blue+ DP-specific genes that are regulated by Bim.
# To this end, create two more DEG lists by comparing
# Choose the reference and target groups
#(1) comparison between DP Blue+Red- and DP Timer negative from Bim KO mice
# reference <- "SP_N_KO"
# target <- "SP_B_KO"
# res <- results(dds, contrast = c("group", target, reference))
# result <- as.data.frame(res)#convert to a data frame object
# lg1 <- !is.na(result$log2FoldChange)
# lg2 <- !is.na(result$padj)
# result <- result[lg1&lg2,]
# lg3 <- result$padj < 0.05
# lg4 <- abs(result$log2FoldChange) > 0
# result <- result[lg3&lg4,]
# SP_B_N_KO = rownames(result)
# (2) comparison between DP Blue+Red+ and DP Blue+ from Bim KO mice
# reference <- "SP_B_KO"
# target <- "SP_BR_KO"
# res <- results(dds, contrast = c("group", target, reference))
# result <- as.data.frame(res)#convert to a data frame object
# lg1 <- !is.na(result$log2FoldChange)
# lg2 <- !is.na(result$padj)
# result <- result[lg1&lg2,]
# lg3 <- result$padj < 0.05
# lg4 <- abs(result$log2FoldChange) > 0
# result <- result[lg3&lg4,]
# SP_BR_B_KO = rownames(result)
# reference <- "SP_BR_KO"
# target <- "SP_R_KO"
# res <- results(dds, contrast = c("group", target, reference))
# result <- as.data.frame(res)#convert to a data frame object
# lg1 <- !is.na(result$log2FoldChange)
# lg2 <- !is.na(result$padj)
# result <- result[lg1&lg2,]
# lg3 <- result$padj < 0.05
# lg4 <- abs(result$log2FoldChange) > 0
# result <- result[lg3&lg4,]
# SP_BR_R_KO = rownames(result)
# To create a list object using the 3 DEG lists for Venn Diagram Analysis
# venn_items <- list(
 # SP_B_N_KO = SP_B_N_KO,
 # SP_BR_B_KO = SP_B_BR_KO,
 # SP_BR_R_KO = SP_BR_R_KO
# )

# Venn diagram plot
# custom_colors = c("dodgerblue", "darkblue", "darkorange", "firebrick2", "darkred")
# venn.diagram(
  # x = venn_items,
  # category.names = names(venn_items),
  # filename = paste0("output4/venn_diagram_KO.png"),
  # output = TRUE,
  # imagetype = "png",
  # lty = "blank",
  # fill = custom_colors[1:length(venn_items)],
  # cex = 0.75,
  # fontface = "bold",
  # cat.cex = 0.75,
  # cat.fontface = "bold",
  # cat.default.pos = "outer",
  # margin = 0.01
# )
# The function calculate.overlap will provide individual gene lists for all Venn diagram areas.
# venn.diagram.result <- calculate.overlap(venn_items)


## OVERVIEW
# summarized the steps in an analysis below:
# 1) Obtaining gene-level counts from Salmon using tximport
# Run tximport
# txi <- tximport(files, type="salmon", tx2gene=t2g, countsFromAbundance = "lengthScaledTPM")
# "files" is a vector wherein each element is the path to the salmon quant.sf file, and each element is named with the name of the sample.
# "t2g" is a 2 column data frame which contains transcript IDs mapped to geneIDs (in that order)
# 2) Creating the dds object:
# Check that the row names of the metadata equal the column names of the **raw counts** data
# all(colnames(txi$counts) == rownames(metadata))
# Create DESeq2Dataset object
# dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)
# 3) Exploratory data analysis (PCA & hierarchical clustering) - identifying outliers and sources of variation in the data:
# Transform counts for data visualization
# rld <- rlog(dds, blind=TRUE)
# Plot PCA 
# plotPCA(rld, intgroup="condition")
# Extract the rlog matrix from the object and compute pairwise correlation values
# rld_mat <- assay(rld)
# rld_cor <- cor(rld_mat)
# Plot heatmap
# pheatmap(rld_cor, annotation = metadata)
# 4) Run DESeq2:
# **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis in include other sources of variation using "dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ covaraite + condition)"
# Run DESeq2 differential expression analysis
# dds <- DESeq(dds)
# **Optional step** - Output normalized counts to save as a file to access outside RStudio using "normalized_counts <- counts(dds, normalized=TRUE)"
# 5) Check the fit of the dispersion estimates:
# Plot dispersion estimates
# plotDispEsts(dds)
# 6) Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:
# Specify contrast for comparison of interest
# contrast <- c("condition", "level_to_compare", "base_level")
# Output results of Wald test for contrast
# res <- results(dds, contrast = contrast, alpha = 0.05)
# Shrink the log2 fold changes to be more accurate
# res <- lfcShrink(dds, coef = "sampletype_group1_vs_group2", type = "apeglm")	 
# The coef will be dependent on what your contras was. and should be identical to what is stored in resultsNames()
# 7) Output significant results:
# Set thresholds
# padj.cutoff < - 0.05
# Turn the results object into a tibble for use with tidyverse functions
# res_tbl <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
# Subset the significant results
# sig_res <- dplyr::filter(res_tbl, padj < padj.cutoff)
# 8) Visualize results: volcano plots, heatmaps, normalized counts plots of top genes, etc.
# 9) Perform analysis to extract functional significance of results: GO or KEGG enrichment, GSEA, etc.
# 10) Make sure to output the versions of all tools used in the DE analysis:
# sessionInfo()
# For better reproducibility, it can help to create RMarkdown reports

