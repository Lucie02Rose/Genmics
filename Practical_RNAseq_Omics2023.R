
#########################################
#OMICS: RNA-seq practical sample code
#Computational Omics, Year 2, Biology
#Imperial College London
#Dr Masahiro Ono    12 May 2023
########################################
library(AnnotationDbi)
library(DESeq2) # differential expression analysis
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(org.Mm.eg.db)
library(fgsea) # gene set enrichment analysis
library(ReactomePA) # reactome pathway analysis
library(VennDiagram) # venn diagram analysis

# Set the working directory ###############
#
#See the instruction
#
###########################################
# load count data
load('/project/data/ono/Tutorial/data/countdata')
metadata_DP <- read.table(file = '/project/data/ono/Tutorial/data/metadata_DP.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE, row.names=1)
metadata_SP <- read.table(file = '/project/data/ono/Tutorial/data/metadata_SP.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE, row.names=1)

colnames(metadata_SP)=colnames(metadata_DP)=c("sample", "group","cell.population","Timer","Bim") #re-name the column names for clarity

#Inspect loaded data
dim(countdata)
countdata[1:10,1:6]
summary(countdata)
dim(metadata_DP)
dim(metadata_SP)

# Choose DP sample data from countdata
counts_SP <- countdata[,metadata_SP[,1]]

#Inspect filtered data
dim(counts_SP)
counts_DP[1:10,1:6]
summary(counts_SP)

# Filter out genes which are barely detected:
logic <- rowSums(counts_DP) > 10
counts_SP <- counts_DP[logic,]

##DESeq2####################################
library(DESeq2) # differential expression analysis
# Run DESeq2 using the count data##########
dds1 <- DESeqDataSetFromMatrix(countData = counts_sp, colData = metadata_DP, design = ~ group)

#Examine dds##########
dds1
dim(dds1)
class(dds1)
colData(dds1) #

# Principal Component Analysis (PCA)#################
#a fast and easy way to perform PCA by DESeq2
#Obtain a gene expression matrix

dir.create("output4")#create a new folder for output files

#Obtain normalised gene expression data
vsd1 <- vst(dds, blind = FALSE)

#Perform PCA using plotPCA() of DESeq2 for visualisation
library(ggplot2)#package for visualisation
pca_plot <- plotPCA(vsd, intgroup = c('group'), returnData = FALSE)
ggsave(filename = 'output4/SP_pca_plot.pdf', plot = pca_plot, device = cairo_pdf, dpi = 300)

#Obtain a gene expression matrix
exprs.data=t(assay(vsd)) #note that the matrix is transposed

#Inspect exprs.data
exprs.data[1:3,1:3]

#Perform PCA by the function prcomp()
pca_prcomp = prcomp(exprs.data, scale=FALSE)

#Plot PCA result
pdf("output4/PCA3.pdf")
par(mfrow=c(2,2))
plot(pca_prcomp$x[,c(2,3)], main='Sample plot', col=colData(vsd)$group, pch=19)
text(pca_prcomp$x[,c(2,3)], labels=rownames(pca_prcomp$x[,c(2,3)]), pos=4, cex=0.8)
abline(v=0, h=0, col=8)

plot( pca_prcomp$rotation[,c(2,3)], main='Gene loading')
abline(v=0, h=0, col=8)

eig = (pca_prcomp$sdev)^2; names(eig) = 1:length(eig); eig= 100* eig/sum(eig)
barplot(eig, main='Scree plot', xlab= "PC", ylab = 'Percentage of explained variances')
dev.off()

#Your task: Plot PC1-PC3 and PC2-PC3 as well


#####################################
##Hierachical clustering and heatmap
##
library(gplots)
library(RColorBrewer)

#Calculate the variance of each gene expression
exprs.data=assays(vsd)[[1]]
var_genes=apply(exprs.data, 1, var)
col.group=brewer.pal(10, "Set3")[as.factor(colData(vsd)$group)]

# Get the names for the top 100 most highly variable genes and subset data
select_var=names(sort(var_genes, decreasing=TRUE))[1:100]
hi_var=exprs.data[select_var,]
dim(hi_var) # Get gene number and sample number as output

pdf("output4/heatmap2.pdf")
heatmap.2(as.matrix(hi_var), col=colorRampPalette(brewer.pal(10, "RdYlBu"))(256), trace="none",
srtCol = 45, margins = c(12,9), cexCol = 1.2,
density.info = "none", lhei = c(2, 8), keysize = 0.9,
main="Top 100 Most Variable Genes", scale="row", ColSideColors=col.group)
dev.off()



######################################
# Run DESeq2 with several comparisons

#Estimation of size factors and dispersion
dds <- DESeq(dds)

###Differentially expressed gene (DEG) analysis
#
#Hypothesis: Nr4a3-Timer Blue+ cells capture DP cells under negative selection
#
#Choose the reference sample group - use the exact name of a sample group
reference <- "DP_N_HET"
#Choose the target group
target <- "DP_B_HET"

##DESeq2 will calculate log2 FC target - reference and perform a Walt test##############

res <- results(dds, contrast = c("group", target, reference))

##Extract DEGs from res##############################
result <- as.data.frame(res)#convert to a data frame object

##Remove NA
lg1 <- !is.na(result$log2FoldChange)
lg2 <- !is.na(result$padj)
result <- result[lg1&lg2,]

##Filter results by padj (adjusted p-value) and log2 fold change (absolute level)
lg3 <- result$padj < 0.05
lg4 <- abs(result$log2FoldChange) > 0

result <- result[lg3&lg4,]

##Order the result according to padj
ord <- order(result$padj, decreasing = FALSE)
result <- result[ord,]

#Inspect
result[1:3,]
summary(result)



######################################
#Volcano plot
library(ggplot2)
library(ggrepel)
#Define the colour of genes by adjusted p-value and log2 fold change
test = ifelse(result$padj < 0.01, TRUE, FALSE)
upregulated = ifelse(test = result$log2FoldChange > 1, TRUE, FALSE)
downregulated = ifelse(test = result$log2FoldChange < -1, TRUE, FALSE)

#ifelse is used to define the colour - check the logical operators
state = ifelse(test == TRUE & upregulated == TRUE, 'Upregulated', ifelse(test == TRUE & downregulated == TRUE, 'Downregulated','Not statistically significant'))

genes = rownames(result)[test & (upregulated|downregulated)]

#Define which gene names are shown
point_label = rownames(result)
point_label[!(rownames(result) %in% genes)] = ""

volcano_plot <- ggplot(data = result, mapping = aes(x = log2FoldChange, y = -log10(padj))) +
geom_point(mapping = aes(color = state, alpha = 0.9), size = 2) +
geom_text_repel(mapping = aes(label = point_label), max.overlaps=20) +
guides(alpha = "none") +
labs(x = 'log2(fold change)',
y = '-log10(adjusted p-value)') +
theme_minimal() +
ggtitle(paste0(target,' minus ', reference))


ggsave(filename = paste0('output4/',target,'_minus_', reference, "_volcanoplot_SP_B_N.pdf"), plot = volcano_plot, device = cairo_pdf, dpi = 300)

###Gene Set Enrichment Analysis (fgsea)###############
library(fgsea) # gene set enrichment analysis

##Load pathway data##########
load('ono_t_cells/pathways')
#the data (immunologic signatures) can be downloaded from WEHI  http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata
#https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#C7
#
#Remove genes without Entrez ID from the DEG result
#Note: fgsea will use the object pathway which defines pathways by entrez ID.
#####################################
###Use the annotation data org.Mm.eg.db for retrieving Entrez ID and gene name and append Entrez ID (and gene name) to the DESeq2 output, result.

library(org.Mm.eg.db)
result$entrez <- mapIds(org.Mm.eg.db,
keys = rownames(result),
column = 'ENTREZID',
keytype = 'SYMBOL',
multiVals = 'first')

##add Gene name
result$gene_name <- mapIds(org.Mm.eg.db,
keys = rownames(result),
column = 'GENENAME',
keytype = 'SYMBOL',
multiVals = 'first')

fgsea_data <- result[!is.na(result$entrez),]

#Extract rank data
ranks <- fgsea_data$log2FoldChange
names(ranks) <- fgsea_data$entrez

#Perform GSEA
fgsea_res <- fgseaMultilevel(pathways,  stats = ranks)

ord <- order(fgsea_res$padj)
fgsea_res <- fgsea_res[ord,]

#Examine GSEA result
fgsea_res[1:3,]
summary(fgsea_res)

#Remove redundant characters from pathway names
fgsea_res$pathway <- sub(fgsea_res$pathway, pattern= 'HALLMARK_', replacement='')#

#Examine the change
fgsea_res[1:3,]

#Show pathways with Normalised enrichment score (NES) greater than 0.5

fgsea_plot <- fgsea_res[abs(fgsea_res$NES) > 0.5,]
ord <- order(fgsea_plot$NES)
fgsea_plot <- fgsea_plot[ord,]

fp <- ggplot(fgsea_plot, mapping = aes(x = reorder(pathway, NES), y = NES)) +
geom_col(mapping = aes(fill = padj > 0.05)) +
coord_flip() +
labs(x = 'Pathway', y = 'Normalized Enrichment Score') +
theme(legend.position = "none")

ggsave(filename = paste0('output4/',target,'_minus_', reference, "_fgsea1.pdf"), plot = fp, device = cairo_pdf, dpi = 300)



#Enrichment plot for significant pathways############
logic <- fgsea_plot$padj < 0.05

pathway_tmp <- as.vector(as.data.frame(fgsea_plot)[logic,1])

#how many pathways have been identified as sigificant?
length(pathway_tmp)

#The above will give zero if there is no significant pathway
#To obtain entrez IDs for the significant pathways

i <- 1
lg <- grepl(pattern=pathway_tmp[i], names(pathways))
pe <- plotEnrichment(pathways[lg][[1]], stats = ranks) + labs(title='')

ggsave(filename = paste0('output4/',target,'_minus_', reference, '_', names(pathways[lg]),".pdf"), plot = pe, device = cairo_pdf, dpi = 300)

#####################################
##Hierachical clustering and heatmap by key DEGs
##
library(gplots)

##Aim to visualise the expression of immunologically important genes in the DEG list from your Venn diagram analysis

##Import immunologically important genes:
key_genes_1 = read.table(file='/project/data/ono/Practical/data/key_genes_to_display.txt', sep=',')
key_genes_1=unique(as.character(key_genes_1[,1]))

##Intersect of your DEGs and the immunological genes
key_genes_1 = unique(intersect(genes,key_genes_1))#genes are DEGs in the first comparison.

exprs.data=assays(vsd)[[1]]
col.group=brewer.pal(10, "Set3")[as.factor(colData(vsd)$group)]

select_var=names(sort(var_genes, decreasing=TRUE))[1:100]
mat=exprs.data[intersect(key_genes_1, rownames(exprs.data)),]
dim(mat) #

#Heatmap with hierarchical clustering
pdf("output4/heatmap6.pdf")
heatmap.2(as.matrix(mat), col=colorRampPalette(brewer.pal(10, "RdYlBu"))(256), trace="none",
srtCol = 45, margins = c(12,9), cexCol = 1.2,
density.info = "none", lhei = c(2, 8), keysize = 0.9,
main="Key immunological DEGs", scale="row", ColSideColors=col.group)
dev.off()

##########################################
#Check expression of some DEG genes
#by Gene expression plots (boxplots)
custom_colors1 = c("dodgerblue", "darkblue", "darkorange", "firebrick2", "darkred", "blueviolet", "chocolate", "cadetblue", "aquamarine", "grey", "brown", "yellow")

g='Tnfrsf4'
tmp <- plotCounts(dds, gene = g, intgroup = "group", returnData = TRUE)
tmp$group <- metadata_SP$group
tmp$bim <- metadata_SP$bim
tmp$timer <- metadata_SP$timer
tmp$gene <- g

#
gene_plots <- ggplot(data = tmp, aes(x = group, y = count, fill = group)) +
geom_boxplot(
size = 0.25,
outlier.size = 0.75,
outlier.alpha = 0.75
) +
scale_y_log10() +
scale_fill_manual(values = custom_colors1) +
labs(
x = "",
y = "Normalised count"
) +
theme_bw() +
ggtitle(g)

ggsave(filename = paste0("output4/gene_expr_boxplots_", g, ".pdf"), plot = gene_plots, device = cairo_pdf, dpi = 300)



###############################
# Venn diagram analysis
library(VennDiagram)

##To obtain the list of DEGs upregulated in DP Blue+ (in comparison to DP Negative) in heterozygotes.
reference <- "SP_N_HET"
target <- "SP_B_HET"
res <- results(dds, contrast = c("group", target, reference))
result <- as.data.frame(res)#convert to a data frame object
lg1 <- !is.na(result$log2FoldChange)
lg2 <- !is.na(result$padj)
result <- result[lg1&lg2,]
lg3 <- result$padj < 0.05
lg4 <- abs(result$log2FoldChange) > 0
result <- result[lg3&lg4,]
SP_B_N_HET = rownames(result)

#Here we aim to obtain Blue+ DP-specific genes that are regulated by Bim.
#To this end, create two more DEG lists by comparing

#Choose the reference and target groups


reference <- "SP_B_HET"
target <- "SP_BR_HET"
res <- results(dds, contrast = c("group", target, reference))
result <- as.data.frame(res)#convert to a data frame object
lg1 <- !is.na(result$log2FoldChange)
lg2 <- !is.na(result$padj)
result <- result[lg1&lg2,]
lg3 <- result$padj < 0.05
lg4 <- abs(result$log2FoldChange) > 0
result <- result[lg3&lg4,]
SP_B_BR_HET = rownames(result)

reference <- "SP_BR_HET"
target <- "SP_R_HET"
res <- results(dds, contrast = c("group", target, reference))
result <- as.data.frame(res)#convert to a data frame object
lg1 <- !is.na(result$log2FoldChange)
lg2 <- !is.na(result$padj)
result <- result[lg1&lg2,]
lg3 <- result$padj < 0.05
lg4 <- abs(result$log2FoldChange) > 0
result <- result[lg3&lg4,]
SP_BR_R_HET = rownames(result)

#(1) comparison between DP Blue+Red- and DP Timer negative from Bim KO mice
reference <- "SP_N_KO"
target <- "SP_B_KO"
res <- results(dds, contrast = c("group", target, reference))
result <- as.data.frame(res)#convert to a data frame object
lg1 <- !is.na(result$log2FoldChange)
lg2 <- !is.na(result$padj)
result <- result[lg1&lg2,]
lg3 <- result$padj < 0.05
lg4 <- abs(result$log2FoldChange) > 0
result <- result[lg3&lg4,]
SP_B_N_KO = rownames(result)

# (2) comparison between DP Blue+Red+ and DP Blue+ from Bim KO mice
reference <- "SP_B_KO"
target <- "SP_BR_KO"
res <- results(dds, contrast = c("group", target, reference))
result <- as.data.frame(res)#convert to a data frame object
lg1 <- !is.na(result$log2FoldChange)
lg2 <- !is.na(result$padj)
result <- result[lg1&lg2,]
lg3 <- result$padj < 0.05
lg4 <- abs(result$log2FoldChange) > 0
result <- result[lg3&lg4,]
SP_BR_B_KO = rownames(result)

reference <- "SP_BR_KO"
target <- "SP_R_KO"
res <- results(dds, contrast = c("group", target, reference))
result <- as.data.frame(res)#convert to a data frame object
lg1 <- !is.na(result$log2FoldChange)
lg2 <- !is.na(result$padj)
result <- result[lg1&lg2,]
lg3 <- result$padj < 0.05
lg4 <- abs(result$log2FoldChange) > 0
result <- result[lg3&lg4,]
SP_BR_R_KO = rownames(result)

##To create a list object using the 3 DEG lists for Venn Diagram Analysis
venn_items <- list(
SP_B_N_KO = SP_B_N_KO,
SP_BR_B_KO = SP_B_BR_KO,
SP_BR_R_KO = SP_BR_R_KO
)

###Venn diagram plot############
custom_colors = c("dodgerblue", "darkblue", "darkorange", "firebrick2", "darkred")

venn.diagram(
  x = venn_items,
  category.names = names(venn_items),
  filename = paste0("output4/venn_diagram_KO.png"),
  output = TRUE,
  imagetype = "png",
  lty = "blank",
  fill = custom_colors[1:length(venn_items)],
  cex = 0.75,
  fontface = "bold",
  cat.cex = 0.75,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  margin = 0.01
)

#The function calculate.overlap will provide individual gene lists for all Venn diagram areas.
venn.diagram.result <- calculate.overlap(venn_items)


