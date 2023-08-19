###############
##OMICS: RNA-seq practical
##Tutorial 2023
#Imperial College London
#Dr Masahiro Ono    12 May 2023
#
# Create and set the working directory ###############
mkdir tmp
cd tmp

ls
ls /project/data/ono/Tutorial/data/

#Inspect fastq file
zcat /project/data/ono/Tutorial/data/DP_B_Het_R1_S1_L001_R1_001.fastq.gz | head -10
#Mac users may need to do the following instead
#zcat < /project/data/ono/Tutorial/data/DP_B_Het_R1_S1_L001_R1_001.fastq.gz |head -10


#check other files as well

###Load modules
module load R/ono
module load salmon

##Mus_musculus.GRCm38.cdna.all_plusTimer.fa.gz is the transcriptome file to be used by salmon indexing
#Samon indexing
#https://salmon.readthedocs.io/en/latest/salmon.html
salmon index -t /project/data/ono/Tutorial/data/Mus_musculus.GRCm38.cdna.all_plusTimer.fa.gz -i transcripts_index --type quasi -k 31

#Check the output folder transcripts_index

#Quantifying in mapping-based mode
#data are paired-end
salmon quant -i transcripts_index -l A --validateMappings -1 /project/data/ono/Tutorial/data/DP_B_Het_R1_S1_L001_R1_001.fastq.gz  -2 /project/data/ono/Tutorial/data/DP_B_Het_R1_S1_L001_R2_001.fastq.gz -o DP_B_Het_R1


salmon quant -i transcripts_index -l A --validateMappings -1 /project/data/ono/Tutorial/data/DP_N_Het_R1_S2_L001_R1_001.fastq.gz  -2 /project/data/ono/Tutorial/data/DP_N_Het_R1_S2_L001_R2_001.fastq.gz -o DP_N_Het_R1


#Inspect output files
head DP_B_Het_R1/quant.sf
head DP_N_Het_R1/quant.sf

# R   ######################################################################################################
#
R

library(gplots)
library(RColorBrewer)
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(tximport)

#Save R output files in the new folder output#########
dir.create('output')

#Load annotation data
load('/project/data/ono/Tutorial/data/tx2gene')

#Inspect tx2gene
class(tx2gene)
tx2gene[1:30,]

#
files <- file.path(".", c('DP_B_Het_R1', 'DP_N_Het_R1'), "quant.sf")
names(files) <- c('DP_B_Het_R1', 'DP_N_Het_R1')

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion =TRUE)
head(txi.salmon$counts)

##################################################################################################################

# load count data
counts <- read.table(file = '/project/data/ono/Tutorial/data/counts.csv', header = TRUE, sep = ',')
#Import meta data
metadata <- read.table(file = '/project/data/ono/Tutorial/data/metadata.csv', header = TRUE, sep = ',')

# Filter out genes which were barely detected:
logic <- rowSums(counts) > 10
counts <- counts[logic,]

#Inspect loaded data
dim(counts)
counts[1:10,1:6]
summary(counts)
dim(metadata)
metadata[1:4,1:3]

##DESeq2#######################################################################################################
library(DESeq2) # differential expression analysis
# Run DESeq2 using the count data##########
dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ sample)

#Examine dds##########
dds
dim(dds)
class(dds)
colData(dds) #

# Get the normalized count matrix##############
dds <- DESeq(dds)

# Principal Component Analysis (PCA)###############################################################################
#a fast and easy way to perform PCA by DESeq2
#Obtain a gene expression matrix
vsd <- vst(dds, blind = FALSE)

library(ggplot2)
pca_plot <- plotPCA(vsd, intgroup = c('sample'))
ggsave(filename = 'output/pca_plot.pdf', plot = pca_plot, device = cairo_pdf, dpi = 300)

#Obtain a gene expression matrix
exprs.data=assay(vsd)

#Inspect exprs.data
exprs.data[1:3,1:3]

#Perform PCA by the function prcomp()
pca_prcomp = prcomp(t(exprs.data), scale=FALSE)

#Inspect the prcomp output
summary(pca_prcomp$x)
summary(pca_prcomp$rotation)
summary(pca_prcomp$sdev)

#Plot PCA result
pdf("output/PCA.pdf")
par(mfrow=c(2,2))
plot(pca_prcomp$x[,1:2], main='Sample plot', col=colData(vsd)$sample, pch=19)
text(pca_prcomp$x[,1:2], labels=rownames(pca_prcomp$x[,1:2]), pos=4, cex=0.8)
abline(v=0, h=0, col=8)

plot( pca_prcomp$rotation[,1:2], main='Gene loading')
abline(v=0, h=0, col=8)

eig = (pca_prcomp$sdev)^2; names(eig) = 1:length(eig); eig= 100* eig/sum(eig)
barplot(eig, main='Scree plot', xlab= "PC", ylab = 'Percentage of explained variances')
dev.off()

#Your task: Plot PC1-PC3 and PC2-PC3 as well


####################################################################
##Hierachical clustering and heatmap ##
library(gplots)
library(RColorBrewer)

exprs.data <- assay(vsd)
## Obtain pearson correlation between samples
cor.data = cor(exprs.data, method='pearson')

## Visualise sample correlations by hierarchical clustering and heatmap
col.group=rainbow(16)[as.factor(colData(vsd)$sample)]#Convert groups into colours

pdf('output/correlation.pdf')
heatmap.2(as.matrix(cor.data), col=colorRampPalette(brewer.pal(10, "RdYlBu"))(256), trace="none",
srtCol = 45, margins = c(12,9), cexCol = 0.8,
density.info = "none", lhei = c(2, 8), keysize = 0.9,
main="Sample correlations", scale="none", ColSideColors=col.group)
dev.off()



