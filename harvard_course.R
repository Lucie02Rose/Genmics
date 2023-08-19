#22/07/2023, useful info from R book for bio
#delete prev stuff
rm(list = ls())
#import and install packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(purrr)  # Load the purrr

install.packages("knitr")
install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("pheatmap")
library(knitr)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db) 
library(pathview)
library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v75)

#dplyr can select, slice, arrange, 
#mutate, filter data, slice rows

#fetching files
getwd()
setwd("C:/Users/Acer/Desktop/FILES_LONDON_07_22/R_course/code")
getwd()

species <- c("ecoli", "human", "corn")
glengths <- c(4.6, 3000, 50000)
combined <- c(species, glengths)
combined
df <- data.frame(species, glengths)
df

expression <- c("low", "high", "medium", "high", "low", "medium", "high")
expression <- factor(expression)
expression

samplegroup <- c("OE","CTL","KO","CTL","OE","OE","CTL","KO","KO")
samplegroup <- factor(samplegroup)
samplegroup

titles <- c("Catch-22", "Pride and Prejudice", "Nineteen Eighty Four")
pages <- c(453, 432, 328)
df1 <- data.frame(titles, pages)
df1

list1 <- list(species, df, species)
list1

temp_conv <- function(temp_f) {
  temp_c = (temp_f - 32) * 5 / 9; temp_k = temp_c + 273.15
  return(round(temp_k))
}
temp_conv(70)

age <- c(15, 22, 45, 52, 73, 81)
idx_num <- which(age > 50 | age < 18)
age[idx_num]
idx <- which(metadata$celltype == "typeA")
sub_meta <- metadata[which(metadata$replicate > 1), ]
#[row, column]
A %in% B
any(A %in% B)
all(A %in% B)

first <- c("A","B","C","D","E")
second <- c("B","D","E","A","C")  # same letters but different order
match(first,second)
reorder_idx <- match(first,second) 
second_reordered <- second[reorder_idx] 

#samplemeans <- map_dbl(rpkm_ordered, mean) 
#map() creates a list.
#map_lgl() creates a logical vector.
#map_int() creates an integer vector.
#map_dbl() creates a "double" or numeric vector.
#map_chr() creates a character vector.


#ggplot
#points (geom_point, geom_jitter for scatter plots, dot plots, etc)
#lines (geom_line, for time series, trend lines, etc)
#boxplot (geom_boxplot, for, well, boxplots!)

#ggplot(new_metadata) +
#geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
               #shape=celltype), size=2.25) +
  #theme_bw() +
  #theme(axis.title = element_text(size=rel(1.5)))	) 

#personal_theme <- function(){
#theme_bw() +
  #theme(axis.title=element_text(size=rel(1.5))) +
  #theme(plot.title=element_text(size=rel(1.5), hjust=0.5))
#}

#boxplot
#scale_color_manual(values=c("purple","orange")) with scale_fill_manual(values=c("purple","orange")).
#ggplot(new_metadata) +
#geom_point(aes(x=age_in_days, y=samplemeans, color=genotype, shape=celltype), size=rel(3.0)) +
  #xlab("Age (days)") +
  #ylab("Mean expression") +
  #ggtitle("Expression with Age") +
  #personal_theme()

#write(glengths, file="data/genome_lengths.txt", ncolumns = 1)
#pdf("figures/scatterplot.pdf")
### Closing the device is essential to save the temporary file created by pdf()/png()
#dev.off()

sessionInfo()
### Running more than one command with piping - tidyverse
sqrt(83) %>% round(digits = 2)
#rownames_to_column() and column_to_rownames()
#as_tibble(name_of_df)

# Read in the functional analysis results
#functional_GO_results <- read_delim(file = "data/gprofiler_results_Mov10oe.tsv", delim = "\t" )

#bp_oe <- functional_GO_results %>%
#filter(domain == "BP")

## Selecting columns to keep
#bp_oe <- bp_oe %>%
  #select(term.id, term.name, p.value, query.size, term.size, overlap.size, intersection)

## Selecting columns to remove
#bp_oe <- bp_oe %>%
  #select(-c(query.number, significant, recall, precision, subgraph.number, relative.depth, domain))

## Order by adjusted p-value ascending
#bp_oe <- bp_oe %>%
  #arrange(p.value)
#descending would be arrange(desc(p.value))

# Provide better names for columns
#bp_oe <- bp_oe %>% 
  #dplyr::rename(GO_id = term.id, 
                #GO_term = term.name)

# Create gene ratio column based on other columns in dataset
#bp_oe <- bp_oe %>%
  #mutate(gene_ratio = overlap.size / query.size)
