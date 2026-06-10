# first read all library 
library(tidyverse)
library(DESeq2)
library(airway)
library("tximport")

# setup the data 

count_data <- read.csv(file.choose('C:/Users/Mahima/Downloads/DESeq/GSE329380_count_.csv'))
head(count_data)

# READ IN SIMPLE INFO
colData <-read.csv(file.choose('C:/Users/Mahima/Downloads/SraRunTable.csv'))
#if you delet some col.
metadata <- colData %>%
  select(1,32) %>%
  head()
dim(metadata)
str(metadata)
metadata <- colData %>%
  dplyr::select(1, 32)

#making sure to all coldata and count data having matched
all(colnames(count_data) %in% rownames(metadata))

head(colnames(count_data))
head(rownames(metadata))

rownames(metadata) <- metadata$Run

all(colnames(count_data) %in% rownames(metadata))
head(colnames(count_data))
head(rownames(metadata))
head(metadata$Run)
#HERE YOU SEE SOME COL. DIDN'T MATCH so i put some to to rewrite there name 

tail(colnames(count_data))
dim(count_data)
dim(metadata)

# Keep first two columns unchanged (X and gene_id)
colnames(count_data)[3:ncol(count_data)] <- metadata$Run
head(colnames(count_data))

ncol(count_data)
length(metadata$Run)
#add here matrix txt
colData <-read.csv(file.choose('C:/Users/Mahima/Downloads/DESeq/GSE329380_matrix.txt'))


# here my some data is duplicate so i remove all duplicates here 

# Set gene IDs as row names
rownames(count_data_unique) <- count_data_unique$gene_id

# Remove X and gene_id columns
count_matrix <- count_data_unique[, -c(1,2)]

# Convert to matrix
count_matrix <- as.matrix(count_matrix)

# Convert to integer counts
mode(count_matrix) <- "integer"
dim(metadata)
dim(count_matrix)

rownames(metadata) <- metadata$Run

metadata <- metadata[colnames(count_matrix), , drop = FALSE]

ncol(count_matrix)
nrow(metadata)

#next step to construct theDESeq data objects

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ tissue)
#pre filtering thats having low gene count 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds
# set the factor level
dds$tissue <- relevel(dds$tissue , ref = "non-small cell lung cancer tumor tissue")

# run DESeq ------------------------

dds <- DESeq(dds)

res <- results(dds)

res
# Exploring result
summary(res)
res0.5 <- results(dds , alpha = 0.5)
summary(res0.5)

#MA plotting 
plotMA(res)

