library(DESeq2)
library(tidyverse)

# Load the txi object
txi <- readRDS("RObjects/txi.rds")

# Load sample info table

sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", 
                       col_types = "cccc")

all(colnames(txi$counts) == sampleinfo$SampleName)

# The simple model

simple.model <- as.formula(~ TimePoint)  

# The model matrix

model.matrix(simple.model, data = sampleinfo)

# Exercise 1

simple.model <- as.formula(~ Status)

model.matrix(simple.model, data = sampleinfo)

# Switching "uninfected" to be the intercept

sampleinfo <- mutate(sampleinfo, 
        Status = fct_relevel(Status, "Uninfected"))

model.matrix(simple.model, data = sampleinfo)

# Building a DESeq2DataSet

simple.model <- as.formula(~ Status)
sampleinfo <- mutate(sampleinfo, 
       Status = fct_relevel(Status, "Uninfected"))

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)

# Filter out the unexpressed

keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep, ]

# Differential expression analysis

## Deseq2 workflow

### Estimate Size Factors

ddsObj <- estimateSizeFactors(ddsObj.filt) 

normalizationFactors(ddsObj.filt)

logcounts <- log2(counts(ddsObj, normalized = FALSE) + 1)

limma::plotMA(logcounts, array = 5, ylim = c(-5, 5))
abline(h = 0, col = "red")

logNormalizedCounts <- log2(counts(ddsObj,
                                   normalized = TRUE) + 1)

limma::plotMA(logNormalizedCounts, 
              array = 5, 
              ylim = c(-5, 5))
abline(h = 0, col = "red")

### Estimate Dispersions

ddsObj <- estimateDispersions(ddsObj)
plotDispEsts(ddsObj)

### Modelling and Wald test

ddsObj <- nbinomWaldTest(ddsObj)


# The DESeq command

ddsObj <- DESeq(ddsObj.filt)

# Generate a table of differential expression results

results.simple <- results(ddsObj, alpha = 0.05)
results.simple

# Exercise 2

## a) Upregulation

sum(results.simple$log2FoldChange > 0 & results.simple$padj < 0.05, 
    na.rm = TRUE)

## b) Downregulation 

sum(results.simple$log2FoldChange < 0 & results.simple$padj < 0.05, 
    na.rm = TRUE)

# The Additive model

additive.model <- as.formula(~ TimePoint + Status)

# Exercise 3

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = additive.model)

ddsObj.filt <- ddsObj.raw[keep, ]
design(ddsObj.filt) <- additive.model

## 1. Run the DESeq2 workflow

ddsObj <- DESeq(ddsObj.filt)

## 2. extract results

results.additive <- results(ddsObj, alpha = 0.05)

### a) How many coefficients

model.matrix(additive.model, data = sampleinfo)

### d) What contrast does `results.additive`

results.additive

### e) How many genes with padj < 0.05

sum(results.additive$padj < 0.05, na.rm = TRUE)

# How `results` picks default contrast

resultsNames(ddsObj)

# Rename results

results.InfectedvUninfected <- results.additive
rm(results.additive)

# Getting top 100 genes

topGenesIvU <- as.data.frame(results.InfectedvUninfected) %>%
  rownames_to_column("GeneID") %>%
  top_n(100, wt = -padj)

# Exercise 5

## 1. d33 v d11

results.d33vd11 <- results(ddsObj,
                           name = "TimePoint_d33_vs_d11",
                           alpha = 0.05)

## 2. How many sig genes

sum(results.d33vd11$padj < 0.05, na.rm = TRUE)

# Look at exploratory data analysis to determine if we want
# an interaction

vstcounts <- vst(ddsObj.raw, blind = TRUE)
plotPCA(vstcounts, intgroup = c("Status", "TimePoint"))

# Exercise 3

## 1 . Build model

interaction.model <- ~ TimePoint + Status + TimePoint:Status

## 2. Run DESeq2

design(ddsObj.filt) <- interaction.model
ddsObj.interaction <- DESeq(ddsObj.filt)

## 3. Extract results

results.interaction <- results(ddsObj.interaction, alpha = 0.05)
sum(results.interaction$padj < 0.05, na.rm = TRUE)

# Extract a specific contrast from an interaction model

resultsNames(ddsObj.interaction)

results.interaction.d11 <- results(ddsObj.interaction,
                 name = "Status_Infected_vs_Uninfected",
                 alpha = 0.05)


results.interaction.d33 <- results(ddsObj.interaction,
                  contrast = list(c("Status_Infected_vs_Uninfected",
                                    "TimePointd33.StatusInfected")),
                  alpha = 0.05)

sum(results.interaction.d11$padj < 0.05, na.rm = TRUE)
sum(results.interaction.d33$padj < 0.05, na.rm = TRUE)

# Exercise 7

## 1. d33 v d11 - Infected

results.interaction.Inf <- results(ddsObj.interaction, 
                      contrast = list(c("TimePoint_d33_vs_d11",
                                        "TimePointd33.StatusInfected")),
                      alpha = 0.05)
sum(results.interaction.Inf$padj < 0.05, na.rm = TRUE)

## 2. d33 v d11 - Uninfected

results.interaction.Uninf <- results(ddsObj.interaction,
                                     name = "TimePoint_d33_vs_d11",
                                     alpha = 0.05)
sum(results.interaction.Uninf$padj < 0.05, na.rm = TRUE)




