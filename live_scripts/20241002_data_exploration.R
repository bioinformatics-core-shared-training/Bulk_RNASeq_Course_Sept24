library(tximport) 
library(DESeq2)
library(tidyverse)

#Load our data

#Read in metadata
sampleinfo <- read_tsv("data/samplesheet.tsv", col_types = c("cccc"))
arrange(sampleinfo, Status, TimePoint, Replicate)

#creating vector for all filepaths
files <- file.path("salmon", sampleinfo$SampleName, "quant.sf")
files <- set_names(files, sampleinfo$SampleName)
tx2gene <- read_tsv("references/tx2gene.tsv")
head(tx2gene)

#Loading in salmon count data
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
str(txi)
head(txi$counts)

#Save txi count matrix as an R object
saveRDS(txi, file='salmon_outputs/txi.rds')

#Create a raw counts matrix
rawCounts <- round(txi$counts, 0)

#Check dimensions of the matrix
dim(rawCounts)

#create logical vector for filtering lowly expressed genes
keep <- rowSums(rawCounts) > 5
#create a new table with filtered genes
table(keep, useNA="always")
filtCounts <- rawCounts[keep,]

#check dimensions of new table
dim(filtCounts)

summary(filtCounts)

#boxplot of raw counts
boxplot(filtCounts, main="Raw Counts", las = 2)

#plot standard deviation against mean counts
plot(rowMeans(filtCounts), rowSds(filtCounts),
     main="Raw counts: sd against mean",
     xlim = c(0, 10000),
     ylim = c(0, 5000))

#Data Transformation
#log2 transformation
logCounts <- log2(filtCounts + 1)

#boxplot for log transformed data
boxplot(logCounts, main="Log2 Counts", las = 2)

#rlog transformation
rlogcounts <- rlog(filtCounts)
boxplot(rlogcounts, main ="rlog counts", las = 2)

#PCA
library(ggfortify) #library to plot PCA

#run PCA
pcDat <- prcomp(t(rlogcounts))

#plotting PCA
library(ggrepel)
autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
  geom_text_repel(aes(x = PC1, y = PC2, label = SampleName), 
                  box.padding = 0.8
                  )

#correct incorrectly labelled sample info
library(dplyr)
sampleinfo <- mutate(sampleinfo, Status = case_when(
    SampleName == "SRR7657873" ~ "Infected",
    SampleName == "SRR7657882" ~ "Uninfected",
    TRUE ~ Status))

write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")











