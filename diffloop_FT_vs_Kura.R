## loading packages
library(diffloop)
library(diffloopdata)
library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(DESeq2)
library(dplyr)

## loading FitHiChip Q-value bed files
data1 <- read.delim("FitHiChIP.interactions_Kura-A_FitHiC_Q0.01.bed")
data2 <- read.delim("FitHiChIP.interactions_Kura-B_FitHiC_Q0.01.bed")
data3 <- read.delim("FitHiChIP.interactions_UWB-A_FitHiC_Q0.01.bed")
data4 <- read.delim("FitHiChIP.interactions_UWB-B_FitHiC_Q0.01.bed")
data5 <- read.delim("FitHiChIP.interactions_FT246-A_FitHiC_Q0.01.bed")
data6 <- read.delim("FitHiChIP.interactions_FT246-B_FitHiC_Q0.01.bed")
data7 <- read.delim("FitHiChIP.interactions_FT33-A_FitHiC_Q0.01.bed")
data8 <- read.delim("FitHiChIP.interactions_FT33-B_FitHiC_Q0.01.bed")

## subsetting the dataframe to required columns
data1.1 <- data1[,c(1:7)]
data2.2 <- data2[,c(1:7)]
data3.3 <- data3[,c(1:7)]
data4.4 <- data4[,c(1:7)]
data5.5 <- data5[,c(1:7)]
data6.6 <- data6[,c(1:7)]
data7.7 <- data7[,c(1:7)]
data8.8 <- data8[,c(1:7)]

## Removing the chr character infront of chromosome number
data1.1$chr1 <- gsub("chr","",as.character(data1.1$chr1))
data1.1$chr2 <- gsub("chr","",as.character(data1.1$chr2))
data2.2$chr1 <- gsub("chr","",as.character(data2.2$chr1))
data2.2$chr2 <- gsub("chr","",as.character(data2.2$chr2))
data3.3$chr1 <- gsub("chr","",as.character(data3.3$chr1))
data3.3$chr2 <- gsub("chr","",as.character(data3.3$chr2))
data4.4$chr1 <- gsub("chr","",as.character(data4.4$chr1))
data4.4$chr2 <- gsub("chr","",as.character(data4.4$chr2))
data5.5$chr1 <- gsub("chr","",as.character(data5.5$chr1))
data5.5$chr2 <- gsub("chr","",as.character(data5.5$chr2))
data6.6$chr1 <- gsub("chr","",as.character(data6.6$chr1))
data6.6$chr2 <- gsub("chr","",as.character(data6.6$chr2))
data7.7$chr1 <- gsub("chr","",as.character(data7.7$chr1))
data7.7$chr2 <- gsub("chr","",as.character(data7.7$chr2))
data8.8$chr1 <- gsub("chr","",as.character(data8.8$chr1))
data8.8$chr2 <- gsub("chr","",as.character(data8.8$chr2))

## Inserting dot as column in dataframe
data1.1[,'dot'] <- "."
data2.2[,'dot'] <- "."
data3.3[,'dot'] <- "."
data4.4[,'dot'] <- "."
data5.5[,'dot'] <- "."
data6.6[,'dot'] <- "."
data7.7[,'dot'] <- "."
data8.8[,'dot'] <- "."

## Rearranging the dataframe columns
kura_1 <- data1.1[,c(1:6,8,7)]
kura_2 <- data2.2[,c(1:6,8,7)]
uwb_1 <- data3.3[,c(1:6,8,7)]
uwb_2 <- data4.4[,c(1:6,8,7)]
ft246_1 <- data5.5[,c(1:6,8,7)]
ft246_2 <- data6.6[,c(1:6,8,7)]
ft33_1 <- data7.7[,c(1:6,8,7)]
ft33_2 <- data8.8[,c(1:6,8,7)]

## Writing bedpe files
write.table(kura_1, file = "kura_1.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(kura_2, file = "kura_2.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(uwb_1, file = "uwb_1.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(uwb_2, file = "uwb_2.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(ft246_1, file = "ft246_1.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(ft246_2, file = "ft246_2.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(ft33_1, file = "ft33_1.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)
write.table(ft33_2, file = "ft33_2.loop_counts.bedpe", row.names = FALSE,
            col.names = FALSE)

## diffloop downstream analysis with bedpe files
bed_dir <- file.path("/Users/ayalurik/Documents/HiChIP_diffloop/diffloop")
bed_dir

samples <- c("kura_1","kura_2","uwb_1","uwb_2",
             "ft246_1","ft246_2","ft33_1","ft33_2")
full <- loopsMake(bed_dir, samples)
celltypes <- c("kura", "kura", "uwb", "uwb", "ft246", "ft246", "ft33", "ft33")
full <- updateLDGroups(full, celltypes)
head(full, 4)
dim(full)

## Quality control
# remove and loops that merged together from import
full1 <- subsetLoops(full, full@rowData$loopWidth >= 5000)
dim(full1)

# Filtering loops that are FDR > 0.01
noCNV_corrected <- mangoCorrection(full1, FDR = 0.01)
dim(noCNV_corrected)

# filtering loops that are present strongly (>= 5 PETs) in one replicate
# but absent (== 0 PETs) in the other replicate
cm <- noCNV_corrected@counts
k_dis <- ((cm[,1]>=5 &cm[,2]==0)|(cm[,2]>=5&cm[,1]==0))
m_dis <- ((cm[,3]>=5 &cm[,4]==0)|(cm[,4]>=5&cm[,3]==0))
qc_filt <- subsetLoops(noCNV_corrected, !(k_dis | m_dis))
dim(qc_filt)

p1 <- loopDistancePlot(qc_filt)
p1

## counts the number of loops for each sample and returns whether they are single,
## self, unique or none
loopMetrics(qc_filt)

# PC plot
pcp1dat <- qc_filt
pcp1dat@colData$sizeFactor <- 1
samples <- c("kura_1", "kura_2", "uwb_1", "uwb_2",
             "ft246_1", "ft246_2", "ft33_1", "ft33_2")
pcp1 <- pcaPlot(pcp1dat) + geom_text_repel(aes(label=samples)) + 
  ggtitle("PC Plot with no
  Size Factor Correction") + theme(legend.position="none")
pcp1

samples <- c("kura_1", "kura_2", "uwb_1", "uwb_2",
             "ft246_1", "ft246_2", "ft33_1", "ft33_2")
pcp2 <- pcaPlot(qc_filt) + geom_text_repel(aes(label=samples)) + 
  ggtitle("PC Plot with Size Factor Correction") +
  theme(legend.position = "none")
pcp2


## Differential Loop Calling
km_filt <- qc_filt
dim(km_filt)


# First model using edgeR over-dispersed Poisson regression
km_res <- loopAssoc(km_filt, method = "edgeR", coef = 2)
head(km_res@rowData)


# # second model using limma-voom empirical Bayes analysis
# km_res2 <- quickAssocVoom(km_res)
# 
# fdrmat <- -log10(km_res2@rowData[,c(8,13)])
# fdrmatdf <- setNames(data.frame(fdrmat), c("edgeR", "Voom"))
# qplot(fdrmatdf$Voom, fdrmatdf$edgeR)+labs(title = "diffloop: edgeR versus Voom", 
#                                           x = "-log10(FDR) Voom", y = "-log10(FDR) edgeR") + scale_y_continuous(limits = c(0, 25)) +        
#   scale_x_continuous(limits = c(0, 25)) + theme_bw()
# 
# cor(fdrmatdf$Voom, fdrmatdf$edgeR)
# cor(fdrmatdf$Voom, fdrmatdf$edgeR, method = "spearman")

## Epigenetic Annotation
h3dir <- file.path("/Users/ayalurik/Documents/HiChIP_diffloop")
kura <- paste0(h3dir, "/", "kura_H3k27ac_hg19_narrowpeak.bed")
uwb <- paste0(h3dir, "/", "UWB_H3k27ac_hg19_narrowpeak.bed")
ft246 <- paste0(h3dir, "/", "FT33_H3k27ac_hg19_narrowpeak.bed")
ft33 <- paste0(h3dir, "/", "FT246_H3k27ac_hg19_narrowpeak.bed")

h3k27ac.k <- rmchr(padGRanges(bedToGRanges(kura), pad = 1000))
h3k27ac.u <- rmchr(padGRanges(bedToGRanges(uwb), pad = 1000))
h3k27ac.f2 <- rmchr(padGRanges(bedToGRanges(ft246), pad = 1000))
h3k27ac.f3 <- rmchr(padGRanges(bedToGRanges(ft33), pad = 1000))
ku <- GenomicRanges::union(h3k27ac.k, h3k27ac.u)
ft <- GenomicRanges::union(h3k27ac.f2, h3k27ac.f3)
enhancer <- GenomicRanges::union(ku,ft)
promoter <- padGRanges(getHumanTSS(), pad = 1000)
km_anno <- annotateLoops(km_res, enhancer = enhancer, promoter = promoter)

## Visulaization
plotTopLoops(km_anno, organism = "h", FDR = 0.01, n = 10)


