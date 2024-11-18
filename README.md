
# Alternative splicing analysis of Arabidopsis under Arsenic stress

```

library(GenomicFeatures)  # for the exonicParts() function
library(Rsamtools)
library(txdbmaker)
library(GenomicAlignments)
library(GenomicRanges)
library(DEXSeq)
BiocManager::install("DEXSeq", force=TRUE)


download.file(
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.60.gtf.gz",
  destfile="Arabidopsis_thaliana.TAIR10.60.gtf.gz")
txdb = makeTxDbFromGFF("Arabidopsis_thaliana.TAIR10.60.gtf.gz")
file.remove("Arabidopsis_thaliana.TAIR10.60.gtf.gz")

flattenedAnnotation <- exonicParts(txdb, linked.to.single.gene.only = TRUE)
flattenedAnnotation


names(flattenedAnnotation) <- sprintf("%s:E%0.3d", 
                                      flattenedAnnotation$gene_id,
                                      flattenedAnnotation$exonic_part)
names(flattenedAnnotation)

bamFiles_1 <- c("AC1.bam",
             "AC2.bam",
             "AC3.bam",
             "CC1.bam",
             "CC2.bam",
             "CC3.bam")

bamFiles_2 <- c("AOE1.bam",
                "AOE2.bam",
                "AOE3.bam",
                "COE1.bam",
                "COE2.bam",
                "COE3.bam")

bamFiles_3 <- c(
  "AFB1.bam",
  "AFB2.bam",
  "AFB3.bam",
  "CFB1.bam",
  "CFB2.bam",
  "CFB3.bam")

bamFiles1 <- BamFileList(bamFiles_1)

bamFiles2 <- BamFileList(bamFiles_2)

bamFiles3 <- BamFileList(bamFiles_3)



bamFiles

if (any(!file.exists(path(bamFiles1)))) {
  stop("One or more BAM files do not exist.")
}

seqlevelsStyle(flattenedAnnotation) <- "NCBI"

seqlevelsStyle(flattenedAnnotation)

se1 <- summarizeOverlaps(
  features = flattenedAnnotation,
  reads = bamFiles1,
  mode = "Union", 
  ignore.strand = TRUE,
)

se2 <- summarizeOverlaps(
  features = flattenedAnnotation,
  reads = bamFiles2,
  mode = "Union", 
  singleEnd = FALSE,
  ignore.strand = TRUE,
  fragments = TRUE
)


se3 <- summarizeOverlaps(
  features = flattenedAnnotation,
  reads = bamFiles3,
  mode = "Union", 
  singleEnd = FALSE,
  ignore.strand = TRUE,
  fragments = TRUE
  
)

combined_se <- cbind(se1, se2, se3)

combined_se

BiocManager::install("DEXSeq", force=TRUE)

library(DEXSeq)


combined_se@assays@data@listData

colData(combined_se)$condition = factor(c("Arsenic",
                                 "Arsenic",
                                 "Arsenic",
                                 "Control",
                                 "Control",
                                 "Control",
                                 "Arsenic",
                                 "Arsenic",
                                 "Arsenic",
                                 "Control",
                                 "Control",
                                 "Control",
                                 "Arsenic",
                                 "Arsenic",
                                 "Arsenic",
                                 "Control",
                                 "Control",
                                 "Control"
                                 ))

colData(combined_se)$libType = factor(c("paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end",
                               "paired-end"))

dxd <- DEXSeqDataSetFromSE(combined_se, design = ~sample + exon + condition:exon)
  
dxd

colData(dxd)

head(counts(dxd), 5)

split( seq_len(ncol(dxd)), colData(dxd)$exon )

head( featureCounts(dxd), 5 )

head( rowRanges(dxd), 3 )

dxd = estimateSizeFactors( dxd )

dxd = estimateDispersions( dxd )

plotDispEsts(dxd)

dxd = testForDEU(dxd)

dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", denominator = "Control")

dxr1 = DEXSeqResults( dxd )

dxr1

dfexon <- as.data.frame(dxr1)

filtered<- subset(dfexon,padj < 0.05 & abs(log2fold_Arsenic_Control) > 1)

filtered <- filtered[!is.na(filtered$padj), ]


write_xlsx(filtered, path = "DEXSeq_resultsAs0.05.csv")



###Plotting#


plotDEXSeq( dxr1, "AT4G04840", cex.axis=1, cex=1, lwd=3, FDR=0.01, names = TRUE , 
            splicing=TRUE, legend = TRUE, expression = TRUE)

plotDEXSeq( dxr1, "AT2G38330", cex.axis=1, cex=1, lwd=3, FDR=0.01, names = TRUE , 
              splicing=TRUE, legend = TRUE, expression = TRUE)


plotDEXSeq( dxr1, "AT1G15520", cex.axis=1, cex=1, lwd=3, FDR=0.01, names = TRUE , 
            splicing=TRUE, legend = TRUE, expression = TRUE)


plotDEXSeq( dxr1, "AT5G37500", cex.axis=1, cex=1, lwd=3, FDR=0.01, names = TRUE , 
            splicing=TRUE, legend = TRUE, expression = TRUE)

```
