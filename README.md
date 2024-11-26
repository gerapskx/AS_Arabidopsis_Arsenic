
# Alternative splicing analysis in _Arabidopsis thaliana_

```

library(GenomicFeatures)
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

filtered <- filtered[order(filtered$log2fold_Arsenic_Control, decreasing = TRUE), ]

filtered$DESCRIPTION <- mapIds(org.At.tair.db,
                                              keys = filtered$groupID,
                                              column = "GENENAME",
                                              keytype = "TAIR",
                                              multiVals = "first")

filtered$gene_name <- mapIds(org.At.tair.db,
                                            keys = filtered$groupID,
                                            column = "SYMBOL",
                                            keytype = "TAIR",
                                            multiVals = "first")

write_xlsx(filtered, path = "DEXSeq_resultsAs0.05_AsvsControl.xlsx")



###Plotting#


plotDEXSeq( dxr1, "AT4G04840", cex.axis=1, cex=1, lwd=3, FDR=0.01, names = TRUE , 
            splicing=TRUE, legend = TRUE, expression = TRUE)

plotDEXSeq( dxr1, "AT2G38330", cex.axis=1, cex=1, lwd=3, FDR=0.01, names = TRUE , 
            splicing=TRUE, legend = TRUE, expression = TRUE)


plotDEXSeq( dxr1, "AT1G15520", cex.axis=1, cex=1, lwd=3, FDR=0.01, names = TRUE , 
            splicing=TRUE, legend = TRUE, expression = TRUE)


plotDEXSeq(dx, "AT2G16220", cex.axis=1, cex=1, lwd=3, FDR=0.01, names = TRUE , 
            splicing=TRUE, legend = TRUE, expression = TRUE)

###Morevariables

colData(combined_se)$genotype = factor(c("WT",
                                          "WT",
                                          "WT",
                                          "WT",
                                          "WT",
                                          "WT",
                                          "OE",
                                          "OE",
                                          "OE",
                                          "OE",
                                          "OE",
                                          "OE",
                                          "KO",
                                          "KO",
                                          "KO",
                                          "KO",
                                          "KO",
                                          "KO"
))

design = ~ genotype + condition + genotype:condition + exon
remove(dxd_COMPLEX)
dxd_COMPLEX <- DEXSeqDataSetFromSE(combined_se, design = ~sample + exon + genotype:condition )
sample + exon + condition:exon
dxd_COMPLEX

colData(dxd_COMPLEX)

head(counts(dxd_COMPLEX), 5)





split( seq_len(ncol(dxd_COMPLEX)), colData(dxd_COMPLEX)$exon )

head( featureCounts(dxd_COMPLEX), 5 )

head( rowRanges(dxd_COMPLEX), 3 )

dxd_COMPLEX = estimateSizeFactors(dxd_COMPLEX)

dxd_COMPLEX = estimateDispersions( dxd_COMPLEX )

plotDispEsts(dxd_COMPLEX)


dxd_COMPLEX = testForDEU(dxd_COMPLEX)

dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", denominator = "Control")

dxr1 = DEXSeqResults( dxd )

dxr1

dfexon <- as.data.frame(dxr1)

filtered<- subset(dfexon,padj < 0.05 & abs(log2fold_Arsenic_Control) > 1)

filtered <- filtered[!is.na(filtered$padj), ]

filtered <- filtered[order(filtered$log2fold_Arsenic_Control, decreasing = TRUE), ]

filtered$DESCRIPTION <- mapIds(org.At.tair.db,
                               keys = filtered$groupID,
                               column = "GENENAME",
                               keytype = "TAIR",
                               multiVals = "first")

filtered$gene_name <- mapIds(org.At.tair.db,
                             keys = filtered$groupID,
                             column = "SYMBOL",
                             keytype = "TAIR",
                             multiVals = "first")

write_xlsx(filtered, path = "DEXSeq_resultsAs0.05_AsvsControl.xlsx")





####IndividualF_box_-___###

bamFiles_3 <- c(
  "AFB1.bam",
  "AFB2.bam",
  "AFB3.bam",
  "CFB1.bam",
  "CFB2.bam",
  "CFB3.bam")

se3@assays@data@listData

colData(se3)$condition = factor(c("Arsenic",
                                          "Arsenic",
                                          "Arsenic",
                                          "Control",
                                          "Control",
                                          "Control"
))

colData(se3)$libType = factor(c("paired-end",
                                        "paired-end",
                                        "paired-end",
                                        "paired-end",
                                        "paired-end",
                                        "paired-end"))

dxd_fbox <- DEXSeqDataSetFromSE(se3, design = ~sample + exon + condition:exon)

dxd_fbox

colData(dxd)

head(counts(dxd_fbox), 5)

split( seq_len(ncol(dxd_fbox)), colData(dxd_fbox)$exon )

head( featureCounts(dxd_fbox), 5 )

head( rowRanges(dxd_fbox), 3 )

dxd_fbox = estimateSizeFactors( dxd_fbox )

dxd_fbox = estimateDispersions( dxd_fbox )

plotDispEsts(dxd_fbox)

dxd_fbox = testForDEU(dxd_fbox)

dxd_fbox = estimateExonFoldChanges( dxd_fbox, fitExpToVar="condition", denominator = "Control")

dxd_fbox1 = DEXSeqResults( dxd_fbox )

BiocManager::install("org.At.tair.db")
library(org.At.tair.db)

dxd_fbox1

dxd_fbox1_df <- as.data.frame(dxd_fbox1)

dxd_fbox1_df_filtered <- subset(dxd_fbox1_df, padj < 0.05 & abs(log2fold_Arsenic_Control) > 1)

dxd_fbox1_df_filtered <- dxd_fbox1_df_filtered[!is.na(dxd_fbox1_df_filtered$padj), ]

dxd_fbox1_df_filtered <- dxd_fbox1_df_filtered[order(dxd_fbox1_df_filtered$log2fold_Arsenic_Control, decreasing = TRUE), ]

dxd_fbox1_df_filtered$DESCRIPTION <- mapIds(org.At.tair.db,
                                 keys = dxd_fbox1_df_filtered$groupID,
                                 column = "GENENAME",
                                 keytype = "TAIR",
                                 multiVals = "first")

write_xlsx(dxd_fbox1_df_filtered, path = "DEXSeq_resultsAs0.05_Fbox.xlsx")


###OE############


se2@assays@data@listData

colData(se2)$condition = factor(c("Arsenic",
                                  "Arsenic",
                                  "Arsenic",
                                  "Control",
                                  "Control",
                                  "Control"
))

colData(se2)$libType = factor(c("paired-end",
                                "paired-end",
                                "paired-end",
                                "paired-end",
                                "paired-end",
                                "paired-end"))

dxd_fboxOE <- DEXSeqDataSetFromSE(se2, design = ~sample + exon + condition:exon)

dxd_fboxOE

colData(dxd_fboxOE)

head(counts(dxd_fboxOE), 5)

split( seq_len(ncol(dxd_fboxOE)), colData(dxd_fboxOE)$exon )

head( featureCounts(dxd_fboxOE), 5 )

head( rowRanges(dxd_fboxOE), 3 )

dxd_fboxOE = estimateSizeFactors( dxd_fboxOE )

dxd_fboxOE = estimateDispersions( dxd_fboxOE )

plotDispEsts(dxd_fboxOE)


dxd_fboxOE = testForDEU(dxd_fboxOE)

dxd_fboxOE = estimateExonFoldChanges( dxd_fboxOE, fitExpToVar="condition", denominator = "Control")

dxd_fboxOE1 = DEXSeqResults(dxd_fboxOE)

dxd_fboxOE1

dxd_fboxOE1_df <- as.data.frame(dxd_fboxOE1)


dxd_fboxOE1_df_filtered <- subset(dxd_fboxOE1_df, padj < 0.05 & abs(log2fold_Arsenic_Control) > 1)

dxd_fboxOE1_df_filtered <- dxd_fboxOE1_df_filtered[!is.na(dxd_fboxOE1_df_filtered$padj), ]

dxd_fboxOE1_df_filtered <- dxd_fboxOE1_df_filtered[order(dxd_fboxOE1_df_filtered$log2fold_Arsenic_Control, decreasing = TRUE), ]

dxd_fboxOE1_df_filtered$DESCRIPTION <- mapIds(org.At.tair.db,
                                            keys = dxd_fboxOE1_df_filtered$groupID,
                                            column = "GENENAME",
                                            keytype = "TAIR",
                                            multiVals = "first")

dxd_fboxOE1_df_filtered$gene_name <- mapIds(org.At.tair.db,
                                              keys = dxd_fboxOE1_df_filtered$groupID,
                                              column = "SYMBOL",
                                              keytype = "TAIR",
                                              multiVals = "first")

write_xlsx(dxd_fboxOE1_df_filtered, path = "DEXSeq_resultsAs0.05_Fbox_OE.xlsx")


###Colombia############


se1@assays@data@listData

colData(se1)$condition = factor(c("Arsenic",
                                  "Arsenic",
                                  "Arsenic",
                                  "Control",
                                  "Control",
                                  "Control"
))

colData(se1)$libType = factor(c("paired-end",
                                "paired-end",
                                "paired-end",
                                "paired-end",
                                "paired-end",
                                "paired-end"))

dxd_Col <- DEXSeqDataSetFromSE(se1, design = ~sample + exon + condition:exon)

dxd_Col

colData(dxd_Col)

head(counts(dxd_Col), 5)

split( seq_len(ncol(dxd_Col)), colData(dxd_Col)$exon )

head( featureCounts(dxd_Col), 5 )

head( rowRanges(dxd_Col), 3 )

dxd_Col = estimateSizeFactors( dxd_Col )

dxd_Col = estimateDispersions( dxd_Col )

plotDispEsts(dxd_Col)


dxd_Col = testForDEU(dxd_Col)

dxd_Col = estimateExonFoldChanges( dxd_Col, fitExpToVar="condition", denominator = "Control")

dxd_Col1 = DEXSeqResults(dxd_Col)

dxd_Col1

dxd_Col1_df <- as.data.frame(dxd_Col1)


dxd_Col1_df_filtered <- subset(dxd_Col1_df, padj < 0.05 & abs(log2fold_Arsenic_Control) > 1)

dxd_Col1_df_filtered <- dxd_Col1_df_filtered[!is.na(dxd_Col1_df_filtered$padj), ]

dxd_Col1_df_filtered <- dxd_Col1_df_filtered[order(dxd_Col1_df_filtered$log2fold_Arsenic_Control, decreasing = TRUE), ]

dxd_Col1_df_filtered$DESCRIPTION <- mapIds(org.At.tair.db,
                                              keys = dxd_Col1_df_filtered$groupID,
                                              column = "GENENAME",
                                              keytype = "TAIR",
                                              multiVals = "first")

dxd_Col1_df_filtered$gene_name <- mapIds(org.At.tair.db,
                                            keys = dxd_Col1_df_filtered$groupID,
                                            column = "SYMBOL",
                                            keytype = "TAIR",
                                            multiVals = "first")

write_xlsx(dxd_Col1_df_filtered, path = "DEXSeq_resultsAs0.05_Col.xlsx")


########Plotting


#arsenic
plotDEXSeq(dxr1, "AT1G15520", 
           cex.axis=1, cex=1, lwd=3, FDR=0.1, 
           legend = TRUE, expression = TRUE, splicing = TRUE)
#colombia
plotDEXSeq(dxd_Col1, "AT1G15520", 
           cex.axis=1, cex=1, lwd=3, FDR=0.1, 
           legend = TRUE,names = TRUE, expression = TRUE, splicing = TRUE)
#KO
plotDEXSeq(dxd_fbox1, "AT1G15520",
           cex.axis=1, cex=1, lwd=3, FDR=0.1, 
           legend = TRUE, names = TRUE, expression = TRUE, splicing = TRUE)
#OE
plotDEXSeq(dxd_fboxOE1, "AT1G15520", 
           cex.axis=1, cex=1, lwd=3, FDR=0.1,
           legend = TRUE, names = TRUE, expression = TRUE, splicing = TRUE)



```

