
### Here, I will follow the example workflow from the Granie package vignette.

library(readr)
library(GRaNIE)

# the url links for the example datasets in the vignette :
file_peaks <- "https://www.embl.de/download/zaugg/GRaNIE/countsATAC.filtered.tsv.gz"
file_RNA <- "https://www.embl.de/download/zaugg/GRaNIE/countsRNA.filtered.tsv.gz"
file_sampleMetadata <- "https://www.embl.de/download/zaugg/GRaNIE/sampleMetadata.tsv.gz"

countsRNA.df <- read_tsv(file_RNA, col_types = cols())
countsPeaks.df <- read_tsv(file_peaks, col_types = cols()) 
sampleMetadata.df <- read_tsv(file_sampleMetadata, col_types = cols())


# Viewing the datasets:

countsRNA.df
countsPeaks.df
sampleMetadata.df

# we save the names of the columns including the respective ids as variables.
# we will later use them when creating a GRN object.

idColumn_peaks <- "peakID"
idColumn_RNA <- "ENSEMBL"

# we set the name of the version of human genome to be used as a variable:

genomeAssembly <- "hg38" 

# Important note from the vignette:

# Either hg19, hg38 or mm10. Both enhancers and RNA data must have the same 
# genome assembly.

# An optional list which contains the metadata to be stored in the GRN object:

objectMetadata.l <- list(name = paste0("Macrophages_infected_primed"),
                        file_peaks = file_peaks,
                        file_rna = file_RNA,
                        file_sampleMetadata = file_sampleMetadata,
                        genomeAssembly = genomeAssembly)

## we will save the output to the current working directory.

dir_output = "."    

## Let's initialize our GRN object:

GRN  <- initializeGRN(objectMetadata = objectMetadata.l,
                    outputFolder = dir_output,
                    genomeAssembly = genomeAssembly)

## Checking if the GRN object was properly formed:

GRN

# Adding the example datasets to the GRN object:

GRN <- addData(GRN,
               counts_peaks = countsPeaks.df, 
               normalization_peaks = "DESeq_sizeFactor", 
               idColumn_peaks = idColumn_peaks,
               counts_rna = countsRNA.df, 
               normalization_rna = "quantile", 
               idColumn_RNA = idColumn_RNA,
               sampleMetadata = sampleMetadata.df,
               forceRerun = TRUE)

# Note: we can view the parameters of the function by "?addData".

# Note that we used different normalization methods for the RNA and peaks 
# data.

### For brevity, we quality check only a type of RNA data. There are raw and 
### normalized types of each RNA and Chip (Enhancer) data. Here as an example, 
### we check the separation of samples in PC space by using 500 mostly varied
### genes.

GRN <- plotPCA_all(GRN, data = c("rna"), topn = 500, type = "normalized", 
plotAsPDF = FALSE, pages = c(2,3,14), forceRerun = TRUE)

### I saved the output images to the plot folder. From Scree plot, 
### we see that 10 PC components explain the 80 % of variation within the 29
### samples. Some metadata features lowly correlate with some PC components 
### (max coefficient : 0.31), which seems fair. In PC1 and PC2 component space,
### samples are not distributed based on their mitochondrial gene fractions. 

### Next, transcription factor binding site data will be downloaded
### It contains data for 6 TFs.

folder_TFBS_6TFs <- "https://www.embl.de/download/zaugg/GRaNIE/TFBS_selected.zip"

# I download the zip of all TFBS files. It might take somke time. It can be also
# downloaded through web browser.

download.file(folder_TFBS_6TFs, file.path("TFBS_selected.zip"), quiet = FALSE)

### Optional: saveRDS(object = GRN, file = "GRN.RDS")

unzip(file.path("TFBS_selected.zip"), overwrite = TRUE)

### View the files in the directory:

list.files(path = "TFBS_selected")

# [1] "BATF3.0.B_TFBS.bed.gz" "E2F6.0.A_TFBS.bed.gz"  "E2F7.0.B_TFBS.bed.gz" 
# [4] "EGR1.0.A_TFBS.bed.gz"  "EGR2.0.A_TFBS.bed.gz"  "ETV5.0.C_TFBS.bed.gz" 
# [7] "translationTable.csv" 

### Save the full path into a variable:
 
motifFolder <- tools::file_path_as_absolute("TFBS_selected")

### Now, I add the TFBS info to my GRN object:
  
GRN <- addTFBS(GRN, motifFolder = motifFolder, 
               TFs = "all", 
               filesTFBSPattern = "_TFBS", 
               fileEnding = ".bed.gz", 
               forceRerun = TRUE)
  
GRN <- overlapPeaksAndTFBS(GRN, nCores = 1, forceRerun = TRUE)

# We specify the chromosomes to be kept for the peaks. 
# For that, we form a new vector including the  chromosome names:

chrToKeep_peaks <- c(paste0("chr", 1:22), "chrX")

# We will filter out the RNA and peak data and also specify the chromosomes
# for which the peaks will be kept.

GRN <- filterData(GRN, minNormalizedMean_peaks = 5, 
                  minNormalizedMeanRNA = 1, 
                  chrToKeep_peaks = chrToKeep_peaks, 
                  maxSize_peaks = 10000, 
                  forceRerun = TRUE)

# 

GRN <- addConnections_TF_peak(GRN, plotDiagnosticPlots = FALSE,
                             connectionTypes = c("expression"),
                             corMethod = "pearson", forceRerun = TRUE)

# TF-enhancer diagnostic plots for EGR1.0.A :

GRN <- plotDiagnosticPlots_TFPeaks(GRN, dataType = c("real", "permuted"), 
                                   plotAsPDF = FALSE, pages = 4)

GRN <- AR_classification_wrapper(GRN, significanceThreshold_Wilcoxon = 0.05,
                                 outputFolder = "plots",
                                 plot_minNoTFBS_heatmap = 100, 
                                 plotDiagnosticPlots = TRUE,
                                 forceRerun = TRUE)


GRN_file_outputRDS <- paste0(dir_output, "/GRN.RDS")
 
saveRDS(GRN, GRN_file_outputRDS)

## add Peak Gene Connections

GRN <- addConnections_peak_gene(GRN,
                               corMethod = "pearson", 
                               promoterRange = 10000, 
                               TADs = NULL,
                               nCores = 1, 
                               plotDiagnosticPlots = FALSE, 
                               plotGeneTypes = list(c("all")), 
                               forceRerun = TRUE)

## Enhancer-gene diagnostic plots:

GRN <- 
  plotDiagnosticPlots_peakGene(GRN, 
                          gene.types = list(c("protein_coding", "lincRNA")), 
                          plotAsPDF = FALSE, pages = 1)

## combine And Filter, 

GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = 0.2,
                               peak_gene.fdr.threshold = 0.2, 
                               peak_gene.fdr.method = "BH", 
                               gene.types = c("protein_coding", "lincRNA"),
                               allowMissingTFs = FALSE, 
                               allowMissingGenes = FALSE)


GRN <- add_TF_gene_correlation(GRN, corMethod = "pearson", 
                               nCores = 1, forceRerun = TRUE)


# Optional convenience function to delete intermediate data from the function 
# AR_classification_wrapper and summary statistics that may occupy 
# a lot of space: 

#  GRN <- deleteIntermediateData(GRN)

#  saveRDS(GRN, GRN_file_outputRDS)

#  saveRDS(GRN, GRN_file_outputRDS)

GRN_connections.all <- getGRNConnections(GRN, type = "all.filtered", 
                                         include_TF_gene_correlations = TRUE, 
                                         include_geneMetadata = TRUE)

GRN_connections.all

## connectionSummary, 
GRN <- generateStatsSummary(GRN,
                           TF_peak.fdr = c(0.05, 0.1, 0.2), 
                           TF_peak.connectionTypes = "all",
                           peak_gene.fdr = c(0.1, 0.2),
                           peak_gene.r_range = c(0,1),
                           allowMissingGenes = c(FALSE, TRUE),
                           allowMissingTFs = c(FALSE),
                           gene.types = c("protein_coding", "lincRNA"),
                           forceRerun = TRUE)


GRN <- plot_stats_connectionSummary(GRN, type = "heatmap", 
                                    plotAsPDF = FALSE, pages = 3)

GRN <- plot_stats_connectionSummary(GRN, type = "boxplot", 
                                    plotAsPDF = FALSE, pages = 1)

## build Graph:

GRN <- build_eGRN_graph(GRN, forceRerun = TRUE) 


## ----visualizeGRN, an example visualisation.

GRN <- visualizeGRN(GRN, plotAsPDF = FALSE)


## ----allNetworkAnalyses, 

#  GRN <- performAllNetworkAnalyses(GRN, ontology = c("GO_BP"), 
#  outputFolder = ".", forceRerun = TRUE)

## General network statistics for the filtered network:

GRN <- plotGeneralGraphStats(GRN, plotAsPDF = FALSE, pages = c(1,6))


## generalEnrichment, 
#  
#  GRN <- calculateGeneralEnrichment(GRN, ontology = "GO_BP")
#  

## plotGeneralEnrichment, General network enrichment for the filtered network

GRN <- plotGeneralEnrichment(GRN, plotAsPDF = FALSE, pages = 1)

GRN <- calculateCommunitiesStats(GRN)
GRN <- calculateCommunitiesEnrichment(GRN, ontology = "GO_BP")


## General statistics for the communities from the filtered network:

GRN <- plotCommunitiesStats(GRN, plotAsPDF = FALSE, pages = c(1,3))


##  Community enrichment for 3 different communities,

GRN <- plotCommunitiesEnrichment(GRN, plotAsPDF = FALSE, pages = c(1,2,3))


## Summary of the community enrichment,

GRN<- plotCommunitiesEnrichment(GRN, plotAsPDF = FALSE, pages = c(5))



## TFEnrichmentCal,

GRN <- calculateTFEnrichment(GRN, ontology = "GO_BP")

## Enrichment summary for EGR1.0.A

GRN <- plotTFEnrichment(GRN, plotAsPDF = FALSE, n = 3, pages = c(1))

## TFEnrichment2, Enrichment summary for selected TFs and the whole eGRN network

GRN <- plotTFEnrichment(GRN, plotAsPDF = FALSE, n = 3, pages = c(5))

#  GRN = deleteIntermediateData(GRN)
#  saveRDS(GRN, GRN_file_outputRDS)


sessionInfo()




writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

### to be continued....


