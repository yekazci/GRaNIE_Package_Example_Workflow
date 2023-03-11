
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

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

### to be continued....


