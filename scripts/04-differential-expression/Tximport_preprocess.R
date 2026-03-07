#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("readr")
#BiocManager::install("tximport")
#BiocManager::install("tximportData")
#BiocManager::install("DESeq2")
#BiocManager::install("biomaRt")

library(dplyr)
library(tidyr)
library(readr)
library(tximport)
library(tximportData)
library(DESeq2)
library(biomaRt)

## Importing sample metadata for raw transcript abundance files

# Samples metadata includes names of salmon output files 
salmon_dir <- "/Users/erinvoss/Documents/Peromyscus-RNA/Salmon_Quant_Output"
samples <- read.table(file.path(salmon_dir, "Salmon_Quant_Files_Experimental_Groups.txt"), header = FALSE)
colnames(samples) <- c("sample", "name", "species", "tissue", "replicate")
head(samples)

# Save individual species metadata
samples.PEBO <- samples[1:18,]
samples.PECA <- samples[19:38,]
samples.PEMA <- samples[39:58,]


## Read in counts, lengths, and metadata as created above for all species, all tissues, no filtering
## Specify directory to save deseq results to  

## Import quant.sf salmon files with read counts and lengths (by transcript)
files <- file.path(salmon_dir, samples$sample)
names(files) <- paste0(samples$name)

files.PEBO <- file.path(salmon_dir, samples.PEBO$sample) # create subset of P. boylii quant.sf salmon files
names(files.PEBO) <- paste0(samples.PEBO$name)

files.PECA <- file.path(salmon_dir, samples.PECA$sample) # create subset of P. californicus quant.sf salmon files
names(files.PECA) <- paste0(samples.PECA$name)
                                 
files.PEMA <- file.path(salmon_dir, samples.PEMA$sample) # create subset of P. maniculatus quant.sf salmon files
names(files.PEMA) <- paste0(samples.PEMA$name)

all(file.exists(files)) # check that all files we named exist



## Use tximport to read in tx2gene files 
## Key preprocessing step: this maps individual transcripts to their parent genes for each species transcriptome
dir <- "./DeSeq_Gene_IDs_ENSPEMG"

tx2gene.PEBO <- read.csv(file.path(dir, "PEBO_transcriptome_annot_ENSPEMG_final.tx2gene.tsv"), sep = "\t", header = FALSE)
colnames(tx2gene.PEBO) <- c("Transcript_ID", "Gene_ID")

tx2gene.PECA <- read.csv(file.path(dir, "PECA_transcriptome_annot_ENSPEMG_final.tx2gene.tsv"), sep = "\t", header = FALSE)
colnames(tx2gene.PECA) <- c("Transcript_ID", "Gene_ID")

tx2gene.PEMA <- read.csv(file.path(dir, "PEMA_transcriptome_ENSPEMG_all_isoforms.tx2gene.tsv"), sep = "\t", header = FALSE)
colnames(tx2gene.PEMA) <- c("Transcript_ID", "Gene_ID")

head(tx2gene.PEBO)



## Reading in salmon quant files using tximport and tx2gene map 

txi.PEBO <- tximport(files.PEBO, type = "salmon", tx2gene = tx2gene.PEBO)
txi.PECA <- tximport(files.PECA, type = "salmon", tx2gene = tx2gene.PECA)
txi.PEMA <- tximport(files.PEMA, type = "salmon", tx2gene = tx2gene.PEMA)

## formatting txi counts as dataframes
txi.PEBO.counts <- as.data.frame(txi.PEBO$counts)
txi.PECA.counts <- as.data.frame(txi.PECA$counts)
txi.PEMA.counts <- as.data.frame(txi.PEMA$counts)

## Merging to create one count dataframe for all species for downstream analyses
## Transcript abundance was estimated separately for each species to their specific transcriptome
## but now we need to compare across species 
## Since we have used the same gene naming conventions (ENSPEMG - P. maniculatus ENSEMBL genes) across all three transcriptomes
## We can now combine counts by gene name 
## Don't worry, we will do a bunch of downstream normalizing to make sure that this is accounted for in differential expression analyses

txi.PEBO_PECA.counts <- merge(txi.PEBO.counts, txi.PECA.counts, by = 0)
row.names(txi.PEBO_PECA.counts) <- txi.PEBO_PECA.counts$Row.names
txi.PEBO_PECA.counts <- txi.PEBO_PECA.counts[,2:39]

txi.PECA_PEMA.counts <- merge(txi.PECA.counts, txi.PEMA.counts, by = 0)
row.names(txi.PECA_PEMA.counts) <- txi.PECA_PEMA.counts$Row.names
txi.PECA_PEMA.counts <- txi.PECA_PEMA.counts[,2:41]

txi.allspecies.counts <- merge(txi.PEBO_PECA.counts, txi.PEMA.counts, by = 0)
row.names(txi.allspecies.counts) <- txi.allspecies.counts$Row.names
txi.allspecies.counts <- txi.allspecies.counts[,2:59]

txi.allspecies.counts <- as.matrix(txi.allspecies.counts) # Convert all species counts dataframe to matrix 
txi.allspecies.counts <- round(txi.allspecies.counts) 

write.table(txi.allspecies.counts,file="all_species.counts.ENSPEMG.tsv", sep = "\t") # Write all counts table to .tsv file


## Creating txi length dataframes and merging to create one length DF with all species
## We will use this transcript length data to normalize counts across species (transcriptomes)

txi.PEBO.length <- as.data.frame(txi.PEBO$length)
txi.PECA.length <- as.data.frame(txi.PECA$length)
txi.PEMA.length <- as.data.frame(txi.PEMA$length)

txi.PEBO_PECA.length <- merge(txi.PEBO.length, txi.PECA.length, by = 0)
row.names(txi.PEBO_PECA.length) <- txi.PEBO_PECA.length$Row.names
txi.PEBO_PECA.length <- txi.PEBO_PECA.length[,2:39]

txi.PECA_PEMA.length <-  merge(txi.PECA.length, txi.PEMA.length, by = 0)
row.names(txi.PECA_PEMA.length) <- txi.PECA_PEMA.length$Row.names
txi.PECA_PEMA.length <- txi.PECA_PEMA.length[,2:41]

txi.allspecies.length <- merge(txi.PEBO_PECA.length, txi.PEMA.length, by = 0)
row.names(txi.allspecies.length) <- txi.allspecies.length$Row.names
txi.allspecies.length <- txi.allspecies.length[,2:59]

txi.allspecies.length <- as.matrix(txi.allspecies.length)
txi.allspecies.length <- round(txi.allspecies.length)

write.table(txi.allspecies.length,file="all_species.lengths.ENSPEMG.tsv", sep = "\t") # save for normalization step 




# Creating metadata table for all species samples 

sample_names.allSpecies <- as.factor(colnames(txi.allspecies.counts))

tissue.allSpecies <- c(rep("Testis", 5), rep("SemVes", 4), rep("Epididymis", 4), rep("Liver", 5), rep("Testis", 5), rep("SemVes", 5), rep("Epididymis", 5), rep("Liver", 5), rep("Testis", 5), rep("SemVes", 5), rep("Epididymis", 5), rep("Liver", 5))

species.allSpecies <- c(rep("PEBO", 18), rep("PECA", 20), rep("PEMA", 20))

samplesTable.allSpecies <- data.frame(tissue.allSpecies, species.allSpecies)

row.names(samplesTable.allSpecies) <- sample_names.allSpecies
colnames(samplesTable.allSpecies) <- c("tissue", "species")

write.table(samplesTable.allSpecies,file="all_species.samplesTable.tsv",row.names = sample_names.allSpecies, sep="\t")


# Salmon Quant RNAseq data has now been imported into R and is ready to analyze with DeSeq2 (using RStudio)
