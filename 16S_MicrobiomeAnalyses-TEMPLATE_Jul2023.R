###### 16S Microbiome Analyses - Molecular Microbiology & Genomics Lab
### This R-Script includes routine microbiome analyses conducted using
### the 16S rRNA gene (V4 region (primers 515F/806R)).
### Initial data cleaning starts from QIIME2, so files needed are from QIIME2
### resulting files.
### You may change anything in this script to fit your project analyses.
### Please make sure to keep this original copy and to save any edited versions
### with your project data!
### Created by Paisley Samuel - Summer 2023

###### Set Working Directory and Seed ######
setwd() #include the file path to the folder where all your files are stored
set.seed() #choose any random number

###### Producing Relative Abundance Data ######
### NOTE: If you ignore metdata, you will NOT be able to perform statistical
###       analyses on your data, only taxonomical analyses.

## RUN FUNCTIONS FIRST!
# use "calculate_abund.metadta" if you have your metadata file
# use "calculate_abund" if ignoring metadata (NO metadata file)
calculate_abund.metadata <- function(feature_csv, metadata_csv) {
  library(vegan)
  dat <- t(data.matrix(read.csv(feature_csv, header = TRUE, row.names = 1)))
  #transposed so that the row names are now the sample names and the ASVs are the columns
  #BUT it is still a MATRIX
  metadata <- read.csv(metadata_csv, header = TRUE, row.names = 1)
  dat <- as.data.frame(dat)
  typeof(dat) #list
  common.rownames <- intersect(rownames(dat), rownames(metadata)) #returns the same row names found in each dataframe
  dat <- dat[common.rownames,] #subsets dat to only include the same row names as metadata
  metadata <- metadata[common.rownames,] #subsets metadata to only include the same row names as dat
  otu.abund <- which(colSums(dat) > 2) #removes singletons and doubletons
  dat.dom <- dat[, otu.abund] #include dominant taxa
  dat.pa <- decostand(dat.dom, method ="pa") #turns dat.dom into presence/absence data (1/0)
  dat.otus.01per <- which(colSums(dat.pa) > (0.01 * nrow(dat.pa)))
  dat.01per <- dat.dom[, dat.otus.01per] #removed ASVs that occur less than 0.1%
  dat.otus.001per <- which(colSums(dat.pa) > (0.001 * nrow(dat.pa)))
  dat.001per <- dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; increases the number of ASVs - includes more "microdiversity" 
  dat.ra <- decostand(dat.01per, method = "total") #relative abundance of >1% taxa
  dfs_to_return <- list(as.data.frame(dat), as.data.frame(metadata),
                        as.data.frame(dat.dom), as.data.frame(dat.pa),
                        as.data.frame(dat.01per), as.data.frame(dat.001per),
                        as.data.frame(dat.ra))
  names(dfs_to_return) <- c("dat", "metadata","dat.dom","dat.pa","dat.01per","dat.001per","dat.ra")
  return(dfs_to_return)
}
calculate_abund <- function(feature_csv) {
  library(vegan)
  dat <- t(data.matrix(read.csv(feature_csv, header = TRUE, row.names = 1)))
  #transposed so that the row names are now the sample names and the ASVs are the columns
  dat <- as.data.frame(dat)
  otu.abund <- which(colSums(dat) > 2) #removes singletons and doubletons
  dat.dom <- dat[, otu.abund] #include dominant taxa
  dat.pa <- decostand(dat.dom, method ="pa") #turns dat.dom into presence/absence data (1/0)
  dat.otus.01per <- which(colSums(dat.pa) > (0.01 * nrow(dat.pa)))
  dat.01per <- dat.dom[, dat.otus.01per] #removed ASVs that occur less than 0.1%
  dat.otus.001per <- which(colSums(dat.pa) > (0.001 * nrow(dat.pa)))
  dat.001per <- dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; increases the number of ASVs - includes more "microdiversity" 
  dat.ra <- decostand(dat.01per, method = "total") #relative abundance of >1% taxa
  dfs_to_return <- list(as.data.frame(dat), as.data.frame(dat.dom), as.data.frame(dat.pa),
                        as.data.frame(dat.01per), as.data.frame(dat.001per),
                        as.data.frame(dat.ra))
  names(dfs_to_return) <- c("dat","dat.dom","dat.pa","dat.01per","dat.001per","dat.ra")
  return(dfs_to_return)
}

## insert your files into the function
# "feature.csv" = insert the name of your feature table .csv file
# "metadata.csv" = insert the name of your metadata .csv file

abund_metadata <- calculate_abund.metadata("feature.csv","metadata.csv") #if you have metadata
abund <- calculate_abund("feature.csv") #if you have NO metadata

# functions will return a list of data frames:
# dat = your original feature table that is transposed (or flipped)
# metadata = your metadata table
# **dat and metadata should have the same samples!**
# dat.dom = the dominant ASVs found in dat
# dat.pa = the dat.dom table transformed into presence/absence data
# dat.01per = ASVs from dat that have an abundance of 0.1% or higher
# dat.001per = ASVs from dat that have an abundance of 0.01% or higher
# dat.ra = ASVs from dat.01per normalized into relative abundance

# To save any dataframe (df) within the function result list as a single df in
# your workspace, simply assign it to its own variable
# Example saving dat and metadata
feature <- abund_metadata$dat
metadta <- abund_metadata$metadata

###### Check for Batch Correction ######
