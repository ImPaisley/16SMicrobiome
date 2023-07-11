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
metadata <- abund_metadata$metadata

###### Check for Batch Correction ######
## make sure you run the abundance function first so you have your abundance
## tables to input into the following functions (again RUN FUNCTIONS FIRST!)

## Batch correction
## first, create Bray-Curtis dissimilarity distance matrix
## you can either use the dat.01per or dat.ra tables as input,
## they give the same result

ra.bc.dist <- vegdist(dat.ra, method = "bray")
#OR dat01.bc.dist <- vegdist(dat.ra, method = "bray") is using dat01.per abundance

## next, determine if the batch effect is present and significant
# NOTE: if checking for the batch effect, you must have a column named
# "Batch" in your metadata how the batch is defined is dependent on your project. 
# For example, if you had multiple individuals aid in numerous sequence runs
# then it will be good to base your "batches" on the sequence run that
# each sample was a part of.

## function
batch_test <-function(bray_dist_matrix, metadata){
  library(vegan)
  dis.Batch <- betadisper(bray_dist_matrix,metadata$Batch) # betadisper calculates dispersion (variances) within each group 
  test <- permutest(dis.Batch, pairwise=TRUE, permutations=999) #determines if the variances differ by groups
  if (test$tab$`Pr(>F)`[1] <= 0.05){    #differences are SIGNIFICANT - use ANOSIM
    ano_sim <- anosim(ra.bc.dist, metadata$Batch, permutations = 999)
    return(ano_sim)
  }
  else{            #differences are NOT SIGNIFICANT - use PERMANOVA (adonis))
    p_anova <- adonis2(ra.bc.dist~metadata$Batch, permutations = 999)
    return(p_anova)
  }
}

## insert your input into the function
# "bray_dist_matrix" = insert the distance matrix you just created 
# "metadata" = insert your metadata variable (may be already named 'metadata') 
batch_test(ra.bc.dist, metadata)

## If p <= 0.05, then the batch effect was found to be significant and you MUST
## correct the batch effect BEFORE moving on to further analyses

## CONTINUE HERE IF YOUR DATA IS SIGNIFICANT FOR BATCH EFFECT!! SKIP IF NOT!
# Using a package called MMUPHin, we will have to adjust the data so that the 
# batch effect is no longer affecting the statistical outcome of the data.

batch_correct <- function(feature_csv, metadata_csv) {
  library(vegan)
  library(MMUPHin)
  dat <- t(data.matrix(read.csv(feature_csv, header = TRUE, row.names = 1)))
  metadata <- read.csv(metadata_csv, header = TRUE, row.names = 1)
  dat <- as.data.frame(dat)
  typeof(dat)
  common.rownames <- intersect(rownames(dat), rownames(metadata))
  dat <- dat[common.rownames, ]
  metadata <- metadata[common.rownames, ]
  #Adjusting (removing) batch effect
  fit_adjust_batch <- adjust_batch(feature_abd = t(dat), # ASVs should be rows in feature table (MATRIX)
                                   batch = "Batch",
                                   data = metadata)   # samples should be rows in metadata (DATAFRAME)
  feat_abd_adj <- fit_adjust_batch$feature_abd_adj #now adjusted feature table MATRIX
  feat_abd_adj <- as.data.frame(feat_abd_adj) #converting to data frame
  write.csv(feat_abd_adj, "feature_ADJUSTED.csv") #saving as csv
  dfs_to_return <- list(as.data.frame(dat), as.data.frame(metadata),
                        as.data.frame(feat_abd_adj))
  names(dfs_to_return) <- c("dat", "metadata", "adj-feature")
  return(dfs_to_return)
}

## insert your files into the function
# "feature.csv" = insert the name of your feature table .csv file
# "metadata.csv" = insert the name of your metadata .csv file 
# The adjusted feature table is saved as a csv under your working directory as "feature_ADJUSTED.csv"

batchadj <- batch_correct("feature.csv", "metadata.csv")

# function will return a list of data frames:
# dat = your original feature table that is transposed (or flipped)
# metadata = your metadata table             **dat and metadata should have the same samples!**
# adj-feature = the adjusted feature table that is now batch corrected

# YOU NEED TO USE THE BATCH CORRECTED FEATURE TABLE FOR THE REST OF THE
# ANALYSES IF YOUR DATA WAS BATCH CORRECTED!!






