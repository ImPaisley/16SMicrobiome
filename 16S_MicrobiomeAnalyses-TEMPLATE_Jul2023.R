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

###### List of Packages Used ######
library(vegan)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(tidyverse)
library(pgirmess)
library(multcompView)
library(MMUPHin)

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
# OR dat01.bc.dist <- vegdist(dat.ra, method = "bray") is using dat01.per abundance

## next, determine if the batch effect is present and significant
# NOTE: if checking for the batch effect, you must have a column named
# "Batch" in your metadata how the batch is defined is dependent on your project. 
# For example, if you had multiple individuals aid in numerous sequence runs
# then it will be good to base your "batches" on the sequence run that
# each sample was a part of.

## function
batch_test <- function(bray_dist_matrix, metadata){
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

###### Generating Rarefaction Curve #######
library(vegan)
## assigning the batch corrected feature table to its own variable; MAKE SURE TO READ ROWNAMES
rardat<-read.csv("feature_ADJUSTED.csv", header=TRUE, row.names=1, sep=',') 

# samples are in columns and need to be in the rows so we need to flip or transpose the file
# transpose the data to rows; transposing changes the data frame to a matrix!
trans.rardat <- t(rardat)
# check file to make sure it worked 
trans.rardat[1:5,1:5] #shows rows 1 through 5 and the samples should now be the rows
# assign the transformed data matrix into main data variable
rardat <- trans.rardat
# change back into data frame instead of matrix
rardat <-as.data.frame(rardat)
#check data file to make sure it looks okay 
View(rardat)

rowSums(rardat) #sums the value of each row in the data frame; 
# this shows the total sequencing reads for each sample

## Creating the rarefaction curve
# count the number of species within each sample
S <- specnumber(rardat)
raremax <- min(rowSums(rardat)) # takes the sample with the lowest number of sequencing reads

## Plotting the rarefaction curves
# ** auto removes samples that have no reads **

# creating color palette
col <- c("darkred", "forestgreen", "hotpink", "blue")  # feel free to edit the colors
                                                      # keep it between 3 and 5 colors
grp <- factor(sample(seq_along(col), nrow(rardat), replace = TRUE))
cols <- col[grp]

# creating rarefaction curve
# create the curve to estimate where the inflection point (the point at which 
# most of the lines being to plateau) lies, then assign that value to the
# variable "inflection" below
rarecurve(rardat, step = 500, sample=raremax, col = cols, label = TRUE, 
          main="Title", cex= 0.35, cex.axis= 0.95, cex.lab= 1, xlim=c(0,200000), 
          xlab = "# of Sequencing Reads", ylab = "# of ASVs")
inflection <- 10000 # insert your estimated inflection point value
abline(0,0) # creates a vertical line at 0,0
abline(v = inflection, col="black", lwd=1.4) # creates the inflection line at specified value
                                             # feel free to change the line color and thickness (col and lwd)

###### Reproducing Abundance Data using Batch-corrected Data ######
### NOTE: If you ignore metadata, you will NOT be able to perform statistical
###       analyses on your data, only taxonomic analyses.

## RUN FUNCTIONS FIRST!
# use "calculate_abund.metadta" if you have your metadata file
# use "calculate_abund" if ignoring metadata (NO metadata file)
calculate_abund.metadata <- function(feature_csv,metadata_csv){
  library(vegan)
  dat<-t(data.matrix(read.csv(feature_csv, header=TRUE, row.names = 1))) 
  #transposed so that the row names are now the sample names and the ASVs are the columns
  #BUT it is still a MATRIX
  metadata <- read.csv(metadata_csv, header = TRUE, row.names = 1)
  dat <- as.data.frame(dat)
  typeof(dat) #list
  common.rownames <- intersect(rownames(dat), rownames(metadata)) #returns the same row names found in each dataframe
  dat <- dat[common.rownames,] #subsets dat to only include the same row names as metadata
  write.csv(dat, "feature_ADJUSTED-Transposed_matched.csv")
  metadata <- metadata[common.rownames,] #subsets metadata to only include the same row names as dat
  write.csv(dat, "metadata_matched.csv")
  otu.abund<-which(colSums(dat)>2) #removes singletons and doubletons
  dat.dom<-dat[,otu.abund] #include dominant taxa
  dat.pa<-decostand(dat.dom, method ="pa") #turns dat.dom into presence/absence data (1/0)
  dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
  dat.01per<-dat.dom[,dat.otus.01per] #removed ASVs that occur less than 0.1%
  write.csv(dat.01per, "feature_01percent.csv")
  dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
  dat.001per<-dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; increases the number of ASVs - includes more "microdiversity"
  write.csv(dat.001per, "feature_001percent.csv")
  dat.ra<-decostand(dat.01per, method = "total") #relative abundance of >1% taxa
  write.csv(dat.ra, "relative-abundance.csv")
  dfs_to_return <- list(as.data.frame(dat),as.data.frame(metadata),
                        as.data.frame(dat.dom),as.data.frame(dat.pa),
                        as.data.frame(dat.01per),as.data.frame(dat.001per),
                        as.data.frame(dat.ra))
  names(dfs_to_return) <- c("dat", "metadata","dat.dom","dat.pa","dat.01per","dat.001per","dat.ra")
  return(dfs_to_return)
}
calculate_abund <- function(feature_csv){
  library(vegan)
  dat<-t(data.matrix(read.csv(feature_csv, header=TRUE, row.names = 1))) 
  #transposed so that the row names are now the sample names and the ASVs are the columns
  dat <- as.data.frame(dat)
  write.csv(dat, "feature_ADJUSTED-Transposed.csv")
  otu.abund<-which(colSums(dat)>2) #removes singletons and doubletons
  dat.dom<-dat[,otu.abund] #include dominant taxa
  dat.pa<-decostand(dat.dom, method ="pa") #turns dat.dom into presence/absence data (1/0)
  dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
  dat.01per<-dat.dom[,dat.otus.01per] #removed ASVs that occur less than 0.1%
  write.csv(dat.01per, "feature_01percent.csv")
  dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
  dat.001per<-dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; increases the number of ASVs - includes more "microdiversity"
  write.csv(dat.001per, "feature_001percent.csv")
  dat.ra<-decostand(dat.01per, method = "total") #relative abundance of >1% taxa
  write.csv(dat.ra, "relative-abundance.csv")
  dfs_to_return <- list(as.data.frame(dat),as.data.frame(dat.dom),as.data.frame(dat.pa),
                        as.data.frame(dat.01per),as.data.frame(dat.001per),
                        as.data.frame(dat.ra))
  names(dfs_to_return) <- c("dat","dat.dom","dat.pa","dat.01per","dat.001per","dat.ra")
  return(dfs_to_return)
}

## insert your files into the function
# "feature.csv" = insert the name of your ADJUSTED feature table .csv file
# "metadata.csv" = insert the name of your metadata .csv file 

abund_metadata <- calculate_abund.metadata("feature_ADJUSTED.csv","metadata.csv") #if you have metadata
abund <- calculate_abund("feature_ADJUSTED.csv") #if you have NO metadata

# functions will return a list of data frames (must assign to a variable in order to access them):
# dat = your original adjusted feature table that is transposed (or flipped)
# metadata = your metadata table             **dat and metadata should have the same samples!**
# dat.dom = the dominant ASVs found in dat
# dat.pa = the dat.dom table transformed into presence/absence data
# dat.01per = ASVs from dat that have an abundance of 0.1% or higher
# dat.001per = ASVs from dat that have an abundance of 0.01% or higher
# dat.ra = ASVs from dat.01per normalized into relative abundance

## NOTE: The dat, dat.01per, dat.001per, and dat.ra tables are saved as csvs into your working directory

# To save any dataframe (df) within the function result list as a single df in 
# your workspace, simply assign it to its own variable
# Example saving dat and metadata
feature <- abund_metadata$dat
metadta <- abund_metadata$metadata

###### Alpha Diversity - Measures ###### 
## Alpha diversity: the species richness that occurs within a given area within a region
## that is smaller than the entire distribution of the species (Moore, 2013)
## You can use the relative abundance or dat.01per data
adivmeasures <- function(abundance_data) {
  library(vegan)
  # Species richness: the number of species within a region
  S <- as.data.frame(specnumber(abundance_data))
  colnames(S)[1] ="S"
  #No. individuals:
  N <- as.data.frame(rowSums(abundance_data))
  colnames(N)[1] ="N"
  #Shannon-Weiner Diversity:
  ## Shannon index: a measure of the information content of a community rather than of the particular species
  ##               that is present (Moore, 2013) [species richness index]
  ## strongly influenced by species richness and by rare species (so sample size is negligible)
  H <- as.data.frame(diversity(abundance_data), index="shannon")
  colnames(H)[1] ="H"
  #Pielou's Evenness:
  ## Pielou's evenness: an index that measures diversity along with the species richness
  ## Formula - J = H/log(S) (aka Shannon evenness index)
  ## evenness = the count of individuals of each species in an area; 0 is no evenness & 1 is complete evenness 
  J = H/log(S)
  colnames(J)[1] ="J"
  #Simpson's Diversity (1/D) (inverse):
  ## gives the Simpson index the property of increasing as diversity increases (the dominance of
  ## a few species decreases)
  inv.D <- as.data.frame(diversity(abundance_data, index="inv"))
  colnames(inv.D)[1] ="inv.D"
  #Combine data together into a single new data frame, export as CSV
  diversitybysample <- cbind(S, N, H, J,inv.D)
  write.csv(diversitybysample, "AlphaDiversity.csv")
  return(diversitybysample)
}
## This function returns one dataframe that contains the alpha diversity measure metrics for each sample

# Merge diversity with metadata table and export as csv
# DONE IN EXCEL BEFORE RUNNING NEXT CODE: rename the first column of both data 
# tables as the SAME NAME. e.g. both should have the sample column labelled as "Sample"
diversitybysample <- read.csv("AlphaDiversity.csv", row.names = 1)
met <- read.csv("metadata.csv", row.names = 1)
adivmet <- cbind(diversitybysample,met)
write.csv(adivmet,"Metadata-Diversity.csv")

###### Alpha Diversity Statistics ######
library(vegan)
library(ggplot2)

# load in your metadata that has the alpha diversity indices included
metadata <- read.csv("Metadata-Diversity.csv", header = TRUE, row.names = 1)

#### Testing Statistical Significance
## Normality - Shapiro Test (only done on NUMERIC data)
## p <= 0.05 = H0 REJECTED -> DATA IS NOT NORMAL 
## p > 0.05 = H0 ACCEPTED -> DATA IS NORMAL

## If data is NOT normal the first time, try transforming the data using log
## and sqrt and retest for normality.

#Alpha Diversity Variables - Test for Normality
shapiro.test(metadata$S)
shapiro.test(metadata$N)
shapiro.test(metadata$H)
shapiro.test(metadata$J)
shapiro.test(metadata$inv.D)

### CONTINUE HERE IF DATA IS NORMAL (OR TRANSFORMATIONS NORMALIZED THE DATA)
# ANOVA: Parametric Data (normal)
# Tukey Test - calculate pairwise comparisons between group levels

# RUN FUNCTION FIRST!!
p_anova <- function(adiv_var, test_var) {
  library(pgirmess)
  library(multcompView)
  anova_res <- aov(adiv_var ~ test_var) # performs an ANOVA 
  summary <- summary(anova_res) # provides an ANOVA table with p-values
  Tukey <- TukeyHSD(anova_res) # pairwise comparison
  return_items <- list(summary, Tukey)
  names(return_items) <- c("ANOVA summary", "TukeyTest")
  return(return_items)
}
## insert your files into the function
# "adiv_var" = insert the alpha diversity metric you want to test from the metadata file
#              (metadata$S, metadata$N, metadata$J, metadata$H, or metadata$inv.D)
# "test_var" = insert the variable(s) you want to test against the alpha diversity
#              metric (e.g. year (metadata$year), month, year and month (metadata$year+month), etc.)
#         ** MAKE SURE YOUR VARIABLE IS A FACTOR! IF NOT THEN SURROUND IT WITH as.factor() **
p_anova(adiv_var, test_var) # or p_anova(adiv_var, as.factor(test_var))

# Example inputs:
S_anova <- p_anova(metadata$S, as.factor(metadata$Year))
N_anova <- p_anova(metadata$N, metadata$Year)
J_anova <- p_anova(metadata$J, metadata$Year)
H_anova <- p_anova(metadata$H, metadata$Year)
inv.D_anova <- p_anova(metadata$inv.D, metadata$Year)

### CONTINUE HERE IF DATA IS NOT NORMAL AND TRANSFORMATIONS DID NOT WORK
library(pgirmess)
library(multcompView)

# Kruskal Wallis: Nonparametric Data (not normal)
# Pairwise Wilcox Test - calculate pairwise comparisons between group levels 
#                        with corrections for multiple testing (non-parametric)

# RUN FUNCTION FIRST!!
nonp_kruskal <- function(adiv_var, test_var) {
  library(pgirmess)
  library(multcompView)
  kruskal.test(adiv_var ~ test_var)
  pair_WilTest <- pairwise.wilcox.test(adiv_var, test_var, p.adjust.method = "fdr") #pairwise comparisons between the variable levels
  kmc <- kruskalmc(adiv_var ~ test_var) # multiple-comparison test
  # comparisons TRUE= significantly different or FALSE= not significantly different
  # To look for homogeneous groups, and give each group a code (letter):
  test <- kmc$dif.com$stat.signif # select logical vector
  names(test) <- row.names(kmc$dif.com)# add comparison names
  # create a list with "homogeneous groups" coded by letter
  let <- multcompLetters(test, compare="<", threshold=0.05,
                         Letters=c(letters, LETTERS, "."),reversed = FALSE)
  # significant letters for the multiple comparison test
  # if the letter are the SAME, then no significant differences were found
  # between those variables
  returned_items <- list(pair_WilTest,kmc,let)
  names(returned_items) <- c("pairwise", "multiComp","letter-comparisons")
  return(returned_items)
}
## insert your files into the function
# "adiv_var" = insert the alpha diversity metric you want to test from the metadata file
#              (metadata$S, metadata$N, metadata$J, metadata$H, or metadata$inv.D)
# "test_var" = insert the variable(s) you want to test against the alpha diversity
#              metric (e.g. year (metadata$year), month, year and month (metadata$year+month), etc.)
#         ** MAKE SURE YOUR VARIABLE IS A FACTOR! IF NOT THEN SURROUND IT WITH as.factor() **
nonp_kruskal(adiv_var, test_var) # or nonp_kruskal(adiv_var, as.factor(test_var))

# Example inputs:
S_krusk <- nonp_kruskal(metadata$S, metadata$Year)
N_krusk <- nonp_kruskal(metadata$N, metadata$Year)
J_krusk <- nonp_kruskal(metadata$J, metadata$Year)
H_krusk <- nonp_kruskal(metadata$H, metadata$Year)
inv.D_krusk <- nonp_kruskal(metadata$inv.D, metadata$Year)

### Plotting boxplots of alpha diversity by specified variable
## NOTES: You can replace "Year" with your specified variable
##        Adding text to your graph is OPTIONAL but if you are adding it then
##        you'll have to play around with their coordinates, what they say, size, and color

# Creating pdf for the plots to populate
pdf("AlphaDiverisityPlots.pdf")
# plot each boxplot on its own page
par(mar=c(5,6,2,2)+0.1)
boxplot(S~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Species Richness (S)", line=4.25, cex.lab=1.15)
text(y=1500, x=3, labels="b", col="blue", cex=1.2)
text(y=1420, x=2, labels="a", col="red", cex=1.2)        # labeling which groups are significantly different than the other 
text(y=1585, x=1, labels="a", col="red", cex=1.2)

par(mar=c(5,4.5,2,2)+0.1)
boxplot(H~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Shannon Diversity Index (H)", line=2.8, cex.lab=1.15)
text(y=4, x=3, labels="b", col="blue", cex=1.2)
text(y=3.4, x=2, labels="a", col="red", cex=1.2)        
text(y=3.6, x=1, labels="ab", col="purple", cex=1.2)

boxplot(J~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Species Evenness (J)", line=3, cex.lab=1.15)
text(y=0.73, x=3, labels="b", col="blue", cex=1.2)
text(y=0.685, x=2, labels="ab", col="purple", cex=1.2)        
text(y=0.73, x=1, labels="a", col="red", cex=1.2)

par(mar=c(5,6,2,2)+0.1)
boxplot(inv.D~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="inverse Simpson Diversity Index (inv.D)", line=3.6, cex.lab=1.15)
text(y=440, x=3, labels="a", col="red", cex=1.2)
text(y=420, x=2, labels="b", col="blue", cex=1.2)        
text(y=340, x=1, labels="a", col="red", cex=1.2)

boxplot(N~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="No. of Individuals (N)", line=4.25, cex.lab=1.15)
text(y=90000, x=3, labels="b", col="blue", cex=1.2)
text(y=130000, x=2, labels="a", col="red", cex=1.2)        
text(y=160000, x=1, labels="a", col="red", cex=1.2)

# stop saving to pdf 
dev.off()