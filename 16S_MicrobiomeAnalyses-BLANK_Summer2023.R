## 16S Microbiome Analyses - Molecular Microbiology & Genomics Lab @ Nova Southeastern University ##
## This R-Script includes routine microbiome analyses conducted using
## the 16S rRNA gene (V4 region (primers 515F/806R)).
## Initial data cleaning starts from QIIME2, so you will need the following 
## files from the QIIME2 resulting files in your desired working folder: 
## feature.csv and taxonomy.csv (both are originally tsv files but use Excel to convert to csv)
## Feel free to change anything in this script to fit your project analyses!
## Please make sure to keep this original copy and to save any edited versions
## with your project data!

## This is a link to the entire color chart used by R: https://r-charts.com/colors/
## You can also run the code: colors() and R will give you a list of the color names.
## I personally use the link since you can see the colors.

              ## Created by Paisley Samuel - Summer 2023 ##

###### Table of Contents/Order of Processes ######
### 1. List of packages that are used throughout the script 
### 2. Setting the working directory and seed 
### 3. Producing the abundance data using raw sequencing data
### 4. Checking the produced abundance data for the batch effect and correcting if necessary
### 5. Generating a rarefaction curve of the original/adjusted abundance data
### 6. Reproducing the abundance data using the batch-corrected sequencing data
### 7. Taxonomy analyses on abundance data
### 8. Creating a table of calculated alpha diversity metrics
### 9. Performing statistical analyses between metadata variables and alpha diversity metrics
### 10. Creating a Bray-Curtis distance matrix and performing beta diversity statistical analyses
### 11. Creating nMDS plots using BC distance matrix 
### 12. Performing CCA analysis and creating CCA plots using 

###### List of Packages Used ######
## Here is a list of all the packages used throughout this R script.
## Take the necessary steps to install any packages that you do not have
## installed already on your computer (may have to use install.packages()
## or Bioconductor depending on the package)

library(ggplot2)
library(MBiome16S) # created specifically for 16S V4 region analyses in Lopez lab
library(microbiome)
library(MMUPHin)
library(multcompView)
library(pgirmess)
library(phyloseq)
library(psych)
library(tidyverse)
library(vegan)

## Un-comment (i.e. erase the "#") and run the code below if any of the packages listed above needs to
## be installed 

##ggplot2
# install.packages("ggplot2")

##MBiome16S
# You will find a .tar.gz file labelled MBiome16S.tar.gz on the desktop of the 
# lab computer. Go to the 'Packages' tab in the window of RStudio in the bottom-right
# and click on 'Install' on the top-right. In the dropdown labelled 'Install from',
# select 'Package Archive File...' and navigate to the desktop where the tar.gz file is.
# Select the MBiome16S.tar.gz file and hit install.

##microbiome
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("microbiome")

##MMUPHin
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("MMUPHin")

##multcomView
# install.packages("multcompView")

##pgirmess
# install.packages("pgirmess")

##phyloseq
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

##psych
# install.packages("psych")

##tidyverse
# install.packages("tidyverse")

##vegan
# install.packages("vegan")


###### Set Working Directory and Seed - ALWAYS RUN THIS FIRST! ######
setwd_seed("file_path", seed_number)
# "file_path" = insert the file path of your desired folder where all your files will be stored
# seed_number = insert a random number


###### Check for Batch Correction ######
## make sure you run the abundance function first so you have your abundance
## tables to input into the following functions (again RUN FUNCTIONS FIRST!)

## Batch correction

## Determine if the batch effect is present and significant
# NOTE: if checking for the batch effect, you must have a column named
# "Batch" in your metadata how the batch is defined is dependent on your project.
# For example, if you had multiple individuals aid in numerous sequence runs
# then it will be good to base your "batches" on the sequence run that
# each sample was a part of.

batch_test(abundance, metadata)

# "abundance" = insert the abundance table you just created (dat.ra, dat.01per, etc.)
# "metadata" = insert your metadata variable (may be already named 'metadata')


## If p <= 0.05, then the batch effect was found to be significant and you MUST
## correct the batch effect BEFORE moving on to further analyses



## CONTINUE HERE IF YOUR DATA IS SIGNIFICANT FOR BATCH EFFECT!! SKIP IF NOT!
# Using a package called MMUPHin, we will have to adjust the data so that the
# batch effect is no longer affecting the statistical outcome of the data.

batchadj <- batch_correct("feature.csv","metadata.csv")

# "feature.csv" = insert the name of your feature table .csv file
# "metadata.csv" = insert the name of your metadata .csv file
# The adjusted feature table is saved as a csv under your working directory as "feature_ADJUSTED.csv"


# function will return a list of data frames:
# dat = your original feature table that is transposed (or flipped)
# metadata = your metadata table             **dat and metadata should have the same samples!**
# adj-feature = the adjusted feature table that is now batch corrected

# YOU NEED TO USE THE BATCH CORRECTED FEATURE TABLE FOR THE REST OF THE
# ANALYSES IF YOUR DATA WAS BATCH CORRECTED!!


###### Generating Rarefaction Curve #######
library(vegan)
## assigning the batch-corrected or original feature table to its own variable; MAKE SURE TO READ ROWNAMES
rardat<-read.csv("feature_ADJUSTED.csv", header=TRUE, row.names=1, sep=',')
# OR
rardat<-read.csv("feature.csv", header=TRUE, row.names=1, sep=',')

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

# Creating rarefaction curve
# create the curve to estimate where the inflection point (i.e. the point at which
# most of the lines being to plateau) lies, then assign that value to the
# variable "inflection" below
rarecurve(rardat, step = 500, col = cols, label = TRUE,
          main="Title", cex= 0.35, cex.axis= 0.95, cex.lab= 1, xlim=c(0,200000),
          xlab = "# of Sequencing Reads", ylab = "# of ASVs")
inflection <- 10000   # insert your estimated inflection point value
abline(v = inflection, col="black", lwd=1.4) # creates the inflection line at specified value
                                             # feel free to change the line color and thickness (col and lwd)

###### Producing Relative Abundance Data ######
### NOTE: If you ignore metadata, you will NOT be able to perform statistical
###       analyses on your data, only taxonomic analyses.

# use "calculate_abund.metadata" if you have your metadata file
# use "calculate_abund" if ignoring metadata (NO metadata file)

abund_metadata <- calculate_abund.metadata("feature.csv","metadata.csv") #if you have metadata
abund <- calculate_abund("feature.csv") #if you have NO metadata


# OR


# use these functions if you are reproducing the abundance data on your batch-
# corrected data
# use "calculate_abund.metadata_ADJ" if you have your metadata file
# use "calculate_abund_ADJ" if ignoring metadata (NO metadata file)

abund_metadata_ADJ <- calculate_abund.metadata_ADJ("feature_ADJUSTED.csv","metadata.csv") #if you have metadata
abund_ADJ <- calculate_abund_ADJ("feature_ADJUSTED.csv") #if you have NO metadata


# "feature_ADJUSTED.csv" = insert the name of your ADJUSTED feature table .csv file
# "feature.csv" = insert the name of your feature table .csv file
# "metadata.csv" = insert the name of your metadata .csv file


# functions will return a list of data frames:
# dat = your original feature table that is transposed (or flipped)
# metadata = your metadata table             **dat and metadata should have the same samples!**
# dat.dom = the dominant ASVs found in dat
# dat.pa = the dat.dom table transformed into presence/absence data
# dat.01per = ASVs from dat that have an abundance of 0.1% or higher
# dat.001per = ASVs from dat that have an abundance of 0.01% or higher
# dat.ra = ASVs from dat.01per normalized into relative abundance

# To save any dataframe (df) within the function result list as a single df in
# your workspace, simply assign it to its own variable
# Example saving dat and metadata
feature <- abund_metadata$dat
metadta <- abund_metadata_ADJ$metadata


###### Taxonomy Analyses ######
## Creating a phyloseq object
physeq_objs <- create_phyloseq(abundance, "taxonomy.csv", "metadata.csv")

# abundance = insert the abundance table you created (dat.ra or dat.01per)
# "taxonomy.csv" = insert your taxonomy .csv file
# "metadata.csv" = insert your metadata .csv file


## Basic sequencing statistics (USING THE "physeq" OBJECT IN YOUR PHYLOSEQ OBJECT LIST VARIABLE)
# First, assign physeq and physeq_transform to their own variables
physeq <- physeq_objs$physeq
physeq_transform <- physeq_objs$physeq_transform

# Check number of sequencing reads observed in each sample
sample_sums(physeq)

# Calculate the total sequencing reads of your data
sum(sample_sums(physeq))

# Calculate the average number sequencing reads
mean(sample_sums(physeq))

# Find the lowest number of sequencing reads in your data
min(sample_sums(physeq))

# Find the highest number of sequencing reads in your data
max(sample_sums(physeq))

# Calculate the standard deviation of sequencing reads between samples
sd(sample_sums(physeq))

# Find the total amount of ASVs within your data
ntaxa(physeq)

## Retrieves the unique taxonomic ranks observed in the data set
## [#] = rank (starting from Domain and onward DPCOFGS)
get_taxa_unique(physeq, taxonomic.rank=rank_names(physeq)[7], errorIfNULL=TRUE)
#Unique Domains [1] =
#Unique Phyla [2] =
#Unique Classes [3] =
#Unique Orders [4] =
#Unique Families [5] =
#Unique Genera [6] =
#Unique Species [7] =

## Aggregating by Taxonomic level
# This function allows you to aggregate the taxonomy based on Taxonomic
# level, gives you the counts of each level, and saves as a CSV file

agg.tax.level <- taxon_aggregate(phyloseq_object,"tax_level")

# phyloseq_object = insert the transformed phyloseq object you created (physeq_transform)
# "tax_level" = insert the taxonomic level you will like to aggregate on

# For example to aggregate by Class, the code would read:
# Class <- taxon_aggregate(physeq_transform,"Class")


## Calculating the top taxonomic levels
# First, find the top taxa names and abundances
# This function allows you to specify the top taxonomy of your data and
# exports the resulting vector as a csv file

top.n.names <- MBiome16S::top_taxa(phyloseq_object, "tax_level", n)

# phyloseq_object = insert the transformed phyloseq object you created (physeq_transform)
# "tax_level" = insert the taxonomic level you will like to aggregate on
# n = insert the number of taxa you want (e.g. top 5, top 10, etc.)

# For example for top 5 Classes, the code would read:
# top.class.names <- top_taxa(physeq_transform,"Class",5)


# Next, subset your phyloseq object to only include the top taxa you just specified
# cuts down the phyloseq object to only the top n
# Replace the "taxonomic level" with the level you want WITHOUT QUOTES!
top_tax <- subset_taxa(physeq_transform, "taxonomic level" %in% names(top.n.names))
# For example for Class, the code would read:
# top_Class <- subset_taxa(physeq_transform, Class %in% names(top.class.names))

## Creating stacked bar plots from your phyloseq object subset
# Change any taxonomic names to the level that you need ("fill=" and "colour=")

library(ggplot2)
topTaxplot <- plot_bar(top_tax, x="Sample", y="Abundance", fill="Class")
topTaxplot <- topTaxplot +
  geom_bar(aes(fill=Class, colour=Class), stat="identity", position="fill", width = 0.9) +   #width=0.96 removes any space between bars
  ggtitle("Top 5 Classes") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+   #vjust= moves the x-axis text labels
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5)) +   #hjust= 0.5 centers the title
  theme(legend.title = element_text(face="italic"))
topTaxplot # populates the graph into the "plot" panel to the right

###### Alpha Diversity - Measures ######
## Alpha diversity: the species richness that occurs within a given area within a region
## that is smaller than the entire distribution of the species (Moore, 2013)
## You can use the relative abundance or dat.01per data

adivmeasures(abundance_data)
## This function returns one dataframe that contains the alpha diversity 
## measure metrics for each sample


# Merge diversity with metadata table and export as csv
# DONE IN EXCEL BEFORE RUNNING NEXT CODE: rename the first column of both data
# tables as the SAME NAME. e.g. both should have the sample column labelled as "Sample"
diversitybysample <- read.csv("AlphaDiversity.csv", row.names = 1)
met <- read.csv("metadata.csv", row.names = 1)
adivmet <- cbind(diversitybysample,met)
write.csv(adivmet,"Metadata-Diversity.csv")

###### Alpha Diversity Statistics ######
library(vegan)
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
p_anova(adiv_var, test_var) # or p_anova(adiv_var, as.factor(test_var))

# "adiv_var" = insert the alpha diversity metric you want to test from the metadata file
#              (metadata$S, metadata$N, metadata$J, metadata$H, or metadata$inv.D)
# "test_var" = insert the variable(s) you want to test against the alpha diversity
#              metric (e.g. year (metadata$year), month, year and month (metadata$year+month), etc.)

#     ** MAKE SURE YOUR VARIABLE IS A FACTOR! IF NOT THEN SURROUND IT WITH as.factor() **

# Example inputs:
S_anova <- p_anova(metadata$S, as.factor(metadata$Year))
N_anova <- p_anova(metadata$N, metadata$Year)
J_anova <- p_anova(metadata$J, metadata$Year)
H_anova <- p_anova(metadata$H, metadata$Year)
inv.D_anova <- p_anova(metadata$inv.D, metadata$Year)

### CONTINUE HERE IF DATA IS NOT NORMAL AND TRANSFORMATIONS DID NOT WORK
# Kruskal Wallis: Nonparametric Data (not normal)
# Pairwise Wilcox Test - calculate pairwise comparisons between group levels
#                        with corrections for multiple testing (non-parametric)

nonp_kruskal(adiv_var, test_var) # or nonp_kruskal(adiv_var, as.factor(test_var))

# "adiv_var" = insert the alpha diversity metric you want to test from the metadata file
#              (metadata$S, metadata$N, metadata$J, metadata$H, or metadata$inv.D)
# "test_var" = insert the variable(s) you want to test against the alpha diversity
#              metric (e.g. year (metadata$year), month, year and month (metadata$year+month), etc.)
#         ** MAKE SURE YOUR VARIABLE IS A FACTOR! IF NOT THEN SURROUND IT WITH as.factor() **

# Example inputs:
S_krusk <- nonp_kruskal(metadata$S, metadata$Year)
N_krusk <- nonp_kruskal(metadata$N, metadata$Year)
J_krusk <- nonp_kruskal(metadata$J, metadata$Year)
H_krusk <- nonp_kruskal(metadata$H, metadata$Year)
inv.D_krusk <- nonp_kruskal(metadata$inv.D, metadata$Year)

# Looking into the Kruskal-Wallis results for each diversity metric
# The important results to be noted are the KW p-values and the letter comparisons.
# You can also take a look at the individual comparison p-values by looking into the pairwise table.

S_krusk$`Kruskal-Wallis Results` ## Kruskal-Wallis chi-squared = 33.035, df = 2, p-value = 6.708e-08
S_krusk$pairwise
S_krusk$`letter-comparisons` 

### You would then repeat these steps for the 4 other alpha diversity metrics!


### Plotting boxplots of alpha diversity by specified variable
## NOTES: You can replace "Year" with your specified variable
##        Adding text to your graph is OPTIONAL but if you are adding it then
##        you'll have to play around with their coordinates, what they say, size, and color.

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

###### Beta Diversity - Statistics ######


betadiv <- betadiv_stats(abundance_data,metadata$Variable)

# "abundance_data" = insert the abundance data frame you created
# "metadata$Variable" = insert your metadata variable with the Variable you want to test (may be already named 'metadata')


###### Beta Diversity - nMDS plots ######
# Creating nMDS plot - 2D
nmds2d <- metaMDS(betadiv$distance_matrix,k=2,autotransform = F,trymax=20) # you will have to extract bc.dist from the previous analysis (object betadiv)
nmds2d
#Dimensions = 2
#Stress =
stressplot(nmds2d)
#Shepard plot shows scatter around the regression between the inter-point
#distances in the final configuration (i.e., the distances between each pair of communities)
#against their original dissimilarities
nmds.plot <- ordiplot(nmds2d,display="sites") #populates the plot into a variable
## Adding ellipses to group years
ordihull(nmds.plot,groups=metadata$Variable,draw="lines",col=c("tomato3","steelblue3","springgreen3")) # adds ellipses around the point groups (OPTIONAL!)
##adjust colors to match each year, pch=20 makes it bullet points
points(nmds.plot,"sites", pch=20, col= "tomato4", select = metadata$Variable == "Level 1")     # set the metadata to the variable you need it to be
points(nmds.plot,"sites", pch=20, col= "steelblue4", select = metadata$Variable == "Level 2")  # and set a different color to each level of the variable
points(nmds.plot,"sites", pch=20, col= "springgreen4", select = metadata$Variable == "Level 3")
##Add Stress Value
text(1.2,1.5,"2D Stress: ", cex=0.9) # make sure you add the stress value in an empty portion of the graph
##Adding legend
legend("topleft",legend= c("Level 1","Level 2", "Level 3"),   # customize the legend to match the colors and variables
       title = "Variable",
       col=c("tomato4","steelblue4","springgreen4"),
       pch=19, cex=1)
##Adding title
title(main="nMDS of Relative Abundances by Variable") # adds a title to the graph

## If making plots for multiple variables, you'll have to redo the point
## customization and legend customization to match the variable

###### CCA analysis & plots ######
library(vegan)
# you can use the dat.ra table or the dat.01per table
# input the quantitative variables (columns that are numeric) from the metadata
ccamodel <- cca(dat.ra~., metadata[,c(1:5)])
finalmodel<- ordistep(ccamodel, scope=formula(ccamodel))
vif.cca(finalmodel) ## values should be under 10
# If VIF>10, the variable presents colinearity with another or other variables.
# In that case, delete the variable from initial dataset and redo the analysis.
# VIF = 1 for completely independent variables,and values above 10 or 20
# (depending on your taste) are regarded as highly multicollinear (dependent on other variables).

# Test the significance of the entire CCA model
anova.cca(finalmodel) #should be significant! (p < 0.05)


finalmodel
## Note that "Total Inertia" is the total variance in species (observations matrix) distributions.
## "Constrained Inertia" is the variance explained by the environmental variables (gradients matrix).
## The "Proportion" values represent the percentages of variance of species distributions explained
## by Constrained (environmental) and Unconstrained variables. Eigenvalues of constrained and
## unconstrained axes represent the amount of variance explained by each CCA axis (graphs usually
## present the first two constrained axes, so take a look at their values).

#Total Inertia = total variance in species (observed distributions)
#Unconstrained Inertia = the variance explained by the environmental variables



R2.adj.cca <- RsquareAdj(finalmodel)
# adjusting the R-squared value: The adjusted R2 tells you the percentage of
# variation explained by only the independent variables that actually affect
# the dependent variable; also indicates how well terms fit a curve or line,
# but adjusts for the number of terms in a model
R2.adj.cca
# r.squared: 
# adj.r.squared:   -> the R2 value you document

summary(finalmodel)

## Correlation between the variables within the model
# Creates pairs plot to see the correlation statistics between each variable
library(psych)
pairs.panels(metadata[,c(1:5)]) # use the metadata columns that are
                                # in your final CCA model

### Creating CCA plots
cca.p <- plot(finalmodel,type = "none")

# Fitting of the environmental variables to the CCA plot
ef.cca <- envfit(cca.p,metadata[,c(1:5)]) # use the metadata columns that are
                                          # in your final CCA model

# Creating R2 threshold for vectors
ef.cca <- select.envfit(ef.cca, 0.3) # only includes significant variables
                                     # with an R2 of 0.3 or higher (feel free to change the 0.3!)

# Setting up base plot
par(mar=c(5.1, 6.1, 3.1, 4.1))
plot(finalmodel,type = "none")
abline(h = 0, v = 0, col = "white", lwd = 2)
box()

# EDIT POINTS AS IT MATCHES YOUR ANALYSIS!
# Adding the points
points(cca.p,"sites", pch=19, col= "goldenrod3", select = metadata$Variable == "Level 1")
points(cca.p,"sites", pch=19, col= "mediumpurple2", select = metadata$Variable == "Level 2")
points(cca.p,"sites", pch=19, col= "springgreen4", select = metadata$Variable == "Level 3")
# Plotting envfit vectors
plot(ef.cca, col = "black", p.max=0.05)
# Add legend & Title
legend(locator(1),legend=c("Level 1","Level 2", "Level 3"),
       col=c("goldenrod3","mediumpurple2", "springgreen4"), pch=19, cex=1.2,
       title = "Variable")
# locator(1) allows you to choose a place to place the legend within your plot
title(main="Title")





