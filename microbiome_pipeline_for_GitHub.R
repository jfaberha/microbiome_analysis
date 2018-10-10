#############################################################################
#   Microbiome analysis for brooding, release, fed, and starved A. burtoni  #
#            Last updated 10/10/2018 by Josh Faber-Hammond                  #
# See associated manuscript for more details about the analysis and results #
#############################################################################

# Load libraries
library(MASS)
library(plyr)
library(tidyr)
#install.packages("tidyr")
library(dplyr)
#install.packages("dplyr")
library(phangorn)
#install.packages("phangorn")
library(vegan)
#biocLite("vegan")
library("DESeq2")
#packageVersion("DESeq2")
#biocLite("DESeq2")
library(phyloseq)
#packageVersion("phyloseq")
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite(c("phyloseq"))
library(ggplot2)
library("scales")
library("grid")
library(tibble)
#install.packages('tibble', dependencies=TRUE, type="source")
library(FSA)
# install.packages("FSA")
library (labdsv)
# install.packages("labdsv")
library(magrittr)
library (rcompanion)
# install.packages("rcompanion")
library  (multcomp)
# install.packages("multcomp")
require(reshape2)
library(ape)
library(devtools)
#install.packages("devtools")
#source("http://bioconductor.org/biocLite.R")
#biocLite("ggtree")
library("ggtree")
#biocLite("ggtree")
library(plotrix)
#install.packages("plotrix")
#devtools::install_github('reptalex/phylofactor')
library(phylofactor)
require(phylofactor)
library(phytools)
#install.packages("phytools")
??PhyloFactor ### for help topics
#install.packages("inflection")
library(inflection)
#install.packages("inflection")
require(inflection)
#install.packages("chngpt")
library(chngpt)
library(Hmisc)
#install.packages("Hmisc")
library(mgcv)

theme_set(theme_bw())
set.seed(42)

ste <- function(x) sd(x)/sqrt(length(x))

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

#' @title
#' Find the Elbow in a Curve
#'
#' @description
#' This utility function finds the elbow in a curve which is concave
#' relative to a line drawn between the first and last points.
#' The elbow is defined as the point with the greatest
#' orthogonal distance from that line.
#'
#' @param y Numeric vector of y values for the curve.
#' 
#' @param plot Logical. Should a plot be made?
#' 
#' @param returnIndex Logical. Should the return value
#' be the index of the elbow point? 
#' 
#' @return If \code{returnIndex = TRUE}, the index of
#' the elbow point.  If \code{returnIndex = FALSE},
#' a data frame containing an index values (x),
#' the y values passed to the function, and the
#' the orthogonal distances of the y values from
#' the line connecting the first and last points.
#' \code{which.max(data_frame_name$dist)} will give the index of 
#' the elbow point.  
#'
#' @references The concept of this function is based on the
#' clever idea in the
#' Stackoverflow post at stackoverflow.com/a/2022348/633251
#' and relies on code posted at
#' paulbourke.net/geometry/pointlineplane/pointline.r
#' (referenced in the SO post).  Minor modifications
#' to the code were made to that code in order to vectorize it.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @export
#'
#' @importFrom stats lm coef
#'
#' @section Warning:
#' This function makes some simple checks that the data is concave
#' as defined above.  Even so, it may give
#' answers in some cases that are not valid.  Please check
#' on typical data that you encounter to verify that it works
#' in your cases.
#' 
#' @examples
#' tmp <- findElbow(c(0.9, 1.1, 1.1, 1.9, 2.5, 2.8, 4.9, 8.5),
#' 	plot = TRUE) # wandering
#' tmp <- findElbow(c(0.9, 1.0, 1.2, 1.3, 1.5, 1.5, 10.0, 22.0),
#' 	plot = TRUE) # late rise
#' tmp <- findElbow(c(2, 4, 6, 8, 10, 12, 14, 16)^2,
#' 	plot = TRUE) # gradual, no obvious break
#' 
#' # Not the usual way to choose the number of PCs:
#' library("chemometrics")
#' data(glass)
#' pca <- prcomp(glass)
#' eigensum <- sum(pca$sdev * pca$sdev)
#' vv <- 100 * (pca$sdev * pca$sdev/eigensum)
#' cs <- cumsum(vv)
#' tmp <- findElbow(vv, plot = TRUE)
#' tmp <- findElbow(cs, plot = TRUE)
#'

findElbow <- function(y, plot = FALSE, returnIndex = TRUE) {
  
  # The following helper functions were found at
  # paulbourke.net/geometry/pointlineplane/pointline.r
  # via the SO reference below.  The short segment check
  # was modified to permit vectorization.
  
  ##========================================================
  ##
  ##  Credits:
  ##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
  ##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
  ##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
  ##  With grateful thanks for answering our needs!
  ##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
  ##
  ##========================================================
  
  distancePointLine <- function(x, y, slope, intercept) {
    ## x, y is the point to test.
    ## slope, intercept is the line to check distance.
    ##
    ## Returns distance from the line.
    ##
    ## Returns 9999 on 0 denominator conditions.
    x1 <- x-10
    x2 <- x+10
    y1 <- x1*slope+intercept
    y2 <- x2*slope+intercept
    distancePointSegment(x,y, x1,y1, x2,y2)
  }
  
  distancePointSegment <- function(px, py, x1, y1, x2, y2) {
    ## px,py is the point to test.
    ## x1,y1,x2,y2 is the line to check distance.
    ##
    ## Returns distance from the line, or if the intersecting point on the line nearest
    ## the point tested is outside the endpoints of the line, the distance to the
    ## nearest endpoint.
    ##
    ## Returns 9999 on 0 denominator conditions.
    lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
    ans <- NULL
    ix <- iy <- 0   # intersecting point
    lineMag <- lineMagnitude(x1, y1, x2, y2)
    if(any(lineMag < 0.00000001)) { # modified for vectorization by BAH
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)
    if(any((u < 0.00001) || (u > 1))) { # BAH added any to vectorize
      ## closest point does not fall within the line segment, take the shorter distance
      ## to an endpoint
      ix <- lineMagnitude(px, py, x1, y1)
      iy <- lineMagnitude(px, py, x2, y2)
      if(ix > iy)  ans <- iy
      else ans <- ix
    } else {
      ## Intersecting point is on the line, use the formula
      ix <- x1 + u * (x2 - x1)
      iy <- y1 + u * (y2 - y1)
      ans <- lineMagnitude(px, py, ix, iy)
    }
    ans
  }
  
  # End of helper functions by PB
  
  ### Now for the actual findElbow function!
  
  # Find the elbow using the method described in
  # stackoverflow.com/a/2022348/633251
  # but translated to R (see above).
  
  # Add an index to argument values for easy plotting
  DF <- data.frame(x = 1:length(y), y = y)
  fit <- lm(y ~ x, DF[c(1,nrow(DF)),]) # 2 point 'fit'
  m <- coef(fit)[2]
  b <- coef(fit)[1]
  
  # Check to make sure the data is concave as described
  # in the documentation, as arbitrary trends could give
  # misleading answers.  The following approach simply
  # checks to make sure all values are either above or
  # below the reference line.  This allows the values
  # to vary quite a bit and still return an answer.
  
  concave <- FALSE
  use <- 2:(nrow(DF)-1)
  refpts <- m*DF$x[use] + b
  if (all(refpts > DF$y[use]) | all(refpts < DF$y[use])) concave <- TRUE
  if (!concave) stop("Your curve doesn't appear to be concave")
  
  # Calculate the orthogonal distances
  use <- 2:(nrow(DF)-1)
  elbowd <- distancePointLine(DF$x[use], DF$y[use], coef(fit)[2], coef(fit)[1])
  DF$dist <- c(NA, elbowd, NA) # first & last points don't have a distance
  
  if (plot) {
    edm <- which.max(DF$dist)
    plot(DF[,1:2], type = "b", xlab = "index", ylab = "y values",
         main = "Looking for the Elbow")
    segments(DF$x[1], DF$y[1],
             DF$x[nrow(DF)], DF$y[nrow(DF)], col = "red")
    points(DF$x[edm], DF$y[edm], cex = 1.5, col = "red")	
    points(DF$x[edm], DF$y[edm], pch = 20)
  }
  
  if (returnIndex) return(which.max(DF$dist))
  if (!returnIndex) return(DF)
  
} # end of findElbow

############################### Begin by importing input files and normalizing data for comparison of samples ################################

#############################
## Create Phyloseq object  ##
#############################
setwd("/Volumes/Renn_RNAseq_2017/microbiome_2017_stomach/R")
otufile <- read.table("phyloseq_OTUtable_SILVA.txt", header = TRUE, row.names=1)
## just import and then tell phyloseq that these are the files that it is looking for

#mapfile <- read.table("phyloseq_meta.txt", sep = "\t", header = TRUE, row.names=1)
mapfile <- read.csv("phyloseq_meta_small_SILVA.csv", sep = ",", header = TRUE, row.names=1)
# phyloseq needs rownames
meta_table<-sample_data(mapfile)
head(mapfile)
head(meta_table) ## they look the same

ep_tax<-read.table("phyloseq_taxa_SILVA.txt",header=TRUE, sep = "\t", row.names=1) 
# phyloseq needs rownames
ep_tax <- as.matrix(ep_tax)  # coerce to matrix
ep_tax <- tax_table(ep_tax)

ep_tree<-read.tree("filtered_seqs_min5_otus_SILVAannot_fasttree.tre")  # must be a rooted tree for phyloseq
ep_tree_rooted <-midpoint(ep_tree) # this is from phanghorn package that sets root to midpoint


#rarefy data!  ### optional related to normalization ## when different samples
####  have different read depths (these argumenents mattered more for 454)
otur <- rrarefy(t(otufile),6202)  # from vegan package (set number for dataset) (transpose)
# subsample to 5,987 reads, this is the minimum count for samples passing filters.
# make this number match your read count/filtered datasets

otu <- otu_table(t(otur), taxa_are_rows = TRUE)  # transpose back just to use commands

#lastly, turn all your inputs into phyloseq object
phyloseq_data <- phyloseq(otu, meta_table, ep_tax,ep_tree_rooted)

###########################
# For unrarified raw data #
###########################
otu <- otu_table(otufile, taxa_are_rows = TRUE)  # transpose back just to use commands
phyloseq_data_norarify <- phyloseq(otu, meta_table, ep_tax, ep_tree_rooted)


##################################################
#     Variance stabilizing normalization and     #
# differential abundance analysis through DESeq2 #
##################################################

# DESeq2
# *** MAKE SURE LOADED BEFORE PHYLOSEQ ***
otufile_pc <- read.table("phyloseq_OTUtable_SILVA.txt", header = TRUE, row.names=1)

otu_ds <- otu_table(otufile_pc, taxa_are_rows = TRUE) 
phyloseq_data_ds <- phyloseq(otu_ds, meta_table, ep_tax, ep_tree_rooted)

#Conversion to DESeq data
phyloseq_ds = phyloseq_to_deseq2(phyloseq_data_ds, ~ treat_c)
#Run negative bionomial Wald test to test for differential abundance of specific taxa
diagdds = DESeq(phyloseq_ds, test="Wald", fitType="parametric", betaPrior=TRUE)

# To extract results, choose the variable of interest and the pair of groups you want to test
# Cooks Cutoff dictates whether or not you want to remove outliers
res <- results( diagdds, contrast = c("treat_c", "B14", "R14"), cooksCutoff = FALSE )

# Filter for your desired alpha value, and extract relevant taxon names
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_data_ds)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#res[c("OTU_105", "OTU_25", "OTU_33", "OTU_9", "OTU_18", "OTU_5"),]
#res[c("OTU_20", "OTU_21", "OTU_212", "OTU_39", "OTU_2", "OTU_6", "OTU_10", "OTU_4"),]

# List results
sigtab

#### Plot differentially abundant taxa from pair of groups ####
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#####################################################
# Convert DESeq2 table to usable table for phyloseq #
# then substitute into the previous phyloseq object #
#####################################################

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)

# duplicate old phyloseq object for backup
phyloseq_data_ds2 <- phyloseq_data_ds

# Recommended to run the next command to floor data at 0 for input 
# into ordination or other downstream analyses

diagvst[diagvst < 0.0] <- 0.0

# Insert DESeq2 normalized data into phyloseq object
otu_table(phyloseq_data_ds2) <- otu_table(diagvst, taxa_are_rows = TRUE)

#optional output from DESeq2 normalized data
# write.table(diagvst, file = "DESeq2_table_floored.txt", sep = "\t")


####################################################################
# Subset the dataset for downstream analyses by metadata variables #
####################################################################
#new subsets, try subsetting based on other metadata variables for additional data exploration
#you probably just want to analyze the full dataset either raw or normalized, or use these subsets with no missing metadata
noNA <- subset_samples(phyloseq_data,check=="y")
noNA_norarify <- subset_samples(phyloseq_data_norarify,check=="y")
noNA_ds2 <- subset_samples(phyloseq_data_ds2,check=="y")

#define subset to test
microB_subset <- noNA
microB_subset <- noNA_norarify
microB_subset <- noNA_ds2
microB_subset <- phyloseq_data
microB_subset <- phyloseq_data_norarify
microB_subset <- phyloseq_data_ds2

################################################## BETA DIVERSITY #######################################################

#ordinations! try different methods and distance metrics
# primarily used "phyloseq_data_ds2" with values floored at 0 for this analysis, but used others for exploratory purposes.

# WARNING, do not use normalized OTU tables with negative values for anything but unifrac and wunifrac distances, otherwise
# it will likely crash R-studio since other distance methods cannot handle negatives.

ordu <- ordinate(microB_subset,"NMDS", "bray")
ordu2 <- ordinate(microB_subset,"NMDS", "unifrac")
ordu3 <- ordinate(microB_subset,"NMDS", "wunifrac")
ordu4 <- ordinate(microB_subset,"NMDS", "jaccard")
ordu5 <- ordinate(microB_subset,"PCoA", "bray")
ordu6 <- ordinate(microB_subset,"PCoA", "unifrac")
ordu7 <- ordinate(microB_subset,"PCoA", "wunifrac")
ordu8 <- ordinate(microB_subset,"PCoA", "jaccard")
ordu9 <- ordinate(microB_subset,"DPCoA")

plot_ordination(microB_subset, ordu7, type = "scree")
stressplot(ordu)

#for stock subsets
plot_ordination(microB_subset, ordu, color="treat_c", shape = "stock",title=paste("Bray-curtis NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu2, color="treat_c", shape = "treat_tank",title=paste("UniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu3, color="treat_c", shape = "treat_tank",title=paste("WuniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu4, color="treat_c", shape = "stock",title=paste("Jaccard NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu5, color="treat_c", shape = "stock",title=paste("Bray-curtis PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu6, color="treat_c", shape = "stock",title=paste("UniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu7, color="treat_c", shape = "stock",title=paste("WuniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu8, color="treat_c", shape = "stock",title=paste("Jaccard PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu9, color="treat_c", shape = "stock",title=paste("DPCoA")) + geom_point(size=4) + coord_fixed()

# Check plotting of different variables
#for treat_c subsets
plot_ordination(microB_subset, ordu, shape="treat_c", color = "stock") + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu2, shape="treat_c", color = "stock") + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu3, shape="treat_c", color = "stock") + geom_point(size=4) + coord_fixed()

plot_ordination(microB_subset, ordu, color="gen", shape = "stock",title=paste("Bray-curtis NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu2, color="gen", shape = "stock",title=paste("UniFrac NMDS")) + geom_point(size=4) + coord_fixed()

#Check effect of sequencing depth
#Do this after adding a read count column labeled "depth" to metadata table
plot_ordination(microB_subset, ordu, color="depth", shape = "stock",title=paste("Bray-curtis NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu2, color="depth", shape = "stock",title=paste("UniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu3, color="depth", shape = "stock",title=paste("WuniFrac NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu4, color="depth", shape = "stock",title=paste("Jaccard NMDS")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu5, color="depth", shape = "stock",title=paste("Bray-curtis PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu6, color="depth", shape = "stock",title=paste("UniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu7, color="depth", shape = "stock",title=paste("WuniFrac PCoA")) + geom_point(size=4) + coord_fixed()
plot_ordination(microB_subset, ordu8, color="depth", shape = "stock",title=paste("Jaccard PCoA")) + geom_point(size=4) + coord_fixed()


############################################
#beta diversity statistsical tests
#I primarily tested unifrac and wunifrac, so change distance method for each

#Use subset without any missing metadata, for example the variance-stabilizing normalization dataset
microB_subset <- noNA_ds2

metaD = as(sample_data(microB_subset), "data.frame")
phydist = distance(microB_subset, "wunifrac")

#ANOSIM, fit data to single variable groupings/distributions
#This tests whether two or more groups of samples are significantly different
#We prefer adonis over anosim for final analysis
anosim_res = anosim(phydist, metaD$treat_c, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$stock, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$GSI, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$BC, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$treat_tank, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$gen, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$student, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$mass_g, distance = "wunifrac")
anosim_res
anosim_res = anosim(phydist, metaD$length_mm, distance = "wunifrac")
anosim_res


#PERMANOVA by different variable groups
#This tests for differences in centroids in multidimensional matrix
adonis_res <- adonis2(formula = phydist ~ metaD$treat_c * metaD$stock * metaD$treat_tank * metaD$GSI * metaD$BC, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$treat_c, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$stock, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$GSI, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$BC, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$treat_tank, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$gen, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$mass_g, data = metaD) 
adonis_res
adonis_res <- adonis2(formula = phydist ~ metaD$length_mm, data = metaD) 
adonis_res


#This tests for differences in variances of clusters based on variables
beta <- betadisper(phydist, metaD$treat_c)
permutest(beta)

beta <- betadisper(phydist, metaD$stock)
permutest(beta)

beta <- betadisper(phydist, metaD$GSI)
permutest(beta)

beta <- betadisper(phydist, metaD$BC)
permutest(beta)

beta <- betadisper(phydist, metaD$treat_tank)
permutest(beta)

beta <- betadisper(phydist, metaD$gen)
permutest(beta)


################################################### alpha diversity ############################################################

#create rarefaction curve figure reflecting both observed number of OTUs and alpha diversity
rarefaction_curve_data <- calculate_rarefaction_curves(phyloseq_data_norarify, c('Observed', 'Shannon'), rep(c(1, 10, 100, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_se = ste(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(phyloseq_data)), by.x = 'Sample', by.y = 'row.names')

ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_se,
    ymax = Alpha_diversity_mean + Alpha_diversity_se,
    color = treat_c,
    group = Sample
  )
) + geom_line(
) + geom_point(size=1
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
)

#estimate species richness/diversity

#alpha diversity
#calculate measures per sample and join with GSI and BC values to check for correlations
#add more variables to your table to check for more correlations with diversity
microB_subset <- phyloseq_data_norarify
richness <- estimate_richness(microB_subset, split = TRUE, measures = NULL)
richness <- cbind(richness, treat_c = meta_table$treat_c, stock = meta_table$stock, GSI = meta_table$GSI, BC = meta_table$BC)
richness
# write table with diversity stats
write.table(richness, file = "richness.txt", sep = "\t")

#to actually test correlations between all fields in the diversity table, make sure your dataset has no "NA" values
#filter samples by subset function prior to correlation check
richness_mat <- as.matrix(cbind(richness, GSI = meta_table$GSI, BC = meta_table$BC))
richness_mat <- na.omit(richness_mat)
rcorr(richness_mat)

#plot alpha diversity stats on a per sample basis
plot_richness(microB_subset, x = "samples", color = NULL, shape = NULL,
              title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
              measures = NULL, sortby = NULL)

#make boxplot of alpha diversity stats for treatment group
p = plot_richness(microB_subset, x = "treat_c",
                  title = NULL, scales = "free_y", nrow = 1, shsi = NULL,
                  measures = c("Observed", "Shannon"), sortby = NULL)
p + geom_boxplot(data = p$microB_subset, aes(x = treat_c, y = value, color = NULL), alpha = 0.1)





###############################################
# Better way to estimate alpha diverstiy:     #
# Perform multiple iterations of rarification #
# before calulating diversity stats           #
###############################################

# Initialize matrices to store richness and evenness estimates
# Make sure subset is your unrarified dataset
microB_subset_norarify <- phyloseq_data_norarify
min_lib <- min(sample_sums(microB_subset_norarify))
nsamp = nsamples(microB_subset_norarify)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(microB_subset_norarify)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(microB_subset_norarify)

shannon <- matrix(nrow = nsamp, ncol = trials)
row.names(shannon) <- sample_names(microB_subset_norarify)

simpson <- matrix(nrow = nsamp, ncol = trials)
row.names(simpson) <- sample_names(microB_subset_norarify)

chao1 <- matrix(nrow = nsamp, ncol = trials)
row.names(chao1) <- sample_names(microB_subset_norarify)
chao1 <- rbind(chao1, chao1)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(microB_subset_norarify, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
  
  # Calculate shannon
  shan <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))
  shannon[ ,i] <- shan
  
  # Calculate simpson
  simp <- as.numeric(as.matrix(estimate_richness(r, measures = "Simpson")))
  simpson[ ,i] <- simp
  
  # Calculate Chao1
  chao <- as.numeric(as.matrix(estimate_richness(r, measures = "Chao1")))
  chao1[ ,i] <- chao
  
}

# Create a new dataframe to hold the means and standard error of richness estimates
SampleID <- row.names(richness)
mean_alpha <- apply(richness, 1, mean)
se <- apply(richness, 1, ste)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean_alpha, se, measure)

# Create a new dataframe to hold the means and standard error of evenness estimates
SampleID <- row.names(evenness)
mean_alpha <- apply(evenness, 1, mean)
se <- apply(evenness, 1, ste)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean_alpha, se, measure)

# Create a new dataframe to hold the means and standard error of diversity estimates
SampleID <- row.names(simpson)
mean_alpha <- apply(simpson, 1, mean)
se <- apply(simpson, 1, ste)
measure <- rep("Simpson", nsamp)
simp_stats <- data.frame(SampleID, mean_alpha, se, measure)

# Create a new dataframe to hold the means and standard error of diveristy estimates
SampleID <- row.names(shannon)
mean_alpha <- apply(shannon, 1, mean)
se <- apply(shannon, 1, ste)
measure <- rep("Shannon", nsamp)
shan_stats <- data.frame(SampleID, mean_alpha, se, measure)

# Create a new dataframe to hold the means and standard error of diversity estimates
SampleID <- row.names(chao1)
mean_alpha <- apply(chao1, 1, mean)
se <- apply(chao1, 1, ste)
measure <- rep("Chao1", nsamp)
chao_stats <- data.frame(SampleID, mean_alpha, se, measure)


#combine your stats and export tables

alpha <- rbind(rich_stats, even_stats, simp_stats, shan_stats)

s <- data.frame(sample_data(microB_subset))
alphadiv <- merge(alpha, s, by = "SampleID") 

write.csv(alphadiv, "alphadiv.csv")

#chao_stats contains both means and variances based on the initial calculation of chao1, therefore we need to export separate since it is not a 1:1 join with n samples.
#be cautious with how you treat means and se of the original chao1 se calculations.  They may be meaningless depending on your goals.
#chao1 mean stats are written first, followed by se calculations
write.csv(chao_stats, "chao1.csv")




#################################################
# Run ANOVA to test whether alpha diversity and #
# BC/GSI stats differ among treatment groups    #
#################################################

rich2 <- data.frame(merge(rich_stats, s, by = "SampleID"))
even2 <- data.frame(merge(even_stats, s, by = "SampleID"))
simp2 <- data.frame(merge(simp_stats, s, by = "SampleID"))
shan2 <- data.frame(merge(shan_stats, s, by = "SampleID"))

fit <- aov(mean_alpha ~ treat_c, data=rich2)
summary(fit)

fit <- aov(mean_alpha ~ treat_c, data=even2)
summary(fit)

fit <- aov(mean_alpha ~ treat_c, data=simp2)
summary(fit)

fit <- aov(mean_alpha ~ treat_c, data=shan2)
summary(fit)

fit <- aov(GSI ~ treat_c, data=shan2)
summary(fit)

fit <- aov(BC ~ treat_c, data=shan2)
summary(fit)

# Run pairwise tests to see specifically which groups differ
fit <- pairwise.t.test(rich2$mean_alpha, rich2$treat_c, p.adj = "bonf")
fit

fit <- pairwise.t.test(even2$mean_alpha, even2$treat_c, p.adj = "bonf")
fit

fit <- pairwise.t.test(simp2$mean_alpha, simp2$treat_c, p.adj = "bonf")
fit

fit <- pairwise.t.test(shan2$mean_alpha, shan2$treat_c, p.adj = "bonf")
fit

fit <- pairwise.t.test(shan2$GSI, shan2$treat_c, p.adj = "bonf")
fit

fit <- pairwise.t.test(shan2$BC, shan2$treat_c, p.adj = "bonf")
fit


#test assumptions of equal variance for ANOVAs
shapiro.test(shan2$mean_alpha)
shapiro.test(simp2$mean_alpha)
shapiro.test(even2$mean_alpha)
shapiro.test(rich2$mean_alpha)


#comparison and plotting of GSI and body condition among groups
fullmeta <- read.csv("phyloseq_meta_small_SILVA.csv", sep = ",", header = TRUE, row.names=1)

#if you want to run GSI and BC comparisons on Wild Stock only, run the following command
#fullmeta <- subset(fullmeta, stock == "WS")

#find mean, se for BC
AVG<-aggregate(BC ~ treat_c , fullmeta, mean, na.rm = T) 
SE<-aggregate(BC ~ treat_c, fullmeta, ste) 
BC <- data.frame(merge(AVG, SE, by = "treat_c"))

#find mean, se for GSI
AVG<-aggregate(GSI ~ treat_c , fullmeta, mean, na.rm = T) 
SE<-aggregate(GSI ~ treat_c , fullmeta, ste) 
GSI <- data.frame(merge(AVG, SE, by = "treat_c"))

ggplot(BC, aes(x=factor(treat_c), y=BC.x)) + 
  stat_summary(fun.y="mean", geom="bar") + 
  geom_errorbar(aes(ymin=BC$BC.x - BC$BC.y, ymax=BC$BC.x + BC$BC.y),
                width=0.2,                    # Width of the error bars
                position=position_dodge(.9))

ggplot(GSI, aes(x=factor(treat_c), y=GSI.x)) + 
  stat_summary(fun.y="mean", geom="bar") + 
  geom_errorbar(aes(ymin=GSI$GSI.x - GSI$GSI.y, ymax=GSI$GSI.x + GSI$GSI.y),
                width=0.2,                    # Width of the error bars
                position=position_dodge(.9))

fit <- pairwise.t.test(fullmeta$BC, fullmeta$treat_c, p.adj = "bonf")
fit

fit <- pairwise.t.test(fullmeta$GSI, fullmeta$treat_c, p.adj = "bonf")
fit


################################################### PhyloFactor ##############################################################

##########################################################################
# In addition to the DESeq2 analysis to test for differential abundance, #
# this analysis can detect either individual OTUs or entire clades that  #
# explain differences among groups                                       #
##########################################################################

#Import abundance table
#OTUs in table must be sorted to reflect the order of OTUs in the FastTree phylogeny
#This table is also sample-sorted by treatment for ease of making figures
otufile <- read.table("pf_OTUtable_SILVA.txt", header = TRUE, row.names=1)
otu <- data.matrix(otufile)
#This table is also sample-sorted by stock for ease of making figures
#otufile <- read.table("phyloseq_otu_sort_stock.txt", header = TRUE, row.names=1)

# rarefy data!  ### recommended normalization for when different samples
# have different read depths. Keep checking Phylofactor documentation to see
# if this changes, since rarification is not always ideal, but often good enough
otur <- rrarefy(t(otufile),6202)  # from vegan package (set number for dataset) (transpose)
otu <- otu_table(t(otur), taxa_are_rows = TRUE)  # transpose back just to use commands

#Import unrooted 16S tree for OTUs created by FastTree
pf_tree<-read.tree("filtered_seqs_min5_otus_SILVAannot_fasttree.tre")  # must be an unrooted tree for phylofactor


#Import metadata, specifically with the stock and treatment categorizations of samples
mapfile <- read.csv("phyloseq_meta_small_SILVA.csv", sep = ",", header = TRUE, row.names=1)
#optional stock sorted metadata for ease of making figures below
#mapfile <- read.csv("phyloseq_meta_stock.csv", sep = ",", header = TRUE, row.names=1)
#optional for WS only
#mapfile <- read.csv("phyloseq_meta_WS.csv", sep = ",", header = TRUE, row.names=1)
# phyloseq needs rownames
meta_table<-sample_data(mapfile)
head(mapfile)
head(meta_table) ## they look the same
#pull out stock and/or treatment variables to perform phylofactorization
meta_list <- meta_table$treat_c
meta_list2 <- meta_table$stock

#Import taxa designations for OTUs
pf_tax<-data.frame(read.table("pf_taxa_SILVA.txt",header=TRUE, sep = ","))


###OPTIONAL, filter out OTUs that are in more than 5 samples###
# Choose your own number based on portion of samples you're analyzizng
# This set of commands also prunes the phylogeny, you can prune after
# running phylofactor as well to make heatmaps more readable with common
# taxa only
common.otus <- which(rowSums(otu>0)>5)
otu <- otu[common.otus,]
pf_tree <- ape::drop.tip(pf_tree,setdiff(pf_tree$tip.label,rownames(otu))) 
pf_tax <- pf_tax[match(pf_tree$tip.label,pf_tax[,1]),]
dim(otu)


#Perform Phylofactorization
#you can change number of factors to look at additional phylogenetic dimensions within dataset
#you can also change meta_list based on which variables you want to consider in PF analysis

#first, set order of otu to order of tip labels in tree for easy downstream processing
otu <- otu[pf_tree$tip.label,]

#for first pass, set nfactors to a large number, ex. >20
#this will allow you to plot a scree-plot of explained variation
#from each factor, then look for elbow

#PF <- PhyloFactor(otu,pf_tree,meta_list,nfactors=20)

#used 7 factors for final analysis given the detection of the curve elbow using ExpVar
PF <- PhyloFactor(otu,pf_tree,meta_list,nfactors=7)

#for stock
#PF <- PhyloFactor(otu,pf_tree,meta_list2,nfactors=20)
#Using Variance Stabilized transformations
#PF <- PhyloFactor(diagvst,pf_tree,meta_list,nfactors=20)


#Look at summary of PF results, specifically which phylogenetic clades and OTUs explain variable of interest 
#This table has many of the main results from PF like ExpVar, F, and Pr(>F)
PF$factors

#names(PF)

#######################################
# This next set of commands will help #
# you know when to stop factorization #
#######################################

# Find the elbow of the scree plot for ExpVar
# To do so, you will preform loess smoothing of the ExpVar plot
# then look for  the peak in the change of slope

# For this exploratory analysis you need to use a larger number 
# of factors than you expect to be biologically significant to 
# identify a somewhat empirical cutoff (ex: n=20)

elb <- PF$factors
elb$Factors <- 1:nrow(elb)
y <- elb$ExpVar
x <- elb$Factors
par(mfrow=c(1,2))
plot(x,y,type="l")
lo <- loess(y~x)
#xl <- seq(min(x),max(x), (max(x) - min(x))/19) # set the denominator to (nfactors - 1)
out = predict(lo,x)
lines(x, out, col='red', lwd=2)

findElbow(out, plot = TRUE, returnIndex = TRUE)
par(mfrow=c(1,1))

#This function finds the elbow of the curve. Take this value and rerun phylofactor 
#with this as the number of factors to report. This method is not perfect but the 
#best way at getting to the number of factors that may be biologically interesting. 
#The author of phylofactor recognizes that one of the main limitations of his program 
#is that there's no systematic way to know when to stop factorization, and this is 
#our best attempt at solving the problem. One remaining issue with this method is that
#it requires a non-noisy concave curve and our loess smoothing may over-smooth the 
#data to the point where it doesn't quite match the original biological data, so 
#proceed with caution. Still better than nothing though.

#####################################################################################
# creates full data summary for factor of interest, and then allows you to look at  #
# relative abundance ratios of the OTU groups defined by each factor on a sample to #
# sample basis                                                                      #
#####################################################################################

# pick factor to examine more in-depth, it will be labeled as group 1 in results
smry <- pf.summary(PF,pf_tax,factor=1)
pf.tidy(smry)


#Allows you to visualize the abundance ratios between OTU groups for the summarized 
#factor chosen above by "pf.summary" function in each of your sample variable categories

#names(smry)
par(mfrow=c(1,2))
plot(meta_list,smry$ilr,main='Isometric log-ratio',ylab='ILR balance')
plot(meta_list,smry$MeanRatio,main='Ratio of Geometric Means',ylab='Group1/Group2')
#for stock
#names(smry)
#par(mfrow=c(1,2))
#plot(meta_list2,smry$ilr,main='Isometric log-ratio',ylab='ILR balance')
#plot(meta_list2,smry$MeanRatio,main='Ratio of Geometric Means',ylab='Group1/Group2')


#Look at the OTU group from the summarized factors 
#find taxonomy
smry$TaxaSplit %>%
  lapply(.,FUN=function(x) unique(x$TaxaIDs)) 
#find OTU names
Group1.otus <- smry$group1$IDs$otuIDs %>% as.list %>% sapply(.,toString)
Group1.otus

##########################################################
# Create heatmap of OTU and clade abundance vs predicted #
##########################################################
# To look at phylogeny only:
#par(mfrow=c(1,1))
#plot.phylo(pf_tree,use.edge.length=FALSE,main='Community Phylogeny')

#create single heatmap of log relative abundance of OTUs based on phylogenetic relationships.
#it takes a while to run with all OTUs, so pruning tree is recommended for time and readability 
colors<-colorRampPalette(colors=c("white", "blue"))(100)
clr <- function(Matrix) apply(Matrix,MARGIN=2,FUN=function(x) log(x)-mean(log(x)))
par(mfrow=c(1,1))
# heatmap for IRL values (log abundance relative to other OTUs within sample)
phylo.heatmap(pf_tree,colors=colors,clr(PF$Data))
# heatmap for raw abundance, but tends to be skewed from high adundance taxa, not very readable
#phylo.heatmap(pf_tree,colors=colors,clr(PF$Data))

# Create multiple heatmaps, showing raw log relative abundance and averages for factors and groups
# Predictions are the average of factors within groups
# to make smaller heatmap with unfiltered OTU dataset, run PF for full dataset then prune phylogeny for input into phylo.heatmap
clr <- function(Matrix) apply(Matrix,MARGIN=2,FUN=function(x) log(x)-mean(log(x)))
colors<-colorRampPalette(colors=c("white", "blue","blue4", "darkblue", "navy"))(100)
#colors<-colorRampPalette(colors=c("white", "red1", "red", "darkred"))(100)
#colors<-colorRampPalette(c("white","green","green4","violet","purple"))(100)
prediction <- pf.predict(PF)
dim(prediction)
prediction <- pf.predict(PF,factors=7)
par(mfrow=c(2,1))
phylo.heatmap(pf_tree,clr(PF$Data),colors=colors, fsize=0.5)
phylo.heatmap(pf_tree,clr(prediction),colors=colors, fsize=0.5)
#phylo.heatmap(pf_tree,clr(PF$Data), fsize=0.5)
#phylo.heatmap(pf_tree,clr(prediction),Colv=FALSE, fsize=0.5)

################################################ Other useful plots ############################################################

#Bar plots for microbiome composition per individal, per treatment, and per treatment grouped by phylum respectively
#Change taxonomic level or variable as desired for exploration of data
plot_bar(phyloseq_data, fill="Phylum")
plot_bar(phyloseq_data, "treat_c", fill="Phylum")
plot_bar(phyloseq_data, "treat_c", fill="treat_c", facet_grid=~Phylum)

plot_bar(phyloseq_data_ds2, fill="Phylum")
plot_bar(phyloseq_data_ds2, "treat_c", fill="Phylum")
plot_bar(phyloseq_data_ds2, "treat_c", fill="treat_c", facet_grid=~Phylum)

plot_bar(phyloseq_data_norarify, fill="Phylum")
plot_bar(phyloseq_data_norarify, "treat_c", fill="Phylum")
plot_bar(phyloseq_data_norarify, "treat_c", facet_grid=~Phylum)

#Transform to make proportion within treatments from DESeq2 normalized data
phyloseq_data_ds3 = merge_samples(phyloseq_data_ds2, "treat_c")
sample_data(phyloseq_data_ds3)$treat_c <- levels(sample_data(phyloseq_data_ds3)$treat_c)
phyloseq_data_ds3 = transform_sample_counts(phyloseq_data_ds3, function(x) 100 * x/sum(x))
#remove lines from barplot
p = plot_bar(phyloseq_data_ds3, fill="Phylum")
p + geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position="stack")


#Transform to make proportion within treatments from rarified data
phyloseq_data_3 = merge_samples(phyloseq_data, "treat_c")
sample_data(phyloseq_data_3)$treat_c <- levels(sample_data(phyloseq_data_3)$treat_c)
phyloseq_data_3 = transform_sample_counts(phyloseq_data_3, function(x) 100 * x/sum(x))
#remove lines from barplot
p = plot_bar(phyloseq_data_3, fill="Phylum")
p + geom_bar(aes(color=Phylum, fill=Phylum), stat = "identity", position="stack")



