# # installing/loading the package:
# if(!require(installr)) {
#   install.packages("installr"); require(installr)} #load / install+load installr
# 
# # using the package:
# updateR() # this will start the updating process of your R installation.  It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 

packages <- c("ggplot2", "ranger", "svglite", "ggthemes", "phyloseq", "vegan", "Rtsne", "dbscan")

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, dependencies = T, suppressUpdates = T)
  sapply(pkg, require, character.only = TRUE)
}

ipak(packages)


#read in functions from miseq (https://github.com/michberr/MicrobeMiseq/blob/master/R/miseqR.R)
source("C:/Users/Matti/Dropbox/Matti-Alex/R/miseqR.R")


theme_set(theme_tufte(base_family = "sans", base_size = 18) + theme(panel.border = element_rect(colour = "black", fill = NA), 
                                                                    axis.text = element_text(colour = "black", size = 18)))

#read in the taxonomic metaphlan2 and functional humann2 data
metaphlan2_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/metaphlan2_merged.txt", sep = "\t")
humann2_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/humann2_merged_norm.tsv", sep = "\t")

#limit metaphlan2 data to species (tips of the tree) or phyla and humann2 data to just the functions (remove tax information)
metaphlan2_data <- metaphlan2_raw[grep("s__", metaphlan2_raw$ID),]
metaphlan2_data <- metaphlan2_data[-grep("t__", metaphlan2_data$ID),]
metaphlan2_data$ID <- as.character(metaphlan2_data$ID)
metaphlan2_data$ID <- sub("^(.+s__).*?", "", metaphlan2_data$ID)
metaphlan2_rownames <- metaphlan2_data$ID
metaphlan2_data[,"ID"] <- list(NULL)
metaphlan2_data <- data.frame(apply(metaphlan2_data, 2, function(x) as.numeric(as.character(x))))
#summarize taxa classified at lower levels as a single group
metaphlan2_data <- rbind(metaphlan2_data, 100-colSums(metaphlan2_data))
rownames(metaphlan2_data) <- c(metaphlan2_rownames, "other")

metaphlan2_phyla <- metaphlan2_raw[grep("s__", metaphlan2_raw$ID),]


humann2_data <- humann2_raw[-grep("\\|", humann2_raw$X..Pathway),]
humann2_rownames <- humann2_data$X..Pathway
humann2_data[,"X..Pathway"] <- list(NULL)
humann2_data <- data.frame(apply(humann2_data, 2, function(x) as.numeric(as.character(x))))
colnames(humann2_data) <- sub("_Abundance", "", colnames(humann2_data))
rownames(humann2_data) <- humann2_rownames
  
#construct phyloseq objects
sample_mapping_data <- import_qiime_sample_data("D:/VirtualBox/VirtualBox Share/hazen_metagenome/sample_data.csv")
metaphlan2 <- phyloseq(otu_table(metaphlan2_data, taxa_are_rows = T), sample_data(sample_mapping_data))
humann2 <- phyloseq(otu_table(humann2_data, taxa_are_rows = T), sample_data(sample_mapping_data))

#calculate bray-curtis dissimilarities, and present them as NMDS ordinations
metaphlan2_dist <- vegdist(t(metaphlan2_data), distance = "bray")
set.seed(42)
ordu_metaphlan2 <- ordinate(metaphlan2, "NMDS", distance = metaphlan2_dist)
NMDS_data_metaphlan2 <- data.frame(sample_data(metaphlan2))
NMDS_data_metaphlan2$NMDS1 <- ordu_metaphlan2$points[ ,1]
NMDS_data_metaphlan2$NMDS2 <- ordu_metaphlan2$points[ ,2]

humann2_dist <- vegdist(t(humann2_data), distance = "bray")
set.seed(42)
ordu_humann2 <- ordinate(humann2, "NMDS", distance = humann2_dist)
NMDS_data_humann2 <- data.frame(sample_data(humann2))
NMDS_data_humann2$NMDS1 <- ordu_humann2$points[ ,1]
NMDS_data_humann2$NMDS2 <- ordu_humann2$points[ ,2]

ggplot(data = NMDS_data_metaphlan2, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=Site), size = 5) + labs(color = "Sample")

ggplot(data = NMDS_data_humann2, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=Site), size = 5) + labs(color = "Sample")