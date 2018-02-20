# # installing/loading the package:
# if(!require(installr)) {
#   install.packages("installr"); require(installr)} #load / install+load installr
# 
# # using the package:
# updateR() # this will start the updating process of your R installation.  It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 

packages <- c("ggplot2", "ranger", "svglite", "ggthemes", "phyloseq", "vegan", "Rtsne", "dbscan", "mefa", "caret", "data.table")

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

phyla_levels_all <-
  c("AC1",
    "Acetothermia",
    "Acidobacteria",
    "Actinobacteria",
    "Aminicenantes",
    "Aquificae",
    "Armatimonadetes",
    "Atribacteria",
    "Bacteroidetes",
    "Candidatus Berkelbacteria",
    "BJ-169",
    "BRC1",
    "Caldiserica",
    "Chlamydiae",
    "Chlorobi",
    "Chloroflexi",
    "Cloacimonetes",
    "CPR2",
    "Cyanobacteria",
    "Deinococcus-Thermus",
    "Elusimicrobia",
    "FBP",
    "FCPU426",
    "Fibrobacteres",
    "Firmicutes",
    "GAL15",
    "Gemmatimonadetes",
    "Gracilibacteria",
    "Hydrogenedentes",
    "Ignavibacteriae",
    "Latescibacteria",
    "LCP-89",
    "Lentisphaerae",
    "Microgenomates",
    "Nitrospinae",
    "Nitrospirae",
    "Omnitrophica",
    "Parcubacteria",
    "PAUC34f",
    "Peregrinibacteria",
    "Planctomycetes",
    "Proteobacteria",
    "Alphaproteobacteria",
    "Betaproteobacteria",
    "Gammaproteobacteria",
    "Deltaproteobacteria",
    "RBG-1 (Zixibacteria)",
    "Candidatus Saccharibacteria",
    "Saccharibacteria",
    "Spirochaetae",
    "SR1 (Absconditabacteria)",
    "Synergistetes",
    "Tectomicrobia",
    "TM6 (Dependentiae)",
    "Verrucomicrobia",
    "WS2",
    "WS6",
    "WWE3",
    "Bathyarchaeota",
    "Crenarchaeota",
    "Euryarchaeota",
    "Miscellaneous Euryarchaeotic Group(MEG)",
    "Parvarchaeota",
    "Thaumarchaeota",
    "Woesearchaeota (DHVEG-6)",
    "Other")

#color palettes from http://tools.medialab.sciences-po.fr/iwanthue/
darkcols <- c(
  "#ad6a43",
  "#5f2900",
  "#ffc686",
  "#2dc1ff",
  "#b9e08f",
  "#802e1d",
  "#005516",
  "#e0d775",
  "#d69ef4",
  "#015aad",
  "#954e00",
  "#00ca88",
  "#6c2d92",
  "#367900",
  "#fd7353",
  "#fb9b38",
  "#b4d043",
  "#68d86a",
  "#3174ec",
  "#c95d0f",
  "#999999"
)

#read in the taxonomic metaphlan2/matam and functional humann2 data
metaphlan2_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/metaphlan2_merged.txt", sep = "\t")
matam_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/matam_contingency_bugged.txt", sep = "\t")
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

metaphlan2_phyla <- metaphlan2_raw[grep("p__", metaphlan2_raw$ID),]
metaphlan2_phyla <- metaphlan2_phyla[-grep("o__", metaphlan2_phyla$ID),]
metaphlan2_prot <- metaphlan2_phyla[grep("proteobacteria", metaphlan2_phyla$ID),]
metaphlan2_phyla <- metaphlan2_phyla[-grep("c__|Proteobacteria", metaphlan2_phyla$ID),]
metaphlan2_phyla <- rbind(metaphlan2_phyla, metaphlan2_prot)
metaphlan2_phyla_rownames <- metaphlan2_phyla$ID
metaphlan2_phyla[,"ID"] <- list(NULL)
metaphlan2_phyla <- data.frame(apply(metaphlan2_phyla, 2, function(x) as.numeric(as.character(x))))
rownames(metaphlan2_phyla) <- metaphlan2_phyla_rownames


humann2_data <- humann2_raw[-grep("\\|", humann2_raw$X..Pathway),]
humann2_rownames <- humann2_data$X..Pathway
humann2_data[,"X..Pathway"] <- list(NULL)
humann2_data <- data.frame(apply(humann2_data, 2, function(x) as.numeric(as.character(x))))
colnames(humann2_data) <- sub("_Abundance", "", colnames(humann2_data))
rownames(humann2_data) <- humann2_rownames

#handle also the matam data
matam_taxonomy <- as.character(matam_raw$Taxonomy.Samples)
matam_data <- matam_raw[-1]
matam_data <- data.frame(apply(matam_data, 2, as.numeric))
colnames(matam_data) <- sub("X", "", colnames(matam_data))
colnames(matam_data) <- paste0("Lane_", sub("_[0-9]+", "", colnames(matam_data)), "_DNA_", sub("00[0-3]_", "", colnames(matam_data)))
rownames(matam_data) <- paste0("OTU_", 1:nrow(matam_data))
matam_data[is.na(matam_data)] <- 0
matam_tax <- strsplit(matam_taxonomy,";")
matam_tax <- as.matrix(rbind.fill(lapply(matam_tax,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})))
#remove subclasses and suborders from actinobacteria
actinobacteria_tax <- matam_tax[matam_tax[,2] %in% "Actinobacteria",c(-4, -6)]
matam_tax <- matam_tax[,-7:-8]
matam_tax[matam_tax[,2] %in% "Actinobacteria",] <- actinobacteria_tax
colnames(matam_tax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rownames(matam_tax) <- paste0("OTU_", 1:nrow(matam_tax))

#gather sample data
sample_mapping_data <- import_qiime_sample_data("D:/VirtualBox/VirtualBox Share/hazen_metagenome/sample_data.csv")
microprobe_data <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/microprobes.csv",sep="\t")
microprobe_data <- transform(microprobe_data, bin = cut(Depth.mm, breaks = seq(0,35,by=5), include.lowest = T))
microprobe_data_means <- ddply(microprobe_data, .(Core,bin), numcolwise(median))[-c(1:3)]
microprobe_data_means <- cbind(microprobe_data_means, bin=rep(unique(microprobe_data$bin),2))
microprobe_data_sds <- ddply(microprobe_data, .(Core,bin), numcolwise(sd))[,-c(1:3)]
microprobe_data_means <- microprobe_data_means[c(1,3,5,8,11,14),-length(microprobe_data_means)]
microprobe_data_means <- mefa:::rep.data.frame(microprobe_data_means, each = 3)
sample_mapping_data <- cbind(sample_mapping_data[,-1], microprobe_data_means)
sample_mapping_data$Sample <- factor(sample_mapping_data$Sample)

#construct phyloseq objects
metaphlan2 <- phyloseq(otu_table(metaphlan2_data, taxa_are_rows = T), sample_data(sample_mapping_data))
humann2 <- phyloseq(otu_table(humann2_data, taxa_are_rows = T), sample_data(sample_mapping_data))
matam <- phyloseq(otu_table(matam_data, taxa_are_rows = T), tax_table(matam_tax), sample_data(sample_mapping_data))

#remove chloroplasts and mitochondria, and correct some phyla names
matam <- subset_taxa(matam, !(Domain %in% "Eukaryota"))
matam <- subset_taxa(matam, !(Class %in% "Chloroplast"))
matam <- subset_taxa(matam, !(Family %in% "Mitochondria"))
tax_table(matam)[tax_table(matam)[,2] %in% "Cyanobacteria/Chloroplast",2] <- "Cyanobacteria"
tax_table(matam)[tax_table(matam)[,2] %in% "Woesearchaeota",2] <- "Woesearchaeota (DHVEG-6)"

#metaphlan2 bars
bars_metaphlan2 <- data.frame(t(metaphlan2_phyla))
colnames(bars_metaphlan2) <- sub(".*?_{2,3}(.+)$", "\\1", colnames(bars_metaphlan2))
bars_metaphlan2$sample <- sub(".*?([0-9]+)$", "\\1", rownames(bars_metaphlan2))
bars_metaphlan2 <- aggregate(.~sample, data = bars_metaphlan2, mean)
bars_metaphlan2 <- bars_metaphlan2[order(as.numeric(bars_metaphlan2$sample)),]
bars_metaphlan2 <- bars_metaphlan2[-1]
bars_metaphlan2 <- data.frame(apply(bars_metaphlan2, 2, function(x) as.numeric(as.character(x))))
sample_order <- c("Deep Hole 0.5 cm", "Deep Hole 1.5 cm", "Deep Hole 2.5 cm", "Snowgoose Bay 0.5 cm", "Snowgoose Bay 2.0 cm", "Snowgoose Bay 3.5 cm")
rownames(bars_metaphlan2) <- sample_order
bars_metaphlan2 <- melt(as.matrix(bars_metaphlan2), varnames=c("Sample", "Phylum"), value.name="Abundance")

#plot phyla
ggplot(bars_metaphlan2, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_metaphlan2$Sample[nchar(as.character(bars_metaphlan2$Sample))!=1]) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=darkcols, name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("Taxonomy")

#phyla bars for the matam data
phyla_matam <- matam
if(any(tax_table(phyla_matam)[,2] %in% "Proteobacteria")) {
  replacement_proteobacteria <- tax_table(phyla_matam)[tax_table(phyla_matam)[,2] %in% "Proteobacteria",3]
  replacement_proteobacteria[replacement_proteobacteria %in% NA] <- "Proteobacteria"
  tax_table(phyla_matam)[tax_table(phyla_matam)[,2] %in% "Proteobacteria",2] <- replacement_proteobacteria
}
#merge phyla_matam that have less than 1% abundance (or unknown Phylum) in the data set to "other"
phyla_matam <- tax_glom(phyla_matam, "Phylum")
least_abundant_phyla_matam <- names(taxa_sums(phyla_matam))[taxa_sums(phyla_matam)/sum(taxa_sums(phyla_matam))*100 < 1]
phyla_matam <- merge_taxa(phyla_matam, least_abundant_phyla_matam, 1)
tax_table(phyla_matam)[tax_table(phyla_matam)[,2] %in% NA,2] <- "Other"
phyla_matam <- merge_samples(phyla_matam, group="Sample", fun=mean)
phyla_matam <- transform_sample_counts(phyla_matam, function(x) 100 * x/sum(x))
bars_matam <- data.frame(otu_table(phyla_matam))
darkcols_matam <- darkcols[(length(darkcols)+1-ncol(bars_matam)):length(darkcols)]
colnames(bars_matam) <- c(tax_table(phyla_matam)[,2])
phyla_matam_levels <- phyla_levels_all[phyla_levels_all %in% colnames(bars_matam)]
bars_matam <- bars_matam[,phyla_matam_levels]
rownames(bars_matam) <- sample_order
bars_matam <- melt(as.matrix(bars_matam), varnames=c("Sample", "Phylum"), value.name="Abundance")
bars_matam$Phylum <- factor(bars_matam$Phylum, levels = rev(phyla_matam_levels))

#plot phyla
ggplot(bars_matam, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_matam$Sample[nchar(as.character(bars_matam$Sample))!=1]) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(darkcols_matam), name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("Taxonomy")

#calculate bray-curtis dissimilarities, and present them as NMDS ordinations
metaphlan2_dist <- vegdist(t(metaphlan2_data), distance = "bray")
set.seed(42)
ordu_metaphlan2 <- ordinate(metaphlan2, "NMDS", distance = metaphlan2_dist)
ef_metaphlan2 <- envfit(ordu_metaphlan2,sample_data(metaphlan2),permu=10000)
NMDS_data_metaphlan2 <- data.frame(sample_data(metaphlan2))
NMDS_data_metaphlan2$NMDS1 <- ordu_metaphlan2$points[ ,1]
NMDS_data_metaphlan2$NMDS2 <- ordu_metaphlan2$points[ ,2]

humann2_dist <- vegdist(t(humann2_data), distance = "bray")
set.seed(42)
ordu_humann2 <- ordinate(humann2, "NMDS", distance = humann2_dist)
ef_humann2 <- envfit(ordu_humann2,sample_data(humann2),permu=10000)
NMDS_data_humann2 <- data.frame(sample_data(humann2))
NMDS_data_humann2$NMDS1 <- ordu_humann2$points[ ,1]
NMDS_data_humann2$NMDS2 <- ordu_humann2$points[ ,2]

ggplot(data = NMDS_data_metaphlan2, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=Site), size = 5) + labs(color = "Sample")

ggplot(data = NMDS_data_humann2, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=Site), size = 5) + labs(color = "Sample")

#implement caret