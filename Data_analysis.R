# # installing/loading the package:
# if(!require(installr)) {
#   install.packages("installr"); require(installr)} #load / install+load installr
# 
# # using the package:
# updateR() # this will start the updating process of your R installation.  It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 

packages <- c("ggplot2", "ranger", "svglite", "ggthemes", "phyloseq", "vegan", "Rtsne", "dbscan", "mefa", "caret", "data.table", "plyr", "protr", "data.table", "MLPUGS", "compiler")

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, dependencies = T, suppressUpdates = T)
  sapply(pkg, require, character.only = TRUE)
}

ipak(packages)

enableJIT(3)


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
matam_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/matam_contingency_samples_combined.txt", quote = "", sep = "\t")
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
#colnames(matam_data) <- paste0("Lane_", sub("_[0-9]+", "", colnames(matam_data)), "_DNA_", sub("00[0-3]_", "", colnames(matam_data)))
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
matam_mapping_data <- merge_samples(sample_data(sample_mapping_data), group="Sample")
matam <- phyloseq(otu_table(matam_data, taxa_are_rows = T), tax_table(matam_tax), sample_data(matam_mapping_data))

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
#phyla_matam <- merge_samples(phyla_matam, group="Sample", fun=mean)
phyla_matam <- transform_sample_counts(phyla_matam, function(x) 100 * x/sum(x))
bars_matam <- data.frame(t(otu_table(phyla_matam)))
darkcols_matam <- darkcols[(length(darkcols)+1-ncol(bars_matam)):length(darkcols)]
colnames(bars_matam) <- c(tax_table(phyla_matam)[,2])
phyla_matam_levels <- phyla_levels_all[phyla_levels_all %in% colnames(bars_matam)]
bars_matam <- bars_matam[,phyla_matam_levels]
rownames(bars_matam) <- sample_order
bars_matam <- melt(as.matrix(bars_matam), varnames=c("Sample", "Phylum"), value.name="Abundance")
bars_matam$Phylum <- factor(bars_matam$Phylum, levels = rev(phyla_matam_levels))

#write out phyla abundances
phyla_percentages <- data.frame(otu_table(phyla_matam))
rownames(phyla_percentages) <- c(tax_table(phyla_matam)[,2])
phyla_percentages$rowmeans <- rowMeans(phyla_percentages)
write.csv(phyla_percentages, "D:/VirtualBox/VirtualBox Share/hazen_metagenome/phyla_abundances.csv", quote=F)

#plot phyla
ggplot(bars_matam, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_matam$Sample[nchar(as.character(bars_matam$Sample))!=1]) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(darkcols_matam), name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("Taxonomy")

#write quality controlled biom files for FAPROTAX
tax_table <- data.frame(tax_table(matam))
tax0 <- do.call("paste", c(list(rownames(tax_table)), tax_table, sep="; "))
tax0 <- sub("^.+(Bacteria.*?); NA.*", "\\1", tax0)
tax0 <- sub("^.+(Archaea.*?); NA.*", "\\1", tax0)
tax0 <- sub("^OTU_[0-9]*; ", "", tax0)
otu_table <- data.frame(otu_table(matam))
colnames(otu_table) <- sample_names(matam)
otu_table$taxonomy <- tax0
otu_table <- cbind(OTU=rownames(otu_table), otu_table)
write.table(otu_table, "D:/VirtualBox/VirtualBox Share/hazen_metagenome/matam_tax_table_qcd.txt", quote=F, sep="\t", row.names = F)

#the following commands ran at CAC cluster
# biom convert -i ~/matti/Hazen_metagenome/matam_tax_table_qcd.txt -o ~/matti/Hazen_metagenome/matam_tax_table_qcd.biom --table-type="OTU table" --to-json
# python ~/bin/FAPROTAX_1.0/collapse_table.py -i ~/matti/Hazen_metagenome/matam_tax_table_qcd.biom -o ~/matti/Hazen_metagenome/matam_func_table.biom -g ~/bin/FAPROTAX_1.0/FAPROTAX_Hazen.txt --collapse_by_metadata 'taxonomy' --group_leftovers_as 'other' --out_group_overlaps ~/matti/Hazen_metagenome/matam_func_table_overlaps.txt --output_format_group_overlaps classical  -l ~/matti/Hazen_metagenome/matam_FAPROTAX.log --disable_group_set_operations --out_groups2records_table ~/matti/Hazen_metagenome/matam_func_table_groups.txt -v --force
# biom convert -i ~/matti/Hazen_metagenome/matam_func_table.biom -o ~/matti/Hazen_metagenome/matam_func_table.txt --to-tsv

func <- read.csv("D:/VirtualBox/VirtualBox Share/Hazen_metagenome/matam_func_table.txt",sep="\t",skip=1,row.names=1)
colnames(func) <- sub("^.", "", colnames(func))

func <- phyloseq(otu_table(func,taxa_are_rows = T), sample_data(matam))
func <- prune_taxa(taxa_sums(func)>0, func)
func_pruned <- prune_taxa(!(taxa_names(func) %in% "other"),func)
rownames(otu_table(func_pruned)) <- paste(toupper(substr(rownames(otu_table(func_pruned)), 1, 1)), substr(rownames(otu_table(func_pruned)), 2, nchar(rownames(otu_table(func_pruned)))), sep="")
rownames(otu_table(func_pruned)) <- gsub("_", " ", rownames(otu_table(func_pruned)))

if (length(taxa_names(func_pruned)) > 15){
  least_abundant_func <- names(sort(taxa_sums(func_pruned),decreasing=T))[16:length(taxa_sums(func_pruned))]
  func_pruned_bar <- merge_taxa(func_pruned, least_abundant_func, 1)
  rownames(otu_table(func_pruned_bar))[rownames(otu_table(func_pruned_bar)) %in% names(sort(taxa_sums(func_pruned),decreasing=T))[16]] <- "Other classified"
} else {
  func_pruned_bar <- func_pruned
}
func_pruned_bar <- transform_sample_counts(func_pruned_bar, function(x) 100 * x/sum(x))
bars_func <- data.frame(t(otu_table(func_pruned_bar)))
colnames(bars_func) <- gsub("\\.", " ", colnames(bars_func))
if (any(colnames(bars_func) %in% "Other classified")) {
  bars_func <- cbind(bars_func[-which(colnames(bars_func) %in% "Other classified")], bars_func[which(colnames(bars_func) %in% "Other classified")])
}

if (!any(colnames(bars_func) %in% "Other classified")) {
  funcolors <- darkcols[2:(length(bars_func)+1)]
} else {
  funcolors <- darkcols[1:(length(bars_func))]
}
rownames(bars_func) <- sample_order
bars_func <- melt(as.matrix(bars_func), varnames=c("Sample", "Function"), value.name="Abundance")
bars_func <- rbind(bars_func[bars_func$Function %in% "Other classified",], bars_func[!bars_func$Function %in% "Other classified",])
bars_func <- bars_func[order(-1:-nrow(bars_func)),]
bars_func$Function <- factor(bars_func$Function, levels = rev(levels(bars_func$Function)))

#bars_func <- rbind(bars_func, data.frame(Sample = empty_levels, Function = levels(bars_func$Function)[1], Abundance = 0))
ggplot(bars_func, aes(x=factor(Sample), y=Abundance, fill=factor(Function))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_bar(stat="identity") + scale_fill_manual(values=rev(darkcols_matam), name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("Function")

#calculate bray-curtis dissimilarities, and present them as NMDS ordinations
metaphlan2_dist <- vegdist(t(metaphlan2_data), distance = "bray")
set.seed(42)
ordu_metaphlan2 <- ordinate(metaphlan2, "NMDS", distance = metaphlan2_dist)
ef_metaphlan2 <- envfit(ordu_metaphlan2,sample_data(metaphlan2),permu=10000)
NMDS_data_metaphlan2 <- data.frame(sample_data(metaphlan2))
NMDS_data_metaphlan2$NMDS1 <- ordu_metaphlan2$points[ ,1]
NMDS_data_metaphlan2$NMDS2 <- ordu_metaphlan2$points[ ,2]

matam_dist <- vegdist(t(matam_data), distance = "bray")
set.seed(42)
ordu_matam <- ordinate(matam, "NMDS", distance = matam_dist)
ef_matam <- envfit(ordu_matam,sample_data(matam),permu=10000)
NMDS_data_matam <- data.frame(sample_data(matam))
NMDS_data_matam$NMDS1 <- ordu_matam$points[ ,1]
NMDS_data_matam$NMDS2 <- ordu_matam$points[ ,2]

matam_func_data <- data.frame(otu_table(func))
matam_func_dist <- vegdist(t(matam_func_data), distance = "bray")
set.seed(42)
ordu_matam_func <- ordinate(func, "NMDS", distance = matam_func_dist)
ef_matam_func <- envfit(ordu_matam_func,sample_data(func),permu=10000)
NMDS_data_matam_func <- data.frame(sample_data(func))
NMDS_data_matam_func$NMDS1 <- ordu_matam_func$points[ ,1]
NMDS_data_matam_func$NMDS2 <- ordu_matam_func$points[ ,2]

humann2_dist <- vegdist(t(humann2_data), distance = "bray")
set.seed(42)
ordu_humann2 <- ordinate(humann2, "NMDS", distance = humann2_dist)
ef_humann2 <- envfit(ordu_humann2,sample_data(humann2),permu=10000)
NMDS_data_humann2 <- data.frame(sample_data(humann2))
NMDS_data_humann2$NMDS1 <- ordu_humann2$points[ ,1]
NMDS_data_humann2$NMDS2 <- ordu_humann2$points[ ,2]

ggplot(data = NMDS_data_metaphlan2, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=Site), size = 5) + labs(color = "Sample")

ggplot(data = NMDS_data_matam, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=factor(Site)), size = 5) + labs(color = "Sample")

ggplot(data = NMDS_data_humann2, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=Site), size = 5) + labs(color = "Sample")

ggplot(data = NMDS_data_matam_func, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=factor(Site)), size = 5) + labs(color = "Sample")

sequencing <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/sample_sequencing_stats.csv",sep="\t")
sequencing <- sequencing[-which(sequencing$sample %in% "total"),]
sequencing <- sequencing[,-c(3,4)]
sequencing$sample <- sample_order
sequencing <- melt(sequencing)

ggplot(sequencing, aes(x=factor(sample), y=value)) + geom_bar(stat="identity", color="black") + scale_x_discrete(limits=(sample_order)) + facet_wrap(~variable, scales="free") + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Sample") + ylab("Abundance (%)") + ggtitle("Sequencing statistics") + scale_y_continuous(labels = scales::comma)

#APD predictions
#get all the sequences
APD_sequences <- readFASTA("http://aps.unmc.edu/AP/APD_AMPs_fasta_DongChuan.fa")
#get all labels in the 5 categories
APD_categories <- c("Antibacterial", "Antiviral", "Antifungal", "Antiparasitic", "Anticancer")
APD_accessions <- list(NULL)
for (i in 1:length(APD_categories)) {
  page <- readLines(paste0("D:/VirtualBox/VirtualBox Share/hazen_metagenome/", APD_categories[i], ".html"))
  page <- gsub(".*> (AP[0-9]{5}) </a>", "\\1", page[grep("AP[0-9]{5}", page)])
  APD_accessions[[i]] <- page
}
#unique_APD <- unique(unlist(APD_accessions))

#limit the data set
APD_data <- APD_data[which(nchar(APD_data) <= 60)]
APD_data <- APD_sequences[which(nchar(APD_sequences) >= 10)]
#APD_data <- APD_data[which(names(APD_data) %in% unique_APD)]
APD_names <- names(APD_data)

#extract features
APD_AAC <- sapply(APD_data, extractAAC)
APD_DC <- sapply(APD_data, extractDC)
APD_data <- data.table(t(rbind(APD_AAC, APD_DC)))
APD_data$accession <- APD_names

#label the sequences
for (i in 1:5) {
  APD_data[, APD_categories[i] := ifelse(APD_data$accession %in% APD_accessions[[i]], 1, 0)]
}

# rtsne_APD <- data.frame(Rtsne(APD_data[, !names(APD_data) %in% APD_categories, with=F], perplexity = 30)$Y)
# rownames(rtsne_APD) <- APD_data$accession
# rtsne_APD <- cbind(rtsne_APD, APD_data[, names(APD_data) %in% APD_categories, with=F])
# ggplot(data = rtsne_APD, aes(x = X1, y = X2)) +  geom_point(aes(color=factor(Anticancer)))


#implement caret

#split data set to training and testing
set.seed(42)
split=0.66
trainIndex <- APD_data[sample(.N, split*nrow(APD_data))]$accession
APD_train <- APD_data[APD_data$accession %in% trainIndex]
APD_train$accession <- list(NULL)
APD_test <- APD_data[!APD_data$accession %in% trainIndex]
APD_test$accession <- list(NULL)

#fit the ECC model (using random forests as the classifier)
fit <- ecc(APD_train[, !APD_categories, with=F], APD_train[, APD_categories, with=F], m = 7, .f = randomForest::randomForest, replace = TRUE, run_parallel = TRUE, prop_subset = 0.8)

#test the unused data
pugs <- predict(fit, APD_test[, !APD_categories, with=F], burn.in = 500, n.iters = 1500, thin = 15, .f = randomForest:::predict.randomForest, type = "prob", run_parallel = TRUE)
test_pred2 <- summary(pugs, type = "prob")

#validate the model
model_stats <- validate_pugs(pugs, APD_test[, APD_categories, with=F])


candidateAP <- readFASTA("D:/VirtualBox/VirtualBox Share/hazen_metagenome/PROKKA_03022018_AP_cand2.faa")
candidateAP <- candidateAP[which(nchar(candidateAP) <= 80)]
candidateAP <- candidateAP[which(nchar(candidateAP) >= 10)]
candidateAP_names <- names(candidateAP)

candidateAP_AAC <- sapply(candidateAP, extractAAC)
candidateAP_DC <- sapply(candidateAP, extractDC)
candidateAP <- data.table(t(rbind(candidateAP_AAC, candidateAP_DC)))
#candidateAP <- cbind(candidateAP, t(data.frame(row.names=APD_categories, rep(0, length(APD_categories)))))
#a <- rbind(candidateAP, APD_test[1:5])


candidateAP_pugs <- predict(fit, candidateAP[, !APD_categories, with=F], burn.in = 500, n.iters = 1500, thin = 15, .f = randomForest:::predict.randomForest, type = "prob", run_parallel = TRUE)
candidateAP_classifications <- data.frame(summary(candidateAP_pugs, type ="prob"))
rownames(candidateAP_classifications) <- candidateAP_names


