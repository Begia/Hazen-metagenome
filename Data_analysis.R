# # installing/loading the package:
# if(!require(installr)) {
#   install.packages("installr"); require(installr)} #load / install+load installr
# 
# # using the package:
# updateR() # this will start the updating process of your R installation.  It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 

packages <- c("ggplot2", "ranger", "svglite", "ggthemes", "phyloseq", "vegan", "Rtsne", "dbscan", "mefa", "caret", "data.table", "plyr", "protr", 
              "data.table", "MLPUGS", "compiler", "DESeq2", "ape", "gggenes", "car", "edarf", "biomformat", "limma")

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

#veganifyOTU function extracted from internal phyloseq function
veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#function for setting levels
f <- function(x) factor(x, levels = unique(x))

#function for decimal points in pd plots
scaleFUN <- function(x) sprintf("%.2f", x)



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
    "Spirochaetes",
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

sample_order <- c("Deep Hole 0.5cm", "Deep Hole 1.5 cm", "Deep Hole 2.5 cm", "Snowgoose Bay 0.5 cm", "Snowgoose Bay 2.0 cm", "Snowgoose Bay 3.5 cm")
# #read in the taxonomic metaphlan2/matam and functional humann2 data
# metaphlan2_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/metaphlan2_merged.txt", sep = "\t")
matam_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/matam_contingency_samples_combined.txt", quote = "", sep = "\t")
# humann2_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/humann2_merged_norm.tsv", sep = "\t")
# 
# #limit metaphlan2 data to species (tips of the tree) or phyla and humann2 data to just the functions (remove tax information)
# metaphlan2_data <- metaphlan2_raw[grep("s__", metaphlan2_raw$ID),]
# metaphlan2_data <- metaphlan2_data[-grep("t__", metaphlan2_data$ID),]
# metaphlan2_data$ID <- as.character(metaphlan2_data$ID)
# metaphlan2_data$ID <- sub("^(.+s__).*?", "", metaphlan2_data$ID)
# metaphlan2_rownames <- metaphlan2_data$ID
# metaphlan2_data[,"ID"] <- list(NULL)
# metaphlan2_data <- data.frame(apply(metaphlan2_data, 2, function(x) as.numeric(as.character(x))))
# #summarize taxa classified at lower levels as a single group
# metaphlan2_data <- rbind(metaphlan2_data, 100-colSums(metaphlan2_data))
# rownames(metaphlan2_data) <- c(metaphlan2_rownames, "other")
# 
# metaphlan2_phyla <- metaphlan2_raw[grep("p__", metaphlan2_raw$ID),]
# metaphlan2_phyla <- metaphlan2_phyla[-grep("o__", metaphlan2_phyla$ID),]
# metaphlan2_prot <- metaphlan2_phyla[grep("proteobacteria", metaphlan2_phyla$ID),]
# metaphlan2_phyla <- metaphlan2_phyla[-grep("c__|Proteobacteria", metaphlan2_phyla$ID),]
# metaphlan2_phyla <- rbind(metaphlan2_phyla, metaphlan2_prot)
# metaphlan2_phyla_rownames <- metaphlan2_phyla$ID
# metaphlan2_phyla[,"ID"] <- list(NULL)
# metaphlan2_phyla <- data.frame(apply(metaphlan2_phyla, 2, function(x) as.numeric(as.character(x))))
# rownames(metaphlan2_phyla) <- metaphlan2_phyla_rownames
# 
# 
# humann2_data <- humann2_raw[-grep("\\|", humann2_raw$X..Pathway),]
# humann2_rownames <- humann2_data$X..Pathway
# humann2_data[,"X..Pathway"] <- list(NULL)
# humann2_data <- data.frame(apply(humann2_data, 2, function(x) as.numeric(as.character(x))))
# colnames(humann2_data) <- sub("_Abundance", "", colnames(humann2_data))
# rownames(humann2_data) <- humann2_rownames

#matam 16S data
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
sample_mapping_data$O2.mgL[which(sample_mapping_data$O2.mgL < 0)] <- 0

#construct phyloseq objects
#metaphlan2 <- phyloseq(otu_table(metaphlan2_data, taxa_are_rows = T), sample_data(sample_mapping_data))
#humann2 <- phyloseq(otu_table(humann2_data, taxa_are_rows = T), sample_data(sample_mapping_data))
matam_mapping_data <- merge_samples(sample_data(sample_mapping_data), group="Sample")
matam_mapping_data$Sample <- sample_order
matam_mapping_data$Site <- c(rep("Deep Hole",3), rep("Snowgoose Bay", 3))
matam <- phyloseq(otu_table(matam_data, taxa_are_rows = T), tax_table(matam_tax), sample_data(matam_mapping_data))
sample_names(matam) <- matam_mapping_data$Sample

#remove chloroplasts and mitochondria, and correct some phyla names
matam <- subset_taxa(matam, !(Domain %in% "Eukaryota"))
matam <- subset_taxa(matam, !(Class %in% "Chloroplast"))
matam <- subset_taxa(matam, !(Family %in% "Mitochondria"))
tax_table(matam)[tax_table(matam)[,2] %in% "Cyanobacteria/Chloroplast",2] <- "Cyanobacteria"
tax_table(matam)[tax_table(matam)[,2] %in% "Woesearchaeota",2] <- "Woesearchaeota (DHVEG-6)"

# #metaphlan2 bars
# bars_metaphlan2 <- data.frame(t(metaphlan2_phyla))
# colnames(bars_metaphlan2) <- sub(".*?_{2,3}(.+)$", "\\1", colnames(bars_metaphlan2))
# bars_metaphlan2$sample <- sub(".*?([0-9]+)$", "\\1", rownames(bars_metaphlan2))
# bars_metaphlan2 <- aggregate(.~sample, data = bars_metaphlan2, mean)
# bars_metaphlan2 <- bars_metaphlan2[order(as.numeric(bars_metaphlan2$sample)),]
# bars_metaphlan2 <- bars_metaphlan2[-1]
# bars_metaphlan2 <- data.frame(apply(bars_metaphlan2, 2, function(x) as.numeric(as.character(x))))
# sample_order <- c("Deep Hole 0.5 cm", "Deep Hole 1.5 cm", "Deep Hole 2.5 cm", "Snowgoose Bay 0.5 cm", "Snowgoose Bay 2.0 cm", "Snowgoose Bay 3.5 cm")
# rownames(bars_metaphlan2) <- sample_order
# bars_metaphlan2 <- melt(as.matrix(bars_metaphlan2), varnames=c("Sample", "Phylum"), value.name="Abundance")
# 
# #plot phyla
# ggplot(bars_metaphlan2, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_metaphlan2$Sample[nchar(as.character(bars_metaphlan2$Sample))!=1]) + 
#   geom_bar(stat="identity") + 
#   scale_fill_manual(values=darkcols, name = "Phylum", guide = guide_legend(reverse = T)) + 
#   theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
#   xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("Taxonomy")

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
bars_matam <- reshape2::melt(as.matrix(bars_matam), varnames=c("Sample", "Phylum"), value.name="Abundance")
bars_matam$Phylum <- factor(bars_matam$Phylum, levels = rev(phyla_matam_levels))

#write out phyla abundances
phyla_percentages <- data.frame(otu_table(phyla_matam))
rownames(phyla_percentages) <- c(tax_table(phyla_matam)[,2])
phyla_percentages$rowmeans <- rowMeans(phyla_percentages)
write.csv(phyla_percentages, "D:/VirtualBox/VirtualBox Share/hazen_metagenome/phyla_abundances.csv", quote=F)

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
colnames(func) <- sample_order

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
bars_func <- reshape2::melt(as.matrix(bars_func), varnames=c("Sample", "Function"), value.name="Abundance")
bars_func <- rbind(bars_func[bars_func$Function %in% "Other classified",], bars_func[!bars_func$Function %in% "Other classified",])
bars_func <- bars_func[order(-1:-nrow(bars_func)),]
bars_func$Function <- factor(bars_func$Function, levels = rev(levels(bars_func$Function)))

# #calculate bray-curtis dissimilarities, and present them as NMDS ordinations
# metaphlan2_dist <- vegdist(t(metaphlan2_data), distance = "bray")
# set.seed(42)
# ordu_metaphlan2 <- ordinate(metaphlan2, "NMDS", distance = metaphlan2_dist)
# ef_metaphlan2 <- envfit(ordu_metaphlan2,sample_data(metaphlan2),permu=10000)
# NMDS_data_metaphlan2 <- data.frame(sample_data(metaphlan2))
# NMDS_data_metaphlan2$NMDS1 <- ordu_metaphlan2$points[ ,1]
# NMDS_data_metaphlan2$NMDS2 <- ordu_metaphlan2$points[ ,2]

matam_distance <- vegdist(t(matam_data), distance = "bray")
set.seed(42)
ordu_matam <- ordinate(matam, "NMDS", distance = matam_distance)
ef_matam <- envfit(ordu_matam,sample_data(matam),permu=10000)
NMDS_data_matam <- data.frame(sample_data(matam))
NMDS_data_matam$NMDS1 <- ordu_matam$points[ ,1]
NMDS_data_matam$NMDS2 <- ordu_matam$points[ ,2]

matam_func_data <- data.frame(otu_table(func))
matam_func_distance <- vegdist(t(matam_func_data), distance = "bray")
set.seed(42)
ordu_matam_func <- ordinate(func, "NMDS", distance = matam_func_distance)
ef_matam_func <- envfit(ordu_matam_func,sample_data(func),permu=10000)
NMDS_data_matam_func <- data.frame(sample_data(func))
NMDS_data_matam_func$NMDS1 <- ordu_matam_func$points[ ,1]
NMDS_data_matam_func$NMDS2 <- ordu_matam_func$points[ ,2]

# humann2_dist <- vegdist(t(humann2_data), distance = "bray")
# set.seed(42)
# ordu_humann2 <- ordinate(humann2, "NMDS", distance = humann2_dist)
# ef_humann2 <- envfit(ordu_humann2,sample_data(humann2),permu=10000)
# NMDS_data_humann2 <- data.frame(sample_data(humann2))
# NMDS_data_humann2$NMDS1 <- ordu_humann2$points[ ,1]
# NMDS_data_humann2$NMDS2 <- ordu_humann2$points[ ,2]

# ggplot(data = NMDS_data_metaphlan2, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=Site), size = 5) + labs(color = "Sample")

ggplot(data = NMDS_data_matam, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=factor(Site)), size = 5) + labs(color = "Sample")

# ggplot(data = NMDS_data_humann2, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=Site), size = 5) + labs(color = "Sample")

ggplot(data = NMDS_data_matam_func, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=factor(Site)), size = 5) + labs(color = "Sample")

sequencing <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/sample_sequencing_stats.csv",sep="\t")
sequencing <- sequencing[-which(sequencing$sample %in% "total"),]
sequencing <- sequencing[,-c(3,4)]
sequencing$sample <- sample_order
sequencing <- reshape2::melt(sequencing)

ggplot(sequencing, aes(x=factor(sample), y=value)) + geom_bar(stat="identity", color="black") + scale_x_discrete(limits=(sample_order)) + facet_wrap(~variable, scales="free") + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Sample") + ylab("Abundance (%)") + ggtitle("Sequencing statistics") + scale_y_continuous(labels = scales::comma)

#genome bin data and functional data analysis
#(sample mean coverage * bin length) / average length of the reads
#calculate relative abundance of bins and reads as % of total reads in a sample
abundances <- list(NULL)
abundances[[1]] <- read.table("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/mean_coverage.txt",sep="\t",row.names=1,header=T)
abundances[[2]] <- read.table("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/all_gene_coverages.txt",sep="\t",row.names=1,header=F)
colnames(abundances[[2]]) <- colnames(abundances[[1]])

genome_gene_lengths <- list(NULL)
genome_gene_lengths[[1]] <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/bins_summary.txt",sep="\t")
genome_gene_lengths[[2]] <- read.table("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/all_gene_calls.txt",sep=" ")
colnames(genome_gene_lengths[[2]]) <- c("annotation", "bin", "gene_call_id", "split_id", "start_bp", "end_bp", "gene_direction")
genome_gene_lengths[[2]] <- genome_gene_lengths[[2]][c("gene_call_id", "annotation", "bin", "split_id", "start_bp", "end_bp", "gene_direction")]
genome_gene_lengths[[2]]$total_length <- genome_gene_lengths[[2]]$end_bp-genome_gene_lengths[[2]]$start_bp

avg_read_lengths <- read.table("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/sample_average_read_lengths.txt",sep="\t")
avg_read_lengths[,2] <- as.numeric(substr(avg_read_lengths[,2], 10, 15))

sample_reads <- data.frame(sample=c(1,3,5,11,14,17), reads=c(132094500, 102955607, 98825210, 129835700, 122863643, 98243245))

rel_abundances <- list(NULL)
for (i in 1:2) {
rel_abundance <- abundances[[i]]
colnames(rel_abundance) <- as.numeric(gsub("\\D", "", colnames(rel_abundance)))
for (j in 1:6) {
  read_hits <- function (a,b) {
    c <- rownames(rel_abundance)[b]
    d <- unique(genome_gene_lengths[[i]]$total_length[which(genome_gene_lengths[[i]][1] == c)])
    a * d
  }
  index <- seq(1:nrow(rel_abundance))
  rel_abundance[,j] <- mapply(read_hits, rel_abundance[,j], index)
  read_length <- avg_read_lengths[which(avg_read_lengths[,1] == colnames(rel_abundance)[j]),2]
  nreads <- sample_reads[which(sample_reads[,1] == colnames(rel_abundance)[j]),2]
  rel_abundance[,j] <- ((rel_abundance[,j] / read_length) / nreads)
  rel_abundance[,j] <- (rel_abundance[,j] / colSums(rel_abundance)[j]) * 100
}
rel_abundance <- rel_abundance[,order(as.numeric(colnames(rel_abundance)))]
colnames(rel_abundance) <- sample_order
rel_abundances[[i]] <- rel_abundance
}


#construct phyloseq objects for (checkm) taxonomy, and metacyc/kegg pathways

phylogeny <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

checkm_taxonomy_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/tree_qa_lineage.out",sep="\t")
checkm_taxonomy <- as.character(checkm_taxonomy_raw$Taxonomy..contained.)
checkm_taxonomy <- gsub("[a-z]__", "", checkm_taxonomy)
checkm_taxonomy <- strsplit(checkm_taxonomy,";")
checkm_taxonomy <- as.matrix(rbind.fill(lapply(checkm_taxonomy,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})))
rownames(checkm_taxonomy) <- gsub("-contigs", "", checkm_taxonomy_raw$Bin.Id)
colnames(checkm_taxonomy) <- phylogeny
checkm_tree <- read.tree("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/tree_qa.tre")
checkm_tree$tip.label <- gsub("-contigs", "", checkm_tree$tip.label)

rel_abundance_pathways <- rel_abundances[[1]]
rel_abundance_pathways$bin <- rownames(rel_abundances[[1]])
#quality_rel_abundance_pathways <- rel_abundance_pathways[which(rel_abundance_pathways$bin %in% quality_bins$bins),]


colnames_pathways <- c("bin", "accession", "pathway")

kegg_pathways <- read.table("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/kegg_pathways.txt",sep="\t")
colnames(kegg_pathways) <- colnames_pathways
# quality_kegg_pathways <- kegg_pathways
kegg_pathways <- merge(kegg_pathways, rel_abundance_pathways, by="bin")
kegg_pathways2 <- aggregate(kegg_pathways[,4:9], by=list(kegg_pathways$pathway), FUN=sum)
rownames(kegg_pathways2) <- kegg_pathways2$Group.1
kegg_pathways2 <- kegg_pathways2[,-1]
# 
# quality_kegg_pathways <- merge(quality_kegg_pathways, quality_rel_abundance_pathways, by="bin")
# quality_kegg_pathways <- aggregate(quality_kegg_pathways[,4:9], by=list(quality_kegg_pathways$pathway), FUN=sum)
# rownames(quality_kegg_pathways) <- quality_kegg_pathways$Group.1
# quality_kegg_pathways <- quality_kegg_pathways[,-1]

metacyc_pathways <- readLines("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/metacyc_pathways.txt")
metacyc_pathways <- strsplit(metacyc_pathways,"\t")
metacyc_pathways <- do.call(rbind.data.frame, metacyc_pathways)[1:3]
colnames(metacyc_pathways) <- colnames_pathways
# quality_metacyc_pathways <- metacyc_pathways
metacyc_pathways <- merge(metacyc_pathways, rel_abundance_pathways, by="bin")
metacyc_pathways2 <- aggregate(metacyc_pathways[,4:9], by=list(metacyc_pathways$pathway), FUN=sum)
rownames(metacyc_pathways2) <- metacyc_pathways2$Group.1
metacyc_pathways2 <- metacyc_pathways2[,-1]

# quality_metacyc_pathways <- merge(quality_metacyc_pathways, quality_rel_abundance_pathways, by="bin")
# quality_metacyc_pathways <- aggregate(quality_metacyc_pathways[,4:9], by=list(quality_metacyc_pathways$pathway), FUN=sum)
# rownames(quality_metacyc_pathways) <- quality_metacyc_pathways$Group.1
# quality_metacyc_pathways <- quality_metacyc_pathways[,-1]

marker_genes <- genome_gene_lengths[[2]]
marker_gene_annotation <- read.table("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/individual_gene_annotations.txt", sep = "\t")
colnames(marker_gene_annotation) <- c("annotation", "gene", "process", "category")
marker_gene_annotation[] <- lapply(marker_gene_annotation, factor)
marker_gene_annotation[] <- lapply(marker_gene_annotation, f)
marker_gene_abundances <- rel_abundances[[2]]
marker_gene_abundances$gene_call_id <- rownames(marker_gene_abundances)
marker_genes2 <- merge(marker_genes, marker_gene_annotation, by="annotation")
marker_genes2 <- merge(marker_genes2, marker_gene_abundances, by="gene_call_id")
annotation_abundances <- aggregate(marker_genes2[,12:17], by=list(marker_genes2$gene), FUN=sum)
rownames(annotation_abundances) <- annotation_abundances$Group.1
annotation_abundances <- annotation_abundances[-1]
annotation_tax_table <- unique(marker_genes2[c(9:11)])
annotation_tax_table <- annotation_tax_table[order(annotation_tax_table$gene),]
rownames(annotation_tax_table) <- annotation_tax_table$gene
annotation_tax_table <- annotation_tax_table[-1]

resfam_genes <- genome_gene_lengths[[2]]
resfam_gene_annotation_raw <- read.csv("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/resfams_annotations.csv", sep = "\t")
resfam_gene_annotation <- resfam_gene_annotation_raw[,c(1:3, 9)]
colnames(resfam_gene_annotation) <- c("annotation", "gene", "description", "category")
resfam_gene_annotation[] <- lapply(resfam_gene_annotation, factor)
resfam_gene_annotation[] <- lapply(resfam_gene_annotation, f)
resfam_gene_abundances <- rel_abundances[[2]]
resfam_gene_abundances$gene_call_id <- rownames(resfam_gene_abundances)
resfam_genes2 <- merge(resfam_genes, resfam_gene_annotation, by="annotation")
resfam_genes2 <- merge(resfam_genes2, resfam_gene_abundances, by="gene_call_id")
resfam_annotation_abundances <- aggregate(resfam_genes2[,12:17], by=list(resfam_genes2$gene), FUN=sum)
rownames(resfam_annotation_abundances) <- resfam_annotation_abundances$Group.1
resfam_annotation_abundances <- resfam_annotation_abundances[-1]
resfam_annotation_tax_table <- unique(resfam_genes2[c(9:11)])
resfam_annotation_tax_table <- resfam_annotation_tax_table[order(resfam_annotation_tax_table$gene),]
rownames(resfam_annotation_tax_table) <- resfam_annotation_tax_table$gene
resfam_annotation_tax_table <- resfam_annotation_tax_table[-1]


#process_abundances <- aggregate(marker_genes2[,11:16], by=list(marker_genes2$process), FUN=sum)


antismash_clusters <- readLines("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/antismash_clusters.out")
antismash_clusters <- data.frame(stringr::str_split_fixed((sub("\".*", "", antismash_clusters)), " ", 2) )
colnames(antismash_clusters) <- c("bin", "cluster")
antismash_clusters <- merge(antismash_clusters, rel_abundance_pathways, by="bin")
antismash_clusters2 <- aggregate(antismash_clusters[,3:8], by=list(antismash_clusters$cluster), FUN=sum)
rownames(antismash_clusters2) <- antismash_clusters2$Group.1
antismash_clusters2 <- antismash_clusters2[,-1]

checkm_phylo <- phyloseq(otu_table(rel_abundances[[1]],taxa_are_rows=T), tax_table(checkm_taxonomy), sample_data(matam), checkm_tree)
checkm_phylo <- transform_sample_counts(checkm_phylo, function(x) 100 * x/sum(x))
kegg_phylo <- phyloseq(otu_table(kegg_pathways2,taxa_are_rows=T), sample_data(matam))
metacyc_phylo <- phyloseq(otu_table(metacyc_pathways2,taxa_are_rows=T), sample_data(matam))
antismash_phylo <- phyloseq(otu_table(antismash_clusters2,taxa_are_rows=T), sample_data(matam))
marker_phylo <- phyloseq(otu_table(annotation_abundances,taxa_are_rows=T), tax_table(as.matrix(rev(annotation_tax_table))), sample_data(matam))
resfam_phylo <- phyloseq(otu_table(resfam_annotation_abundances,taxa_are_rows=T), tax_table(as.matrix(rev(resfam_annotation_tax_table))), sample_data(matam))

#quality_kegg_phylo<- phyloseq(otu_table(quality_kegg_pathways,taxa_are_rows=T), sample_data(matam))
#quality_metacyc_phylo<- phyloseq(otu_table(quality_metacyc_pathways,taxa_are_rows=T), sample_data(matam))

checkm_distance <- DPCoA(checkm_phylo)
ordu_checkm <- ordinate(checkm_phylo, "NMDS", distance = checkm_distance$RaoDis)
continuous_sample_data <- sample_data(checkm_phylo)[,3:10]
checkm_env_distance <- vegdist(continuous_sample_data, method="euclidean")
checkm_mantel <- mantel(checkm_distance$RaoDis, checkm_env_distance, method="pearson", permutations=10000)
ef_checkm <- envfit(ordu_checkm,sample_data(checkm_phylo),permu=10000)

# quality_checkm_distance <- DPCoA(quality_checkm_phylo)
# ordu_quality_checkm <- ordinate(quality_checkm_phylo, "NMDS", distance = quality_checkm_distance$RaoDis)
# quality_checkm_env_distance <- vegdist(continuous_sample_data, method="euclidean")
# quality_checkm_mantel <- mantel(quality_checkm_distance$RaoDis, quality_checkm_env_distance, method="pearson", permutations=10000)
# ef_quality_checkm <- envfit(ordu_quality_checkm,sample_data(quality_checkm_phylo),permu=10000)

kegg_distance <- vegdist(wisconsin(sqrt(veganifyOTU(kegg_phylo))), distance = "bray")
ordu_kegg <- ordinate(kegg_phylo, "NMDS", distance = kegg_distance)
kegg_env_distance <- vegdist(continuous_sample_data, method="euclidean")
kegg_mantel <- mantel(kegg_distance, kegg_env_distance, method="pearson", permutations=10000)
ef_kegg <- envfit(ordu_kegg,sample_data(kegg_phylo),permu=10000)

# quality_kegg_distance <- vegdist(wisconsin(sqrt(veganifyOTU(quality_kegg_phylo))), distance = "bray")
# ordu_quality_kegg <- ordinate(quality_kegg_phylo, "NMDS", distance = quality_kegg_distance)
# quality_kegg_env_distance <- vegdist(continuous_sample_data, method="euclidean")
# quality_kegg_mantel <- mantel(quality_kegg_distance, quality_kegg_env_distance, method="pearson", permutations=10000)
# ef_quality_kegg <- envfit(ordu_quality_kegg,sample_data(quality_kegg_phylo),permu=10000)

metacyc_distance <- vegdist(wisconsin(sqrt(veganifyOTU(metacyc_phylo))), distance = "bray")
ordu_metacyc <- ordinate(metacyc_phylo, "NMDS", distance = metacyc_distance)
metacyc_env_distance <- vegdist(continuous_sample_data, method="euclidean")
metacyc_mantel <- mantel(metacyc_distance, metacyc_env_distance, method="pearson", permutations=10000)
ef_metacyc <- envfit(ordu_metacyc,sample_data(metacyc_phylo),permu=10000)

# quality_metacyc_distance <- vegdist(wisconsin(sqrt(veganifyOTU(quality_metacyc_phylo))), distance = "bray")
# ordu_quality_metacyc <- ordinate(quality_metacyc_phylo, "NMDS", distance = quality_metacyc_distance)
# quality_metacyc_env_distance <- vegdist(continuous_sample_data, method="euclidean")
# quality_metacyc_mantel <- mantel(quality_metacyc_distance, quality_metacyc_env_distance, method="pearson", permutations=10000)
# ef_quality_metacyc <- envfit(ordu_quality_metacyc,sample_data(quality_metacyc_phylo),permu=10000)

cluster_palette <- c("#03ff00",
                     "#d660cd",
                     "#00e1ff",
                     "#ffb100",
                     "#e65141",
                     "#9800ff")




dataset_names <- c("matam", "matam_func", "checkm", "kegg", "metacyc")
rtsne_sample_data <- sample_data(matam)
rtsne_sample_names <- sample_names(matam)
for (i in 1:length(dataset_names)) {
set.seed(42)
if (dataset_names[i] %in% "checkm") {
  rtsne_data <- data.frame(Rtsne(eval(parse(text = paste0(dataset_names[i], "_distance$RaoDis"))), is_distance = T, perplexity = 1.6)$Y)
} else {
  rtsne_data <- data.frame(Rtsne(eval(parse(text = paste0(dataset_names[i], "_distance"))), is_distance = T, perplexity = 1.6)$Y)
}
rownames(rtsne_data) <- rtsne_sample_names
rtsne_data_clusters <- hdbscan(rtsne_data, xdist = vegdist(rtsne_data, "euclidean"), minPts = 2)
rtsne_data <- cbind(rtsne_data, cluster = rtsne_data_clusters$cluster, site = rtsne_sample_data$Site, depth = rtsne_sample_data$Depth.cm)
image <- ggplot(data = rtsne_data, aes(x = X1, y = X2)) + stat_ellipse(level = 0.68, aes(color = factor(cluster)), linetype = 2, size = 1.5) +
  geom_point(aes(color = factor(cluster), fill = site, size = depth), shape = 21, stroke = 2) +
  scale_colour_manual(name = "Cluster", values = cluster_palette) +
  scale_size_continuous(range = c(2,6), name = "Sediment depth (cm)") +  scale_fill_manual(values = c("#910a80", #deep hole
                                                                                                      "#a67e48" #snowgoose bay
  ), name = "Site") +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 1), alpha = guide_legend(order = 2), color = guide_legend(override.aes = list(shape = 1, size = 4, linetype = "blank"))) + theme(axis.title = element_blank())
ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/rtsne_", dataset_names[i], ".svg"), plot=image, units="mm", width=225, height=150)
}


#phyla bars for the checkm data
phyla_checkm <- checkm_phylo
if(any(tax_table(phyla_checkm)[,2] %in% "Proteobacteria")) {
  replacement_proteobacteria <- tax_table(phyla_checkm)[tax_table(phyla_checkm)[,2] %in% "Proteobacteria",3]
  replacement_proteobacteria[replacement_proteobacteria %in% NA] <- "Proteobacteria"
  tax_table(phyla_checkm)[tax_table(phyla_checkm)[,2] %in% "Proteobacteria",2] <- replacement_proteobacteria
}
if(any(is.na(tax_table(matam)[,2]))) {
  tax_table(phyla_checkm)[which(is.na(tax_table(phyla_checkm)[,2])),2] <- tax_table(phyla_checkm)[which(is.na(tax_table(phyla_checkm)[,2])),1]
}
#merge phyla_checkm that have less than 1% abundance (or unknown Phylum) in the data set to "other"
phyla_checkm <- tax_glom(phyla_checkm, "Phylum")
least_abundant_phyla_checkm <- names(taxa_sums(phyla_checkm))[taxa_sums(phyla_checkm)/sum(taxa_sums(phyla_checkm))*100 < 1]
phyla_checkm <- merge_taxa(phyla_checkm, least_abundant_phyla_checkm, 1)
tax_table(phyla_checkm)[tax_table(phyla_checkm)[,2] %in% NA,2] <- "Other"
#phyla_checkm <- merge_samples(phyla_checkm, group="Sample", fun=mean)
phyla_checkm <- transform_sample_counts(phyla_checkm, function(x) 100 * x/sum(x))
bars_checkm <- data.frame(t(otu_table(phyla_checkm)))

#darkcols_checkm <- darkcols[(length(darkcols)+1-ncol(bars_checkm)):length(darkcols)]
colnames(bars_checkm) <- c(tax_table(phyla_checkm)[,2])
phyla_checkm_levels <- phyla_levels_all[phyla_levels_all %in% colnames(bars_checkm)]
bars_checkm <- bars_checkm[,phyla_checkm_levels]
bars_checkm <- reshape2::melt(as.matrix(bars_checkm), varnames=c("Sample", "Phylum"), value.name="Abundance")
bars_checkm$Phylum <- factor(bars_checkm$Phylum, levels = rev(phyla_checkm_levels))
#match the colors of phyla
darkcols_checkm <- data.frame(Phylum=phyla_matam_levels, Color=darkcols_matam)
matchcolors_checkm <- darkcols_checkm[which(darkcols_checkm$Phylum %in% phyla_checkm_levels),]
matchcolors_checkm <- rbind(matchcolors_checkm, data.frame(Phylum=phyla_checkm_levels[which(!phyla_checkm_levels %in% matchcolors_checkm$Phylum)], Color=c("#638ccc","#87a141")))
matchcolors_checkm <- matchcolors_checkm[match(phyla_checkm_levels, matchcolors_checkm$Phylum),]
darkcols_checkm <- as.vector(matchcolors_checkm$Color)

#plot phyla for matam and checkm
image <- ggplot(bars_matam, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_matam$Sample[nchar(as.character(bars_matam$Sample))!=1]) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(darkcols_matam), name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("MATAM taxonomy")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/matam_taxonomy_bars.svg", plot=image, units="mm", width=400, height=200)

image <- ggplot(bars_checkm, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_checkm$Sample[nchar(as.character(bars_checkm$Sample))!=1]) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(darkcols_checkm), name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("CheckM taxonomy")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/checkm_taxonomy_bars.svg", plot=image, units="mm", width=400, height=200)

#plot bars for FAPROTAX and heatmaps for kegg, metacyc, antismash, and marker genes
#bars_func <- rbind(bars_func, data.frame(Sample = empty_levels, Function = levels(bars_func$Function)[1], Abundance = 0))
image <- ggplot(bars_func, aes(x=factor(Sample), y=Abundance, fill=factor(Function))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_bar(stat="identity") + scale_fill_manual(values=rev(darkcols_matam), name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("FAPROTAX predictions")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/FAPROTAX_function_bars.svg", plot=image, units="mm", width=400, height=200)

kegg_heatmap <- data.frame(otu_table(kegg_phylo))
colnames(kegg_heatmap) <- sample_order
kegg_order <- readLines("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/kegg_pathways.order")
kegg_levels <- kegg_order[which(kegg_order %in% rownames(kegg_heatmap))]
kegg_heatmap <- kegg_heatmap[match(kegg_levels, rownames(kegg_heatmap)),]
kegg_heatmap <- reshape2::melt(log10(as.matrix(kegg_heatmap)), varnames=c("Pathway", "Sample"), value.name="Abundance")
kegg_heatmap$Pathway <- factor(kegg_heatmap$Pathway)

metacyc_heatmap <- data.frame(otu_table(metacyc_phylo))
colnames(metacyc_heatmap) <- sample_order
metacyc_order <- readLines("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/metacyc_pathways.order")
metacyc_levels <- metacyc_order[which(metacyc_order %in% rownames(metacyc_heatmap))]
metacyc_heatmap <- metacyc_heatmap[match(metacyc_levels, rownames(metacyc_heatmap)),]
metacyc_heatmap <- reshape2::melt(log10(as.matrix(metacyc_heatmap)), varnames=c("Pathway", "Sample"), value.name="Abundance")
metacyc_heatmap$Pathway <- factor(metacyc_heatmap$Pathway)

antismash_heatmap <- data.frame(otu_table(antismash_phylo))
colnames(antismash_heatmap) <- sample_order
antismash_order <- factor(rownames(antismash_clusters2))
antismash_levels <- antismash_order[which(antismash_order %in% rownames(antismash_heatmap))]
antismash_heatmap <- antismash_heatmap[match(antismash_levels, rownames(antismash_heatmap)),]
antismash_heatmap <- reshape2::melt(log10(as.matrix(antismash_heatmap)), varnames=c("Cluster", "Sample"), value.name="Abundance")
antismash_heatmap$Cluster <- factor(antismash_heatmap$Cluster)

marker_heatmap <- data.frame(otu_table(marker_phylo))
colnames(marker_heatmap) <- sample_order
#marker_order <- factor(rownames(marker_clusters2))
#marker_levels <- marker_order[which(marker_order %in% rownames(marker_heatmap))]
#marker_heatmap <- marker_heatmap[match(marker_levels, rownames(marker_heatmap)),]
marker_heatmap <- reshape2::melt(log10(as.matrix(marker_heatmap)), varnames=c("Gene", "Sample"), value.name="Abundance")
marker_heatmap$Gene <- factor(marker_heatmap$Gene, levels=rev(levels(marker_heatmap$Gene)))

resfam_heatmap <- data.frame(otu_table(resfam_phylo))
colnames(resfam_heatmap) <- sample_order
#resfam_order <- factor(rownames(resfam_clusters2))
#resfam_levels <- resfam_order[which(resfam_order %in% rownames(resfam_heatmap))]
#resfam_heatmap <- resfam_heatmap[match(resfam_levels, rownames(resfam_heatmap)),]
resfam_heatmap <- reshape2::melt(log10(as.matrix(resfam_heatmap)), varnames=c("Resfam", "Sample"), value.name="Abundance")
resfam_heatmap$Resfam <- factor(resfam_heatmap$Resfam, levels=rev(levels(resfam_heatmap$Resfam)))

image <- ggplot(bars_func, aes(x=factor(Sample), y=Abundance, fill=factor(Function))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_bar(stat="identity") + scale_fill_manual(values=rev(darkcols_matam), name = "FAPROTAX group", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("FAPROTAX predictions")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/FAPROTAX_function_bars.svg", plot=image, units="mm", width=400, height=200)

image <- ggplot(kegg_heatmap, aes(x=factor(Sample), y=factor(Pathway))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.border = element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("KEGG pathways") +labs(fill="Log 10 abundance")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/KEGG_pathways_heatmap.svg", plot=image, units="mm", width=400, height=600)

image <- ggplot(metacyc_heatmap, aes(x=factor(Sample), y=factor(Pathway))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 16), panel.border = element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("MetaCyc pathways") +labs(fill="Log 10 abundance")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/metacyc_pathways_heatmap.svg", plot=image, units="mm", width=400, height=1200)

image <- ggplot(antismash_heatmap, aes(x=factor(Sample), y=factor(Cluster))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.border = element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("Antismash clusters") +labs(fill="Log 10 abundance")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/antismash_clusters_heatmap.svg", plot=image, units="mm", width=400, height=600)

image <- ggplot(marker_heatmap, aes(x=factor(Sample), y=factor(Gene))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.border = element_blank(), axis.text.x = element_text(size=22, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("Marker gene abundances") +labs(fill="Log 10 abundance")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/marker_abudance_heatmap.svg", plot=image, units="mm", width=400, height=600)

image <- ggplot(resfam_heatmap, aes(x=factor(Sample), y=factor(Resfam))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.border = element_blank(), axis.text.x = element_text(size=22, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("Resfams abundances") +labs(fill="Log 10 abundance")
ggsave(file = "D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/resfam_abudance_heatmap.svg", plot=image, units="mm", width=400, height=600)


#genome data here

genomes_qa <- readLines("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/qa_trimmed.out")
genomes_qa <- genomes_qa[-c(1:3,(length(genomes_qa)))]
genomes_qa <- data.frame(matrix(unlist(sapply(genomes_qa, function(x) strsplit(x, "\\s+"))), nrow=length(genomes_qa), byrow=T))
genome_header <- c("Bin.Id", "Marker.lineage", "UID", "n.genomes", "n.markers", "n.marker.sets", "0", "1", "2", "3", "4", "5+", "Completeness", "Contamination", "Strain.heterogeneity")
colnames(genomes_qa) <- genome_header
genomes_qa$Bin.Id <- gsub("-contigs", "", genomes_qa$Bin.Id)
genomes_qa <- genomes_qa[order(genomes_qa$Bin.Id),]
genomes <- as.list(as.data.frame(t(genomes_qa[,c(13:15)])))
names(genomes) <- genomes_qa$Bin.Id

genomes2 <- lapply(1:length(genomes),function(x) {
  df1 <- as.data.frame(t(as.numeric(as.character(genomes[[x]]))))
  colnames(df1) <- names(genomes[[x]])
  binname <- names(genomes)[x]
  y <- which(rownames(checkm_taxonomy)==binname)
  df2 <- as.data.frame(kegg_pathways[which(kegg_pathways$bin==binname),3])
  colnames(df2) <- "kegg_pathways"
  df3 <- as.data.frame(metacyc_pathways[which(metacyc_pathways$bin==binname),3])
  colnames(df3) <- "metacyc_pathways"
  df4 <- as.data.frame(antismash_clusters[which(antismash_clusters$bin==binname),2])
  colnames(df4) <- "antismash_clusters"
  df5 <- as.data.frame(resfam_genes2[which(resfam_genes2$bin==binname),c(9,11)])
  colnames(df5) <- c("resfams", "category")
  df6 <- as.data.frame(marker_genes2[which(marker_genes2$bin==binname),9:10])
  colnames(df6) <- c("marker_gene", "process")
  list(stats=df1,taxonomy=checkm_taxonomy[y,],kegg=df2,metacyc=df3,antismash=df4,resfams=df5,marker_genes=df6)
  })
                                                 
names(genomes2) <- genomes_qa$Bin.Id

#export data for DS-FDR
dsfdr_out <- data.frame(otu_table(checkm_phylo))
dsfdr_out <- make_biom(dsfdr_out)
write_biom(dsfdr_out,"D:/VirtualBox/VirtualBox Share/Hazen_metagenome/processed_files/dsfdr_data.biom")


#data sets of only "good quality" bins
quality_bins <- read.table("D:/VirtualBox/VirtualBox Share/hazen_metagenome/processed_files/bins_summary.txt",sep="\t", header=1)
# long_genomes <- quality_bins[which(as.numeric(as.character(quality_bins$total_length)) > 2000000),]
# quality_long_genomes <- genomes_qa$Bin.Id[which(as.numeric(as.character(genomes_qa$Completeness)) > 10)]
# quality_long_genomes <- as.character(long_genomes$bins[which(as.character(long_genomes$bins) %in% quality_long_genomes)])
quality_genomes <- genomes_qa$Bin.Id[which(as.numeric(as.character(genomes_qa$Completeness)) > 50)]
# quality_genomes <- unique(c(quality_genomes,quality_long_genomes))
quality_rel_abundance <- rel_abundances[[1]][which(rownames(rel_abundances[[1]]) %in% quality_genomes),]

quality_checkm_phylo <- prune_taxa((taxa_names(checkm_phylo) %in% rownames(quality_rel_abundance)), checkm_phylo)

#random forest analysis and partial dependence

#make tables of data available for models
# variables <- colnames(sample_data(matam))[-1]
# datasets <- c(rep("quality_checkm_phylo",length(variables)), rep("kegg_phylo",length(variables)), rep("metacyc_phylo",length(variables)), rep("antismash_phylo",length(variables)), rep("marker_phylo",length(variables)))
# modelstorun_frame <- data.frame(set=datasets, model=variables)
# continuous_variables <- variables[-1]
# 
# 
# plot_pd_data <- list(NULL)
# models <- list(NULL)
# all_full_model_errors <- list(NULL)
# all_model_errors_2 <- list(NULL)
# models_2 <- list(NULL)
# pd_data <- list(NULL)
# plot_pd_data <- list(NULL)
# 
# #run random forests for all available data and select the best model with at least 5 predictors
# for (j in 1:nrow(modelstorun_frame)) {
#   #prepare data for input to random forest
#   predictors <-
#     t(otu_table(eval(parse(
#       text = paste0(modelstorun_frame$set[j])
#     ))))
#   if (modelstorun_frame$model[j] %in% continuous_variables) {
#     response <-
#       eval(parse(
#         text = paste0(
#           "sample_data(",
#           modelstorun_frame$set[j],
#           ")$",
#           modelstorun_frame$model[j]
#         )
#       ))
# 
#   } else {
#     response <-
#       as.factor(eval(parse(
#         text = paste0(
#           "sample_data(",
#           modelstorun_frame$set[j],
#           ")$",
#           modelstorun_frame$model[j]
#         )
#       )))
#   }
#   rf_data <- data.frame(response, predictors)
# 
#   #run random forest model
#   set.seed(42)
#   classify <- ranger(response ~ ., data = rf_data, num.trees=5000, importance="impurity")
#   models[[j]] <- classify
# 
#   if (modelstorun_frame$model[j] %in% continuous_variables) {
# 
#     #calculate MSPE estimate and CI for the full model
#     pred1 <- classify$predictions
#     y <-  rf_data$response
#     n <- length(y)
# 
#     # psi is the mean squared prediction error (MSPE) estimate
#     # sigma2 is the estimate of the variance of the MSPE
#     psi1 <- mean((y - pred1)^2)
#     sigma21 <- 1/n * var((y - pred1)^2)
#     # 95% CI:
#     full_MSPE <- c(psi1 - 1.96 * sqrt(sigma21), psi1, psi1 + 1.96 * sqrt(sigma21))
# 
#     #save errors for the full model
#     all_full_model_errors[[j]] <- full_MSPE
# 
#     #sort the data by feature importance
#     importances <- sort(importance(classify), decreasing = T)
# 
#     #add one feature at a time to the model in order of importance and compare MSPEs (consider only features with > mean importance)
#     #when a minimum in MSPE is reached, stop and return the model that was run
# 
#     #reorder the sets
#     rf_data_2 <- rf_data[c("response", names(importances))]
# 
#     comp_classify <- list(NULL)
#     comp_MSPEs <- list(NULL)
#     for (k in 1:length(importances)) {
#       #if all the importances are below the mean (all are 0?) break the loop
#       if (k == 0) {break}
#       new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
#       set.seed(42)
#       comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
#       pred2 <- comp_classify[[k]]$predictions
#       y2 <- new_data$response
#       n2 <- length(y2)
#       psi2 <- mean((y2 - pred2)^2)
#       sigma22 <- 1/n * var((y2 - pred2)^2)
#       # 95% CI:
#       comp_MSPEs[[k]] <- c(psi2 - 1.96 * sqrt(sigma22), psi2, psi2 + 1.96 * sqrt(sigma22))
#     }
#     #find the minimum MSPE and get the feature names
#     comp_MSPEs <- do.call("rbind", comp_MSPEs)
#     nfeatures <- min(which(comp_MSPEs[,2] %in% min(comp_MSPEs[,2])))
#     if (nfeatures < 5) {nfeatures <- 5}
#     feature_names <- names(importances)[1:nfeatures]
#     models_2[[j]] <- comp_classify[[nfeatures]]
#     all_model_errors_2[[j]] <- comp_MSPEs[nfeatures,]
#   } else {
#     #calculate model Kappa for the full model (inherently imbalanced data sets so we are using Cohen's Kappa to compare the models)
#     pred1 <- classify$predictions
#     kappa1 <- postResample(pred1, rf_data$response)[[2]]
#     all_full_model_errors[[j]] <- kappa1
# 
#     #sort the data by feature importance
#     importances <- sort(importance(classify), decreasing = T)
# 
#     #add one feature at a time to the model in order of importance and compare Kappa values
# 
#     #reorder the sets
#     rf_data_2 <- rf_data[c("response", names(importances))]
# 
#     comp_classify <- list(NULL)
#     kappa2 <- NA
#     for (k in 1:length(importances)) {
#       #if all the importances are below the mean (all are 0?) break the loop
#       if (k == 0) {break}
#       new_data <- data.frame(response=rf_data_2$response, rf_data_2[seq(2,k+1)])
#       set.seed(42)
#       comp_classify[[k]] <- ranger(response ~ ., data = new_data, num.trees=5000, importance="impurity")
#       pred2 <- comp_classify[[k]]$predictions
#       kappa2[k] <- postResample(pred2, new_data$response)[[2]]
#       if (kappa2[k] == 1 && k >= 5) {
#         break
#       }
#     }
#     if (length(kappa2) == length(importances)) {
#       k <- min(which(kappa2 %in% max(kappa2)))
#       if (k < 5) {k <- 5}
#     }
#     #this will store either values of the first "1" or the highest Kappa value of all the models run
#     all_model_errors_2[[j]] <- kappa2[k]
#     models_2[[j]] <- comp_classify[[k]]
# 
#   }
# 
#   #make partial dependence plots of the best model using all of the data
#   nfeatures <- models_2[[j]]$num.independent.variables
#   rf_data2 <- rf_data_2[1:(nfeatures+1)]
#   pd <- partial_dependence(models_2[[j]], vars=colnames(rf_data2)[-1], data=rf_data2, n=c(25,nrow(rf_data2)))
#   pd_data[[j]] <- pd
#   plot_pd_data[[j]] <- plot_pd(pd)$data
# 
# }
# 
# #format and save partial dependency ggplots
# pd_ggplots <- list(NULL)
# for (k in 1:nrow(modelstorun_frame)) {
#     plot_pd_data[[k]]$variable <- factor(gsub("\\.", " ", plot_pd_data[[k]]$variable))
#     taxonomy_glom <- eval(parse(text = paste0(modelstorun_frame$set[k])))
#     taxa_names(taxonomy_glom) <- factor(gsub("\\.", " ", taxa_names(taxonomy_glom)))
#     taxonomy_glom <- prune_taxa(levels(plot_pd_data[[k]]$variable), taxonomy_glom)
#     taxonomy_data.frame <- data.frame(variable=rownames(otu_table(taxonomy_glom)))
#     taxonomy_merged <- merge(plot_pd_data[[k]],taxonomy_data.frame,by="variable")
#     if (modelstorun_frame$model[k] %in% continuous_variables) {
#       pd_ggplots[[k]] <- ggplot(data = taxonomy_merged, aes(value, response)) + geom_line(aes(colour=variable), size= 1) +
#         scale_x_continuous() + labs(x="Normalized abundance", y=modelstorun_frame$model[k], colour="Function") +
#         theme(legend.position="none")
#       if (length(unique(taxonomy_merged$variable)) >= 15) {
#         pd_ggplots[[k]] <- pd_ggplots[[k]] + facet_wrap(~variable, scales="free_x", ncol=8)
#       } else if (length(unique(taxonomy_merged$variable)) >= 5) {
#         pd_ggplots[[k]] <- pd_ggplots[[k]] + facet_wrap(~variable, scales="free_x", ncol=4)
#       } else {
#         pd_ggplots[[k]] <- pd_ggplots[[k]] + facet_wrap(~variable, scales="free_x", ncol=2)
#       }
#     } else {
#       pd_ggplots[[k]] <- ggplot(data = taxonomy_merged, aes(value, prediction*100)) + geom_line(aes(colour=variable), size= 1) +
#         scale_x_continuous(labels=scaleFUN) + labs(x="Normalized abundance", y="Prediction (% chance to be classified)", colour="Function") +
#         facet_grid(class~variable, scales="free_x") + theme(legend.position="none")
#     }
#   if ((length(unique(taxonomy_merged$variable)) >= 20 && !(modelstorun_frame$model[k] %in% continuous_variables)) | 
#       (length(unique(taxonomy_merged$variable)) >= 40 && modelstorun_frame$model[k] %in% continuous_variables)) {
#     ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
#            units="mm", width=1200, height=200)
#   } else if (length(unique(taxonomy_merged$variable)) >= 10) {
#     ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
#            units="mm", width=600, height=200)
#   } else {
#     ggsave(file = paste0("D:/VirtualBox/VirtualBox Share/hazen_metagenome/results/", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
#            units="mm", width=300, height=200)
#   }
# }
# 
# models_out <- modelstorun_frame
# models_out$set <- as.vector(models_out$set)
# models_out$n.before <- NA
# models_out$MSPE.before <- NA
# models_out$CI.95.before <- NA
# models_out$r2.before <- NA
# models_out$Kappa.before <- NA
# models_out$OOBPredErr.before <- NA
# models_out$n.after <- NA
# models_out$MSPE.after <- NA
# models_out$CI.95.after <- NA
# models_out$r2.after <- NA
# models_out$Kappa.after <- NA
# models_out$OOBPredErr.after <- NA
# for (n in 1:nrow(models_out)) {
#   models_out$n.before[n] <- models[[n]]$num.independent.variables
#   models_out$n.after[n] <- models_2[[n]]$num.independent.variables
#   if (models_out$model[n] %in% continuous_variables) {
#     models_out$MSPE.before[n] <- round(all_full_model_errors[[n]][2],2)
#     models_out$CI.95.before[n] <- paste0(round(all_full_model_errors[[n]][1],2),"-",round(all_full_model_errors[[n]][3],2))
#     models_out$r2.before[n] <- round(models[[n]]$r.squared,3)
#     models_out$MSPE.after[n] <- round(all_model_errors_2[[n]][2],2)
#     models_out$CI.95.after[n] <- paste0(round(all_model_errors_2[[n]][1],2),"-",round(all_model_errors_2[[n]][3],2))
#     models_out$r2.after[n] <- round(models_2[[n]]$r.squared,3)
#   } else {
#     models_out$Kappa.before[n] <- round(all_full_model_errors[[n]],2)
#     models_out$OOBPredErr.before[n] <- round(models[[n]]$prediction.error * 100,2)
#     models_out$Kappa.after[n] <- round(all_model_errors_2[[n]],2)
#     models_out$OOBPredErr.after[n] <- round(models_2[[n]]$prediction.error * 100,2)
#   }
# }
# 
# 
write.csv(models_out,"D:/VirtualBox/VirtualBox Share/Hazen_metagenome/results/final/results_table.csv",quote=F)


all_levels <- list(kegg_levels, metacyc_levels, antismash_levels, data.frame(resfam_gene_annotation$gene, resfam_gene_annotation$category), data.frame(marker_gene_annotation$gene, marker_gene_annotation$process))
checkm_model_shared <- list(NULL)
checkm_model_present <- list(NULL)
for (i in 1:length(variables)){
  bins <- data.frame(bin=names(importance(models_2[[i]])), value=rep(NA, length(names(importance(models_2[[i]])))))
  if (i == 1) {
    bin_pd_data <- plot_pd_data[[i]][seq(1:(nrow(bins)*25)),]
  } else {
    bin_pd_data <- plot_pd_data[[i]]
      }
  for (j in 1:nrow(bins)){
    if (i == 1) {
      bin_dots_y <- bin_pd_data$prediction[which(as.character(bin_pd_data$variable) == as.character(bins[j,1]))]
    } else {
        bin_dots_y <- bin_pd_data$response[which(as.character(bin_pd_data$variable) == as.character(bins[j,1]))]
        }
    if (bin_dots_y[1] < bin_dots_y[length(bin_dots_y)]) {bins[j,2] <- "positive"} else {bins[j,2] <- "negative"}
  }
  positives <- as.character(bins$bin[which(bins$value %in% "positive")])
  positives <- genomes2[positives]
  positives2 <- list(NULL)
  positives3 <- list(NULL)
  for (k in 1:7) {
    if (k < 3) {
      positives2[[k]] <- do.call(rbind,lapply(positives, '[[', k))
      positives3[[k]] <- do.call(rbind,lapply(positives, '[[', k))
      }
    else if (k < 6) {
      common_elements <- lapply(positives, '[[', k)
      common_elements <- Reduce(intersect, do.call(c,common_elements))
      positives2[[k]] <- common_elements
      all_elements <- data.frame(table(do.call(rbind,lapply(positives, '[[', k))))
      all_elements <- all_elements[match(all_levels[[k-2]], as.character(all_elements$Var1)),]
      colnames(all_elements) <- c("pathway", "count")
      all_elements <- all_elements[which(all_elements$count > 0),]
      if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
      positives3[[k]] <- all_elements
    } else if (k == 6) {
      common_elements <- lapply(positives, '[[', k)
      common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
      common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
      common_elements <- list(common_elements1, common_elements2)
      names(common_elements) <- c("genes", "category")
      positives2[[k]] <- common_elements
      all_elements <- data.frame(do.call(rbind,lapply(positives, '[[', k)))
      all_elements1 <- data.frame(table(all_elements$resfams))
      all_elements1 <- all_elements1[match(as.character(all_levels[[k-2]]$resfam_gene_annotation.gene), as.character(all_elements1$Var1)),]
      colnames(all_elements1) <- c("gene", "count")
      all_elements1 <- all_elements1[which(all_elements1$count > 0),]
      all_elements2 <- unique(data.frame(table(all_elements$category)))
      all_elements2 <- all_elements2[match(unique(as.character(all_levels[[k-2]]$resfam_gene_annotation.category)), as.character(all_elements2$Var1)),]
      colnames(all_elements2) <- c("category", "count")
      all_elements2 <- all_elements2[which(all_elements2$count > 0),]
      all_elements <- list(all_elements1, all_elements2)
      if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
      names(all_elements) <- c("genes", "category")
      positives3[[k]] <- all_elements
    } 
    else {
      common_elements <- lapply(positives, '[[', k)
      common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
      common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
      common_elements <- list(common_elements1, common_elements2)
      names(common_elements) <- c("genes", "processes")
      positives2[[k]] <- common_elements
      all_elements <- data.frame(do.call(rbind,lapply(positives, '[[', k)))
      all_elements1 <- data.frame(table(all_elements$marker_gene))
      all_elements1 <- all_elements1[match(as.character(all_levels[[k-2]]$marker_gene_annotation.gene), as.character(all_elements1$Var1)),]
      colnames(all_elements1) <- c("gene", "count")
      all_elements1 <- all_elements1[which(all_elements1$count > 0),]
      all_elements2 <- unique(data.frame(table(all_elements$process)))
      all_elements2 <- all_elements2[match(unique(as.character(all_levels[[k-2]]$marker_gene_annotation.process)), as.character(all_elements2$Var1)),]
      colnames(all_elements2) <- c("process", "count")
      all_elements2 <- all_elements2[which(all_elements2$count > 0),]
      all_elements <- list(all_elements1, all_elements2)
      if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
      names(all_elements) <- c("genes", "processes")
      positives3[[k]] <- all_elements
    }
  }
  names(positives2) <- names(genomes2[[1]])
  names(positives3) <- names(genomes2[[1]])
  negatives <- as.character(bins$bin[which(bins$value %in% "negative")])
  negatives <- genomes2[negatives]
  negatives2 <- list(NULL)
  negatives3 <- list(NULL)
  for (k in 1:6) {
    if (k < 3) {
      negatives2[[k]] <- do.call(rbind,lapply(negatives, '[[', k))
      negatives3[[k]] <- do.call(rbind,lapply(negatives, '[[', k))
      }
    else if (k < 6) {
      common_elements <- lapply(negatives, '[[', k)
      common_elements <- Reduce(intersect, do.call(c,common_elements))
      negatives2[[k]] <- common_elements
      all_elements <- data.frame(table(do.call(rbind,lapply(negatives, '[[', k))))
      all_elements <- all_elements[match(all_levels[[k-2]], as.character(all_elements$Var1)),]
      colnames(all_elements) <- c("pathway", "count")
      all_elements <- all_elements[which(all_elements$count > 0),]
      if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
      negatives3[[k]] <- all_elements
    } else {
      common_elements <- lapply(negatives, '[[', k)
      common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
      common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
      common_elements <- list(common_elements1, common_elements2)
      names(common_elements) <- c("genes", "processes")
      negatives2[[k]] <- common_elements
      all_elements <- data.frame(do.call(rbind,lapply(negatives, '[[', k)))
      all_elements1 <- data.frame(table(all_elements$marker_gene))
      all_elements1 <- all_elements1[match(as.character(all_levels[[k-2]]$marker_gene_annotation.gene), as.character(all_elements1$Var1)),]
      colnames(all_elements1) <- c("gene", "count")
      all_elements1 <- all_elements1[which(all_elements1$count > 0),]
      all_elements2 <- unique(data.frame(table(all_elements$process)))
      all_elements2 <- all_elements2[match(unique(as.character(all_levels[[k-2]]$marker_gene_annotation.process)), as.character(all_elements2$Var1)),]
      colnames(all_elements2) <- c("process", "count")
      all_elements2 <- all_elements2[which(all_elements2$count > 0),]
      all_elements <- list(all_elements1, all_elements2)
      if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
      names(all_elements) <- c("genes", "processes")
      negatives3[[k]] <- all_elements
    }
  }
  names(negatives2) <- names(genomes2[[1]])
  result_list1 <- list(positives2, negatives2)
  names(result_list1) <- c("positives", "negatives")
  checkm_model_shared[[i]] <- result_list1
  names(checkm_model_shared)[i] <- as.character(modelstorun_frame$model[i])
  names(negatives3) <- names(genomes2[[1]])
  result_list2 <- list(positives3, negatives3)
  names(result_list2) <- c("positives", "negatives")
  checkm_model_present[[i]] <- result_list2
  names(checkm_model_present)[i] <- as.character(modelstorun_frame$model[i])
}

unknown_checkm_phylo <- subset_taxa(quality_checkm_phylo, !(Phylum %in% phyla_levels_all))
plot_tree(quality_checkm_phylo, nodelabf=nodeplotboot(), ladderize="left", color="Phylum", size="Abundance", shape="Site") + coord_polar(theta="y") + 
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) +
  scale_color_brewer(palette="Set3")
radiation_phylo <- subset_taxa(quality_checkm_phylo, is.na(Phylum))

plot_tree(radiation_phylo, nodelabf=nodeplotboot(), ladderize="left", size="Abundance", shape="Site", color="Sample") +  
  guides(shape = guide_legend(override.aes = list(size = 5)), color = guide_legend(override.aes = list(size = 5))) +
  scale_color_brewer(palette="Set3")

top10_names <- list(NULL)
for (i in 1:length(sample_names(quality_checkm_phylo))) {
  top10 <- subset_samples(quality_checkm_phylo, Sample == sample_names(quality_checkm_phylo)[i])
  top10_names[[i]] <- names(sort(taxa_sums(top10), decreasing = TRUE)[1:50])
}
top10_df <- reshape2::melt(top10_names)
top10_df <- as.data.frame.matrix(table(top10_df))
#colnames(top10_df) <- sample_names(quality_checkm_phylo)
sites_venn <- data.frame(Deep.Hole = apply(top10_df[ , c(1:3)] , 1 , paste , collapse = "" ), Snowgoose.Bay = apply(top10_df[ , c(4:6)] , 1 , paste , collapse = "" ))
sites_venn2 <- data.frame(Deep.Hole = ifelse(sites_venn$Deep.Hole %in% "000", 0, 1), Snowgoose.Bay = ifelse(sites_venn$Snowgoose.Bay %in% "000", 0, 1))
top10_venn <- vennCounts(sites_venn2)
vennDiagram(top10_venn)

top10_DH <- subset_samples(quality_checkm_phylo, Site == "Deep Hole")
top10_DH_names <- names(sort(taxa_sums(top10_DH), decreasing = TRUE)[1:10])
top10_DH <- prune_taxa(top10_DH_names, top10_DH)

top10_SB <- subset_samples(quality_checkm_phylo, Site == "Snowgoose Bay")
top10_SB_names <- names(sort(taxa_sums(top10_SB), decreasing = TRUE)[1:10])
top10_SB <- prune_taxa(top10_SB_names, top10_SB)

hgc_splits <- unique(as.character(marker_genes2$split_id[which(marker_genes2$gene %in% c("hgcA","hgcB"))]))
hgc_splits <- marker_genes2[which(marker_genes2$split_id %in% hgc_splits),]
hgc_splits <- hgc_splits[,c("bin", "split_id", "start_bp", "end_bp", "gene_direction", "gene")]
hgc_splits$start_bp <- as.numeric(hgc_splits$start_bp)
hgc_splits$end_bp <- as.numeric(hgc_splits$end_bp)


hgc_splits$direction <- ifelse(hgc_splits$gene_direction == "f", 1, -1)

#hgc_splits <- hgc_splits[which(hgc_splits$gene %in% c("hgcA","hgcB")),]

ggplot(hgc_splits, aes(xmin = start_bp, xmax = end_bp, y = split_id, fill = gene, forward = direction)) + geom_gene_arrow() + facet_wrap(~ split_id, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") + theme_genes() + ggtitle("Spirochaetales bacteria hgcAB genes")

mer_splits <- unique(as.character(marker_genes2$split_id[which(marker_genes2$gene %in% "merA")]))
mer_splits <- marker_genes2[which(marker_genes2$split_id %in% mer_splits),]
mer_splits <- mer_splits[,c("bin", "split_id", "start_bp", "end_bp", "gene_direction", "gene")]

mer_splits$direction <- ifelse(mer_splits$gene_direction == "f", 1, -1)

mer_splits <- mer_splits[which(mer_splits$gene %in% c("merA","merB")),]

ggplot(mer_splits, aes(xmin = start_bp, xmax = end_bp, y = split_id, fill = gene, forward = direction)) + geom_gene_arrow() + facet_wrap(~ split_id, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") + theme_genes()

ggplot2::ggplot(hgc_splits, ggplot2::aes(xmin = start_bp, xmax = end_bp, y =
                                           split_id, fill = gene)) +
  geom_gene_arrow() +
  ggplot2::facet_wrap(~ split_id, scales = "free", ncol = 1) +
  ggplot2::scale_fill_brewer(palette = "Set3")
# #APD predictions
# #get all the sequences
# APD_sequences <- readFASTA("http://aps.unmc.edu/AP/APD_AMPs_fasta_DongChuan.fa")
# #get all labels in the 5 categories
# APD_categories <- c("Antibacterial", "Antiviral", "Antifungal", "Antiparasitic", "Anticancer")
# APD_accessions <- list(NULL)
# for (i in 1:length(APD_categories)) {
#   page <- readLines(paste0("D:/VirtualBox/VirtualBox Share/hazen_metagenome/", APD_categories[i], ".html"))
#   page <- gsub(".*> (AP[0-9]{5}) </a>", "\\1", page[grep("AP[0-9]{5}", page)])
#   APD_accessions[[i]] <- page
# }
# #unique_APD <- unique(unlist(APD_accessions))
# 
# #limit the data set
# APD_data <- APD_data[which(nchar(APD_data) <= 60)]
# APD_data <- APD_sequences[which(nchar(APD_sequences) >= 10)]
# #APD_data <- APD_data[which(names(APD_data) %in% unique_APD)]
# APD_names <- names(APD_data)
# 
# #extract features
# APD_AAC <- sapply(APD_data, extractAAC)
# APD_DC <- sapply(APD_data, extractDC)
# APD_data <- data.table(t(rbind(APD_AAC, APD_DC)))
# APD_data$accession <- APD_names
# 
# #label the sequences
# for (i in 1:5) {
#   APD_data[, APD_categories[i] := ifelse(APD_data$accession %in% APD_accessions[[i]], 1, 0)]
# }
# 
# # rtsne_APD <- data.frame(Rtsne(APD_data[, !names(APD_data) %in% APD_categories, with=F], perplexity = 30)$Y)
# # rownames(rtsne_APD) <- APD_data$accession
# # rtsne_APD <- cbind(rtsne_APD, APD_data[, names(APD_data) %in% APD_categories, with=F])
# # ggplot(data = rtsne_APD, aes(x = X1, y = X2)) +  geom_point(aes(color=factor(Anticancer)))
# 
# 
# #implement caret
# 
# #split data set to training and testing
# set.seed(42)
# split=0.66
# trainIndex <- APD_data[sample(.N, split*nrow(APD_data))]$accession
# APD_train <- APD_data[APD_data$accession %in% trainIndex]
# APD_train$accession <- list(NULL)
# APD_test <- APD_data[!APD_data$accession %in% trainIndex]
# APD_test$accession <- list(NULL)
# 
# #fit the ECC model (using random forests as the classifier)
# fit <- ecc(APD_train[, !APD_categories, with=F], APD_train[, APD_categories, with=F], m = 7, .f = randomForest::randomForest, replace = TRUE, run_parallel = TRUE, prop_subset = 0.8)
# 
# #test the unused data
# pugs <- predict(fit, APD_test[, !APD_categories, with=F], burn.in = 500, n.iters = 1500, thin = 15, .f = randomForest:::predict.randomForest, type = "prob", run_parallel = TRUE)
# test_pred2 <- summary(pugs, type = "prob")
# 
# #validate the model
# model_stats <- validate_pugs(pugs, APD_test[, APD_categories, with=F])
# 
# 
# candidateAP <- readFASTA("D:/VirtualBox/VirtualBox Share/hazen_metagenome/PROKKA_03022018_AP_cand2.faa")
# candidateAP <- candidateAP[which(nchar(candidateAP) <= 80)]
# candidateAP <- candidateAP[which(nchar(candidateAP) >= 10)]
# candidateAP_names <- names(candidateAP)
# 
# candidateAP_AAC <- sapply(candidateAP, extractAAC)
# candidateAP_DC <- sapply(candidateAP, extractDC)
# candidateAP <- data.table(t(rbind(candidateAP_AAC, candidateAP_DC)))
# #candidateAP <- cbind(candidateAP, t(data.frame(row.names=APD_categories, rep(0, length(APD_categories)))))
# #a <- rbind(candidateAP, APD_test[1:5])
# 
# 
# candidateAP_pugs <- predict(fit, candidateAP[, !APD_categories, with=F], burn.in = 500, n.iters = 1500, thin = 15, .f = randomForest:::predict.randomForest, type = "prob", run_parallel = TRUE)
# candidateAP_classifications <- data.frame(summary(candidateAP_pugs, type ="prob"))
# rownames(candidateAP_classifications) <- candidateAP_names
# 
