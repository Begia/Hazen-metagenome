# # installing/loading the package:
# if(!require(installr)) {
#   install.packages("installr"); require(installr)} #load / install+load installr
# 
# # using the package:
# updateR() # this will start the updating process of your R installation.  It will check for newer versions, and if one is available, will guide you through the decisions you'd need to make.

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 

packages <- c("ggplot2", "ranger", "svglite", "ggthemes", "phyloseq", "vegan", "Rtsne", "dbscan", "mefa", "caret", "data.table", "plyr", "protr", 
              "data.table", "MLPUGS", "compiler", "DESeq2", "ape", "gggenes", "car", "edarf", "biomformat", "limma", "seqinr", "taxize", "ggtree", 
              "tidytree", "dplyr")

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
    "Rhodothermaeota",
    "RBG-1 (Zixibacteria)",
    "Candidatus Berkelbacteria",
    "Candidatus Falkowbacteria",
    "Candidatus Kerfeldbacteria",
    "Candidatus Pacebacteria",
    "Candidatus Saccharibacteria",
    "Candidatus Shapirobacteria",
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

treepalette <- c("#368f23",
                 "#9151cf",
                 "#75c134",
                 "#9e38b2",
                 "#48ce60",
                 "#d16ce7",
                 "#3eab2f",
                 "#4667e6",
                 "#b2c634",
                 "#5452c5",
                 "#7dcc57",
                 "#bd31a2",
                 "#6aa729",
                 "#8c75f1",
                 "#91a824",
                 "#4483f2",
                 "#d3ba37",
                 "#6347a6",
                 "#54ae50",
                 "#d33393",
                 "#43d087",
                 "#e83381",
                 "#3fae6e",
                 "#e65ec0",
                 "#3f862d",
                 "#b25cc1",
                 "#8ac66a",
                 "#8c6ad4",
                 "#eaa731",
                 "#3a63bf",
                 "#ee8d27",
                 "#3295e9",
                 "#e75f29",
                 "#34b9e1",
                 "#e0422f",
                 "#45ccac",
                 "#de255d",
                 "#3c954d",
                 "#e482e1",
                 "#217127",
                 "#ef5fa3",
                 "#72c888",
                 "#a82975",
                 "#96b248",
                 "#9159b5",
                 "#a9a134",
                 "#a983e3",
                 "#ab8e18",
                 "#738ce6",
                 "#e2712e",
                 "#529cde",
                 "#c0461f",
                 "#41c4cc",
                 "#d1293b",
                 "#2b9b77",
                 "#ea4b78",
                 "#279e94",
                 "#eb5d68",
                 "#2b6e3f",
                 "#9e3f91",
                 "#719d47",
                 "#bd77c6",
                 "#728828",
                 "#c799e9",
                 "#48762c",
                 "#e679ba",
                 "#516615",
                 "#df97d9",
                 "#8b7b1c",
                 "#5375bf",
                 "#cd7325",
                 "#67b3e4",
                 "#a43819",
                 "#377dad",
                 "#e29944",
                 "#695ba6",
                 "#cebe5e",
                 "#7a488d",
                 "#d3a648",
                 "#3e5e99",
                 "#eb8d4a",
                 "#7b8fca",
                 "#b4791f",
                 "#ada7eb",
                 "#a75410",
                 "#8e74b6",
                 "#8dc485",
                 "#c05191",
                 "#499065",
                 "#d65283",
                 "#7ac5a6",
                 "#c0384c",
                 "#206e54",
                 "#b1325c",
                 "#759e5f",
                 "#a35a9e",
                 "#b2c077",
                 "#924274",
                 "#9a9c57",
                 "#834c77",
                 "#c4ab67",
                 "#b377aa",
                 "#686219",
                 "#e97090",
                 "#5c6a31",
                 "#e694b8",
                 "#895918",
                 "#c0678d",
                 "#c2934b",
                 "#9a5b7b",
                 "#e2a979",
                 "#90425d",
                 "#ef956d",
                 "#97444d",
                 "#896e38",
                 "#a22e2f",
                 "#b77d4e",
                 "#c56267",
                 "#97481c",
                 "#ea908f",
                 "#914d34",
                 "#ee755c",
                 "#bb6f6b",
                 "#c5673d",
                 "#c55247",
                 "#c06e51")
#treepalette <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

sample_order <- c("Deep Hole 0.5 cm", "Deep Hole 1.5 cm", "Deep Hole 2.5 cm", "Snowgoose Bay 0.5 cm", "Snowgoose Bay 2.0 cm", "Snowgoose Bay 3.5 cm")



#Ribosomal Protein alignments
RP_annotations <- read.table("E:/hazen_metagenome/processed_files/ribosomal_proteins/RP_annotations.txt")
colnames(RP_annotations) <- c("Pfam", "RP")

for (RP in RP_annotations$RP) {
alignment <- read.alignment(paste0("E:/hazen_metagenome/processed_files/ribosomal_proteins/", RP, "_mafft_translatorx.nt_ali.fasta"), format="fasta")

#calculate the codons at each position in frame +1
alignment_codons <- list(NULL)
for (i in 1:alignment$nb) {
  alignment_codons[[i]] <- splitseq(s2c(alignment$seq[[i]]), frame=0, word=3)
}
alignment_codons <- do.call(rbind, alignment_codons)

#calculate consensus codons (most common) at each position
consensus_codons <- apply(alignment_codons, 2, function(x) names(which.max(table(x))))

#check which positions have a gap as the most common codon, but limit the tails of sequences from the examined region
gap_columns <- which(grepl("---", consensus_codons[min(which(!grepl("---", consensus_codons))):max(which(!grepl("---", consensus_codons)))]))
#fix the numbering of the gap columns
gap_columns <- gap_columns + min(which(!grepl("---", consensus_codons)))

#find sequences that are bridging these "gap positions"
bad_sequences <- list(NULL)
i <- 1
for (co in gap_columns) {
  bad_sequences[[i]] <- alignment$nam[which(!grepl("---", alignment_codons[,co]))]
  i <- i+1
}

#write out the bridging sequences for gaps which are bridged by < 10% of all sequences 
bridging_sequences <- unique(unlist(bad_sequences[which(lengths(bad_sequences) <= 0.001*alignment$nb)]))
write.table(bridging_sequences, paste0("E:/hazen_metagenome/processed_files/ribosomal_proteins/", RP, "_bridging_sequences.list"), row.names=F, col.names=F, quote=F)
}

# following code run as arrays at CSC cluster
# remove the outliers, realign and trim with trimal
# cd $WRKDIR/ncbi/alignment_files
# module load biokit
# 
# RP=$(awk '{print $2}' /wrk/ruuskan1/DONOTREMOVE/RP_annotations.txt | sed -n "${SLURM_ARRAY_TASK_ID}p")
# 
# faSomeRecords ${RP}_mafft_translatorx.nt_ali.fasta ${RP}_bridging_sequences.list ${RP}_bridging_sequences.fasta
# faSomeRecords -exclude ${RP}_mafft_translatorx.nt_ali.fasta ${RP}_bridging_sequences.list ${RP}_refined_mafft_translatorx.nt_ali.fasta
# sed -i 's/-//g' ${RP}_refined_mafft_translatorx.nt_ali.fasta
# perl /homeappl/home/ruuskan1/translatorx_vLocal.pl -i ${RP}_refined_mafft_translatorx.nt_ali.fasta -o ${RP}_refined_v2_mafft_translatorx -c 11 -p F
# trimal -in ${RP}_refined_v2_mafft_translatorx.nt_ali.fasta -out ${RP}_trimal.fasta -gappyout

# remove duplicates, correct headers and construct approximate ML trees with FastTree
# cd $WRKDIR/ncbi/alignment_files
# module load biokit
#
# for RP in $(awk '{print $2}' /wrk/ruuskan1/DONOTREMOVE/RP_annotations.txt); do dedupe.sh in=${RP}_trimal.fasta out=${RP}_trimal_deduped.fasta ac=f requirematchingnames ignorejunk; done
# ### make sure removing exact duplicates from MAGs solved the issue ###
# grep ">" *_trimal.fasta | sort | uniq -c | sort -nr | awk '{if ($1 > 1) print $0}' 
# grep ">" *_deduped.fasta | sort | uniq -c | sort -nr | awk '{if ($1 > 1) print $0}'
#
# RP=$(awk '{print $2}' /wrk/ruuskan1/DONOTREMOVE/RP_annotations.txt | sed -n "${SLURM_ARRAY_TASK_ID}p")
# FastTreeMP -gtr -gamma -nt ${RP}_trimal_deduped.fasta > ${RP}_trimal.tre
# treeshrink.py -q 0.01 -i ${RP}_trimal.tre
# tr -s "\t" "\n" < ${RP}_trimal_treeshrink/${RP}_trimal_shrunk_RS_0.01.txt > ${RP}_trimal_treeshrink/${RP}_outliers.list
# faSomeRecords -exclude ${RP}_trimal_unique.fasta ${RP}_trimal_treeshrink/${RP}_outliers.list ${RP}_trimal_v2.fasta
# FastTreeMP -gtr -gamma -nt ${RP}_trimal_v2.fasta > ${RP}_trimal_v2.tre

#plot old and new trees together in one figure
#TO DO!!!!!


#concatenate alignments
raw_files <- list(NULL)
all_samples <- list(NULL)
index <- 1
for (RP in RP_annotations$RP) {
  samples <- readLines(paste0("E:/hazen_metagenome/processed_files/ribosomal_proteins/", RP, "_trimal_v2.fasta"))
  raw_files[[index]] <- samples
  samples <- samples[which(grepl(">", samples))]
  samples <- gsub("(.+?)\\_(.*)", "\\2", samples)
  all_samples[[index]] <- samples
  index <- index + 1
}
all_samples <- unique(unlist(all_samples))
all_samples <- all_samples[order(all_samples)]

files <- list(NULL)
index <- 1
for (RP in RP_annotations$RP) {
  alignment <- raw_files[[index]]
  cat_sequences <- list(NULL)
  seq_breaks <- which(grepl(">", alignment))
  for (i in 1:length(seq_breaks)) {
    start_index <- seq_breaks[i]+1
    if (i == length(seq_breaks)) {end_index <- length(alignment)} else {end_index <- seq_breaks[i+1]-1}
    cat_sequences[[i]] <- paste0(alignment[start_index:end_index], collapse="")
    names(cat_sequences)[i] <- gsub("(.+?)\\_(.*)", "\\2",alignment[seq_breaks[i]])
  }
  alignment <- data.frame(sample=names(cat_sequences), sequence=do.call(rbind, cat_sequences))
  missing_samples <- all_samples[-which(all_samples %in% alignment$sample)]
  missing_sequences <- data.frame(sample=missing_samples, sequence=rep(paste(rep("-",nchar(as.character(alignment$sequence[1]))), collapse=''), length(missing_samples)))
  alignment <- rbind(alignment,missing_sequences)
  alignment <- alignment[order(match(alignment$sample, all_samples)),]
  files[[index]] <- alignment[2]
  index <- index +1
}

all_sequences <- do.call(cbind, files)
all_sequences <- apply(all_sequences, 1, function(row) paste(row, collapse=""))
outfile <- c(matrix(c(paste(">",all_samples,sep=""), all_sequences), 2, byrow = T))
#remove NCBI genomes with more than 25% gaps over the whole alignment, and metagenomic bins with more than 50% gaps
alignment_length <- nchar(outfile[2])
too_short_seq_ncbi <- which(nchar(as.character(outfile)) -nchar( gsub("-", "", outfile, fixed=T)) > alignment_length*0.25)
too_short_seq_ncbi <- too_short_seq_ncbi[too_short_seq_ncbi > max(grep("Bin", outfile)) + 1]
too_short_seq_bins <- which(nchar(as.character(outfile)) -nchar( gsub("-", "", outfile, fixed=T)) > alignment_length*0.5)
too_short_seq_bins <- too_short_seq_bins[too_short_seq_bins <= max(grep("Bin", outfile)) + 1]
too_short_seq <- c(too_short_seq_bins, too_short_seq_ncbi)
outfile <- outfile[-sort(c(too_short_seq - 1, too_short_seq))]
#drop the CPR genome GCA_000989185 because a complete genome exists: GCA_001029715
GCA_000989185_id <- which(grepl("GCA_000989185", outfile))
outfile <- outfile[-c(GCA_000989185_id, GCA_000989185_id+1)]
write.table(outfile, file="E:/hazen_metagenome/processed_files/ribosomal_proteins/concatenated_RP.fasta", quote=F, row.names=F, col.names=F)

#make trees with FastTree: 
#FastTreeMP -gtr -gamma -nt concatenated_RP.fasta > concatenated_RP.tre

#read in the tree and genomes files
concatenated_tree <- read.tree(file="E:/hazen_metagenome/processed_files/ribosomal_proteins/concatenated_RP.tre")
#root on Archaea
concatenated_tree <- root(concatenated_tree, node=7808, resolve.root=T)

ncbi_taxonomy <- read.table("E:/hazen_metagenome/processed_files/ribosomal_proteins/all_genomes_taxonomy.txt")
colnames(ncbi_taxonomy) <- c("label", "ncbi_id")
ncbi_taxonomy <- ncbi_taxonomy[which(ncbi_taxonomy$label %in% concatenated_tree$tip.label),]

#extract taxonomy for the NCBI genomes and save the file (so it can later be read in)
# taxize_classification <- classification(ncbi_taxonomy$ncbi_id, db="ncbi")
# taxize_classification <- do.call(rbind, taxize_classification)
# taxize_classification$label <- gsub(".[0-9]+$", "", rownames(taxize_classification))
# taxize_classification <- taxize_classification[-3]
# taxize_classification <- taxize_classification[-which(taxize_classification$name %in% "cellular organisms"),]
# taxize_classification2 <- reshape(taxize_classification, timevar = "rank", idvar = "label", direction = "wide")[1:9]
# rownames(taxize_classification2) <- taxize_classification2$label
# taxize_classification2 <- taxize_classification2[,-which(names(taxize_classification2) %in% c("name.no rank", "label"))]
# colnames(taxize_classification2) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# 
# write.table(taxize_classification2, file="E:/hazen_metagenome/processed_files/ribosomal_proteins/NCBI_genomes_tax_table.txt", quote=F, sep="\t")

#construct a tax table and otu table for the tree of both ncbi and metagenome genomes
ncbi_genomes_taxonomy <- read.csv("E:/hazen_metagenome/processed_files/ribosomal_proteins/NCBI_genomes_tax_table.txt", quote="", sep="\t")
ncbi_genomes_taxonomy$ncbi_id <- rownames(ncbi_genomes_taxonomy)
ncbi_genomes_taxonomy <- merge(ncbi_genomes_taxonomy, ncbi_taxonomy, by="ncbi_id")
bin_tip_names <- concatenated_tree$tip.label[which(!grepl("GCA_", concatenated_tree$tip.label))]
ncbi_genomes_taxonomy <- rbind.fill(ncbi_genomes_taxonomy, data.frame(label=bin_tip_names, Species=bin_tip_names))
ncbi_genomes_taxonomy$source <- ifelse(grepl("Bin_", ncbi_genomes_taxonomy$label),"metagenome","repository") 

#construct an OTU table
# ncbi_genomes_otu_table <- data.frame(row.names=rownames(ncbi_genomes_taxonomy), NCBI=rep(0, nrow(ncbi_genomes_taxonomy)), metagenome=rep(0, nrow(ncbi_genomes_taxonomy)))
# ncbi_genomes_otu_table$NCBI[which(grepl("GCA", rownames(ncbi_genomes_otu_table)))] <- 1
# ncbi_genomes_otu_table$metagenome[which(grepl("Bin", rownames(ncbi_genomes_otu_table)))] <- 1
# 
# #construct phyloseq object
# ncbi_genomes_taxonomy <- ncbi_genomes_taxonomy[-8]
# ncbi_genomes <- phyloseq(otu_table(ncbi_genomes_otu_table, taxa_are_rows = T), tax_table(as.matrix(ncbi_genomes_taxonomy)), phy_tree(concatenated_tree))

#plot a tree with ggtree
new_ggtree <- as_data_frame(phy_tree(concatenated_tree))
all_taxonomy <- ncbi_genomes_taxonomy
all_taxonomy$fulltax <- apply(all_taxonomy[,3:8], 1, paste, collapse = "_")
all_taxonomy$fulltax <- gsub("NA_", "", all_taxonomy$fulltax)
new_ggtree <- full_join(new_ggtree, all_taxonomy, by="label")

#change the phylum grouping to metagenome for MAGS and to class level for Proteobacteria
new_ggtree$Phylum_grp <- as.character(new_ggtree$Phylum)
new_ggtree$Phylum_grp[grepl("metagenome", new_ggtree$source)] <- "metagenome"
new_ggtree$Phylum_grp[new_ggtree$Phylum_grp %in% "Proteobacteria"] <- as.character(new_ggtree$Class[new_ggtree$Phylum_grp %in% "Proteobacteria"])

new_ggtree_groupinfo <- split(new_ggtree$label, new_ggtree$Phylum_grp)
new_ggtree <- ggtree::groupOTU(new_ggtree, new_ggtree_groupinfo)

new_ggtree2 <- as.treedata(new_ggtree)

treepalette2 <- treepalette
treepalette2[which(levels(new_ggtree$group) == "metagenome")] <- "black"
treepalette2[which(levels(new_ggtree$group) == "0")] <- "darkgrey"

phylum_to_color_map <- data.frame(Phylum=levels(new_ggtree$group), color=treepalette2)

#node_annotation <- read.csv("E:/hazen_metagenome/processed_files/ribosomal_proteins/RP_tree_node_annotation.csv", quote="", sep="\t")

image <- ggtree(new_ggtree2, aes(color=group), ladderize = T, size=0.1) + geom_text2(aes(subset=!isTip, label=paste0(label, "_", node)), hjust=-.3, size=0.2) + 
  geom_tiplab(aes(label=new_ggtree$fulltax), size=0.2) + 
  geom_tippoint(aes(alpha=new_ggtree$source), color="red", shape=20, size=0.2) + 
  scale_alpha_manual(values=c("1", "0")) +
  geom_treescale(linesize=0.2, fontsize=1, offset=1) +
  scale_color_manual(values=c(treepalette2))

ggsave(file = "E:/hazen_metagenome/results/CPR_phyla_nodes_supports.pdf", plot=image, units="cm", width=15, height=70, limitsize = F, scale=2.5)

# nodestoplot <- new_ggtree[c("node", "label")]
# nodestoplot <- nodestoplot[which(is.na(new_ggtree$source)),]
# nodestoplot <- nodestoplot[-which(nodestoplot$node %in% node_annotation$node),]


# image <- ggtree(new_ggtree2, aes(color=group), ladderize = T, size=0.1) + geom_text2(aes(subset=!isTip & !(node %in% node_annotation$node), label=label), hjust=-.3, size=0.2) + 
#   geom_tiplab(aes(label=new_ggtree$fulltax), size=0.2) + 
#   geom_tippoint(aes(alpha=new_ggtree$source), color="red", shape=20, size=0.2) + 
#   scale_alpha_manual(values=c("1", "0")) +
#   geom_treescale(linesize=0.2, fontsize=1, offset=1) +
#   scale_color_manual(values=c(treepalette2))
# 
# for (i in 1:length(node_annotation$node)) {
#   image <- collapse(image, node=node_annotation$node[i])
#   image <- image + geom_cladelabel(node=node_annotation$node[i], node_annotation$clade[i], fontsize=0.3, color="black")
# }
# 
# image <- image + geom_point2(aes(subset=(node %in% node_annotation$node)), size=0.2, color="darkgrey")
image <- ggtree(new_ggtree2, aes(color=group), ladderize = T, size=0.1) + #geom_text2(aes(subset=!isTip, label=label), hjust=-.3, size=0.2) + 
  geom_nodelab(aes(subset=(!isTip & as.numeric(label) > 0.75)), label="\u25CF", size=2, color="black", hjust=6) + 
  geom_tiplab(aes(label=new_ggtree$fulltax), size=0.2) + 
  geom_tippoint(aes(alpha=new_ggtree$source), color="red", shape=20, size=0.2) + 
  scale_alpha_manual(values=c("1", "0")) +
  geom_treescale(linesize=0.2, fontsize=1, offset=1) +
  scale_color_manual(values=c(treepalette2))

ggsave(file = "E:/hazen_metagenome/results/CPR_phyla_supports.pdf", plot=image, units="cm", width=15, height=80, limitsize = F, scale=2.5, device=cairo_pdf)

image <- ggtree(new_ggtree2, layout="circular", aes(color=group), ladderize = T, size=0.1) + 
  #geom_nodelab(aes(subset=(!isTip & as.numeric(label) > 0.9)), label="\u25CF", size=0.2, hjust=0, color="black") + 
  geom_tippoint(aes(alpha=new_ggtree$source), color="red", shape=18, size=0.2) + 
  scale_alpha_manual(values=c("1", "0")) +
  geom_treescale(linesize=0.2, fontsize=1, offset=1, x=-2) +
  scale_color_manual(values=c(treepalette2))

# for (i in 1:length(node_annotation$node)) {
#   image <- collapse(image, node=node_annotation$node[i])
#   image <- image + geom_cladelabel(node=node_annotation$node[i], node_annotation$clade[i], fontsize=0.3, color="black")
# }
# 
# image <- image + geom_point2(aes(subset=(node %in% node_annotation$node)), size=0.2, color="darkgrey")
image <- image + xlim(-0.3, NA)
image <- image %>% scaleClade(11641, 0.25)
image <- image %>% scaleClade(8094, 0.25)
image <- image %>% scaleClade(9248, 0.25)
image <- image %>% scaleClade(9780, 0.25)
image <- image %>% scaleClade(9504, 0.25)
image <- image %>% scaleClade(10449, 0.25)
image <- image %>% scaleClade(10212, 0.25)
image <- image %>% scaleClade(7249, 0.25)
image <- image %>% scaleClade(6841, 0.25)
image <- image %>% scaleClade(6076, 0.25)
image <- image %>% scaleClade(6552, 0.25)
image <- image %>% scaleClade(6213, 0.25)
image <- image %>% scaleClade(6032, 4)
image <- image %>% scaleClade(6705, 4)
image <- image %>% scaleClade(6163, 4)
image <- image %>% scaleClade(6187, 4)
image <- image %>% scaleClade(10430, 4)

ggsave(file = "E:/hazen_metagenome/results/CPR_phyla_circular.pdf", plot=image, units="cm", width=10, height=10, limitsize = F, scale=0.5, device=cairo_pdf)

branchlength_comparison <- new_ggtree[-which(is.na(new_ggtree$source)),]
repository_non_0 <- branchlength_comparison$branch.length[which(branchlength_comparison$source %in% "repository")]
repository_non_0 <- repository_non_0[which(repository_non_0 > 0)]
repository_quantile <- quantile(repository_non_0, probs=0.95)

image <- ggplot(branchlength_comparison, aes(x = branch.length, fill = source)) + geom_density(adjust=1/5,alpha = 0.5) + 
  geom_vline(aes(xintercept=repository_quantile), linetype="dashed", color="red", size=1) + ggtitle("Branch length distribution in genomes") +
  xlab("Branch length") + ylab("Density") + scale_fill_discrete(name="Source", labels=c("Metagenome", "Repository"))

ggsave(file = "E:/hazen_metagenome/results/branchlength_comparison.svg", plot=image, units="cm", width=10, height=7.5, limitsize = F, scale=2)

over_the_line <- branchlength_comparison[which(branchlength_comparison$branch.length > repository_quantile),]
over_the_line_MAG <- over_the_line[which(over_the_line$source %in% "metagenome"),]
over_the_line_ncbi <- over_the_line[which(over_the_line$source %in% "repository"),]

# #read in the taxonomic metaphlan2/matam and functional humann2 data
matam_raw <- read.csv("E:/hazen_metagenome/matam_contingency_samples_combined.txt", quote = "", sep = "\t")

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
sample_mapping_data <- import_qiime_sample_data("E:/hazen_metagenome/processed_files/sample_data.csv")
microprobe_data <- read.csv("E:/hazen_metagenome/processed_files/microprobes.csv",sep="\t")
microprobe_data <- transform(microprobe_data, bin = cut(Depth.mm, breaks = seq(0,35,by=5), include.lowest = T))
microprobe_data_means <- ddply(microprobe_data, .(Core,bin), numcolwise(median))[-c(1:3)]
microprobe_data_means <- cbind(microprobe_data_means, bin=rep(unique(microprobe_data$bin),2))
microprobe_data_sds <- ddply(microprobe_data, .(Core,bin), numcolwise(sd))[,-c(1:3)]
orig_microprobe_data <- microprobe_data_means
orig_microprobe_sds <-  microprobe_data_sds
microprobe_data_means <- microprobe_data_means[c(1,3,5,8,11,14),-length(microprobe_data_means)]
microprobe_data_means <- mefa:::rep.data.frame(microprobe_data_means, each = 3)
sample_mapping_data <- cbind(sample_mapping_data[,-1], microprobe_data_means)
sample_mapping_data$Sample <- factor(sample_mapping_data$Sample)
sample_mapping_data$O2.mgL[which(sample_mapping_data$O2.mgL < 0)] <- 0

#plot chemical data
chemistry_data <- read.csv("E:/hazen_metagenome/processed_files/chemistry_data.csv", sep="\t")
chemistry_data$Depth.cm <- factor(chemistry_data$Depth.cm)
chemistry_sds <- data.frame(cbind(chemistry_data[1:2], matrix(NA, ncol=5, nrow=14)))
colnames(chemistry_sds) <- colnames(chemistry_data)
chemistry_data <- cbind(chemistry_data, orig_microprobe_data[-4])
chemistry_data$O2.mgL[which(chemistry_data$O2.mgL < 0)] <- 0
chemistry_sds <- cbind(chemistry_sds, orig_microprobe_sds)

#transform scales for the plot
chemistry_data$O2.mgL <- chemistry_data$O2.mgL*(450/10)
chemistry_data$pH <- chemistry_data$pH*(450/10)
chemistry_data$TDP.mgL <- chemistry_data$TDP.mgL*(450/60)
chemistry_data$SO42.mgL <- chemistry_data$SO42.mgL*(450/60)
chemistry_data$Cl.mgL <- chemistry_data$Cl.mgL*1000

chemistry_sds$O2.mgL <-chemistry_sds$O2.mgL*(450/10)
chemistry_sds$pH <- chemistry_sds$pH*(450/10)
chemistry_sds$TDP.mgL <- chemistry_sds$TDP.mgL*(450/60)
chemistry_sds$SO42.mgL <- chemistry_sds$SO42.mgL*(450/60)
chemistry_sds$Cl.mgL <- chemistry_sds$Cl.mgL*1000

chemistry_data_plot <- merge(melt(chemistry_data), melt(chemistry_sds), by = c("Site", "variable", "Depth.cm"), all = T)
#chemistry_data_plot$Depth.cm <- factor(sample_data_plot$Depth.cm, levels = rev(levels(sample_data_plot$Depth.cm)))
chemistry_data_plot$sdmin <- replace(chemistry_data_plot$value.x-chemistry_data_plot$value.y, chemistry_data_plot$value.x-chemistry_data_plot$value.y < 0, 0)
chemistry_data_plot$sdmax <- chemistry_data_plot$value.x+chemistry_data_plot$value.y
chemistry_data_plot$Depth.cm <- factor(chemistry_data_plot$Depth.cm, levels = rev(levels(chemistry_data_plot$Depth.cm)))

chemistry_palette=c("#d23416",
                    "#27a720",
                    "#a31faa",
                    "#00a376",
                    "#e80060",
                    "royalblue1",
                    "#46499d",
                    "#000000")

image <- ggplot(chemistry_data_plot, aes(x=Depth.cm, y=value.x, group=Site)) + geom_point(aes(fill = variable, color = variable), shape = 18, size=7, alpha = 0.5) + 
  facet_grid(Site~., scales="free_y", space="free_y") + 
  geom_errorbar(aes(ymin=sdmin, ymax=sdmax, color = variable), width = 0.5) + coord_flip() + scale_y_continuous(breaks = seq(0,450,100), limits = c(0,450)) +
  scale_color_manual(values = chemistry_palette) + ggtitle("Geochemistry") + xlab("Depth from sediment surface (cm)") + 
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.spacing = unit(2, "lines"), strip.text.y = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), 
        legend.position="bottom", plot.margin = unit(c(0.8,0.5,6.2,0.5), "lines"), plot.title = element_text(size=28, hjust = 0.5))
ggsave(file = "E:/hazen_metagenome/results/geochemistry.svg", plot=image, units="cm", width=12, height=12, scale = 2)

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
colnames(bars_matam) <- c(tax_table(phyla_matam)[,2])
phyla_matam_levels <- phyla_levels_all[phyla_levels_all %in% colnames(bars_matam)]
bars_matam <- bars_matam[,phyla_matam_levels]
bars_matam <- reshape2::melt(as.matrix(bars_matam), varnames=c("Sample", "Phylum"), value.name="Abundance")
bars_matam$Phylum <- factor(bars_matam$Phylum, levels = phyla_matam_levels)

#write out phyla abundances
phyla_percentages <- data.frame(otu_table(phyla_matam))
rownames(phyla_percentages) <- c(tax_table(phyla_matam)[,2])
phyla_percentages$rowmeans <- rowMeans(phyla_percentages)
write.csv(phyla_percentages, "E:/hazen_metagenome/phyla_abundances.csv", quote=F)

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
write.table(otu_table, "E:/hazen_metagenome/matam_tax_table_qcd.txt", quote=F, sep="\t", row.names = F)

#the following commands ran at CAC cluster
# biom convert -i ~/matti/Hazen_metagenome/matam_tax_table_qcd.txt -o ~/matti/Hazen_metagenome/matam_tax_table_qcd.biom --table-type="OTU table" --to-json
# python ~/bin/FAPROTAX_1.0/collapse_table.py -i ~/matti/Hazen_metagenome/matam_tax_table_qcd.biom -o ~/matti/Hazen_metagenome/matam_func_table.biom -g ~/bin/FAPROTAX_1.0/FAPROTAX_Hazen.txt --collapse_by_metadata 'taxonomy' --group_leftovers_as 'other' --out_group_overlaps ~/matti/Hazen_metagenome/matam_func_table_overlaps.txt --output_format_group_overlaps classical  -l ~/matti/Hazen_metagenome/matam_FAPROTAX.log --disable_group_set_operations --out_groups2records_table ~/matti/Hazen_metagenome/matam_func_table_groups.txt -v --force
# biom convert -i ~/matti/Hazen_metagenome/matam_func_table.biom -o ~/matti/Hazen_metagenome/matam_func_table.txt --to-tsv

func <- read.csv("E:/Hazen_metagenome/matam_func_table.txt",sep="\t",skip=1,row.names=1)
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

ggplot(data = NMDS_data_matam, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=factor(Site)), size = 5) + labs(color = "Sample")

ggplot(data = NMDS_data_matam_func, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=factor(Site)), size = 5) + labs(color = "Sample")


# #genome bin data and functional data analysis
# #(sample mean coverage * bin length) / average length of the reads
# #calculate relative abundance of bins and reads as % of total reads in a sample
abundances <- list(NULL)
abundances[[1]] <- read.table("E:/hazen_metagenome/processed_files/mean_coverage.txt",sep="\t",row.names=1,header=T)
abundances[[2]] <- read.table("E:/hazen_metagenome/processed_files/all_gene_coverages.txt",sep="\t",row.names=1,header=F)
colnames(abundances[[2]]) <- colnames(abundances[[1]])

genome_gene_lengths <- list(NULL)
genome_gene_lengths[[1]] <- read.csv("E:/hazen_metagenome/processed_files/bins_summary.txt",sep="\t")
genome_gene_lengths[[2]] <- read.table("E:/hazen_metagenome/processed_files/all_gene_calls.txt",sep=" ")
colnames(genome_gene_lengths[[2]]) <- c("annotation", "bin", "gene_call_id", "split_id", "start_bp", "end_bp", "gene_direction")
genome_gene_lengths[[2]] <- genome_gene_lengths[[2]][c("gene_call_id", "annotation", "bin", "split_id", "start_bp", "end_bp", "gene_direction")]
genome_gene_lengths[[2]]$total_length <- genome_gene_lengths[[2]]$end_bp-genome_gene_lengths[[2]]$start_bp

avg_read_lengths <- read.table("E:/hazen_metagenome/processed_files/sample_average_read_lengths.txt",sep="\t")
avg_read_lengths[,2] <- as.numeric(substr(avg_read_lengths[,2], 10, 15))

sample_reads <- data.frame(sample=c(1,3,5,11,14,17), reads=c(132094500, 102955607, 98825210, 129835700, 122863643, 98243245))

# rel_abundances <- list(NULL)
# for (i in 1:2) {
# rel_abundance <- abundances[[i]]
# colnames(rel_abundance) <- as.numeric(gsub("\\D", "", colnames(rel_abundance)))
# for (j in 1:6) {
#   read_hits <- function (a,b) {
#     c <- rownames(rel_abundance)[b]
#     d <- unique(genome_gene_lengths[[i]]$total_length[which(genome_gene_lengths[[i]][1] == c)])
#     a * d
#   }
#   index <- seq(1:nrow(rel_abundance))
#   rel_abundance[,j] <- mapply(read_hits, rel_abundance[,j], index)
#   read_length <- avg_read_lengths[which(avg_read_lengths[,1] == colnames(rel_abundance)[j]),2]
#   nreads <- sample_reads[which(sample_reads[,1] == colnames(rel_abundance)[j]),2]
#   rel_abundance[,j] <- ((rel_abundance[,j] / read_length) / nreads) * 100
# }
# rel_abundance <- rel_abundance[,order(as.numeric(colnames(rel_abundance)))]
# colnames(rel_abundance) <- sample_order
# rel_abundances[[i]] <- rel_abundance
# }
# 
# write.csv(rel_abundances[[1]], "E:/Hazen_metagenome/processed_files/rel_abundance_genomes.csv",quote=F)
# write.csv(rel_abundances[[2]], "E:/Hazen_metagenome/processed_files/rel_abundance_genes.csv",quote=F)

rel_abundances <- list(NULL)
rel_abundances[[1]] <- read.csv("E:/Hazen_metagenome/processed_files/rel_abundance_genomes.csv",sep=",",row.names=1)
colnames(rel_abundances[[1]]) <- sample_order
rel_abundances[[2]] <- read.csv("E:/Hazen_metagenome/processed_files/rel_abundance_genes.csv",sep=",",row.names=1)
colnames(rel_abundances[[2]]) <- sample_order


#construct phyloseq objects for (checkm) taxonomy, and metacyc/kegg pathways

phylogeny <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# checkm_taxonomy_raw <- read.csv("E:/hazen_metagenome/processed_files/tree_qa_lineage.out",sep="\t")
# checkm_taxonomy <- as.character(checkm_taxonomy_raw$Taxonomy..contained.)
# checkm_taxonomy <- gsub("[a-z]__", "", checkm_taxonomy)
# checkm_taxonomy <- strsplit(checkm_taxonomy,";")
# checkm_taxonomy <- as.matrix(rbind.fill(lapply(checkm_taxonomy,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})))
# rownames(checkm_taxonomy) <- gsub("-contigs", "", checkm_taxonomy_raw$Bin.Id)
# colnames(checkm_taxonomy) <- phylogeny
# checkm_tree <- read.tree("E:/hazen_metagenome/processed_files/tree_qa.tre")
# checkm_tree$tip.label <- gsub("-contigs", "", checkm_tree$tip.label)
new_RP_MAG_taxonomy <- read.table("E:/hazen_metagenome/processed_files/RP_MAGs_taxonomy.txt", sep="\t", header=T, na.strings=c("","NA"))
rownames(new_RP_MAG_taxonomy) <- new_RP_MAG_taxonomy$label
new_RP_MAG_tree <- as_data_frame(treeio::drop.tip(new_ggtree2, tip=new_ggtree$label[which(new_ggtree$source %in% "repository")]))
new_RP_MAG_tree <- full_join(new_RP_MAG_tree %>% select(-(5:ncol(new_RP_MAG_tree))), new_RP_MAG_taxonomy, by="label")


rel_abundance_pathways <- rel_abundances[[1]]
rel_abundance_pathways$bin <- rownames(rel_abundances[[1]])
rel_abundance_pathways <- rel_abundance_pathways[which(rel_abundance_pathways$bin %in% new_RP_MAG_tree$label),]
for (i in 1:6) {
  rel_abundance_pathways[,i] <-  sapply(rel_abundance_pathways[,i], function(x) x <- x/sum(rel_abundance_pathways[,i])*100)
}


sequencing <- read.csv("E:/hazen_metagenome/sample_sequencing_stats.csv",sep="\t")
sequencing <- sequencing[-which(sequencing$sample %in% "total"),]
sequencing <- sequencing[,-c(3,4)]
sequencing$sample <- sample_order
# sequencing$binned <- colSums(rel_abundances[[1]])*0.01*sequencing$n.reads
# sequencing$quality.bins <- colSums(rel_abundance_pathways[,-7])*0.01*sequencing$n.reads
sequencing <- reshape2::melt(sequencing)

image <- ggplot(sequencing, aes(x=factor(sample), y=value)) + geom_bar(stat="identity", color="black") + scale_x_discrete(limits=(sample_order)) + facet_wrap(~variable, scales="free") + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Sample") + ylab("Abundance (%)") + ggtitle("Sequencing statistics") + scale_y_continuous(labels = scales::comma)
ggsave(file = "E:/hazen_metagenome/results/sequecing_statistics.svg", plot=image, units="mm", width=300, height=200)

colnames_pathways <- c("bin", "accession", "pathway")

kegg_pathways <- read.table("E:/hazen_metagenome/processed_files/kegg_pathways.txt",sep="\t")
colnames(kegg_pathways) <- colnames_pathways
kegg_pathways <- merge(kegg_pathways, rel_abundance_pathways, by="bin")
kegg_pathways2 <- aggregate(kegg_pathways[,4:9], by=list(kegg_pathways$pathway), FUN=sum)
rownames(kegg_pathways2) <- kegg_pathways2$Group.1
kegg_pathways2 <- kegg_pathways2[,-1]

metacyc_pathways <- readLines("E:/hazen_metagenome/processed_files/metacyc_pathways.txt")
metacyc_pathways <- strsplit(metacyc_pathways,"\t")
metacyc_pathways <- do.call(rbind.data.frame, metacyc_pathways)[1:3]
colnames(metacyc_pathways) <- colnames_pathways
metacyc_annotation <- read.csv("E:/hazen_metagenome/processed_files/metacyc_annotation.csv", sep="\t", header=F)
colnames(metacyc_annotation) <- c("pathway", "process", "category")
metacyc_pathways <- merge(metacyc_pathways, metacyc_annotation, by="pathway")
metacyc_pathways <- merge(metacyc_pathways, rel_abundance_pathways, by="bin")
#remove "suspicious" (archaeal) pathways
metacyc_pathways <- metacyc_pathways[-which(metacyc_pathways$process %in% "Suspicious"),]
metacyc_pathways2 <- aggregate(metacyc_pathways[,6:11], by=list(metacyc_pathways$pathway), FUN=sum)
rownames(metacyc_pathways2) <- metacyc_pathways2$Group.1
metacyc_pathways2 <- metacyc_pathways2[,-1]
metacyc_tax_table <- unique(metacyc_pathways[c(2,4,5)])
metacyc_tax_table <- metacyc_tax_table[order(metacyc_tax_table$pathway),]
rownames(metacyc_tax_table) <- metacyc_tax_table$pathway
metacyc_tax_table <- metacyc_tax_table[-1]

marker_genes <- genome_gene_lengths[[2]]
marker_gene_annotation <- read.table("E:/hazen_metagenome/processed_files/individual_gene_annotations.txt", sep = "\t")
colnames(marker_gene_annotation) <- c("annotation", "gene", "process", "category")
marker_gene_annotation[] <- lapply(marker_gene_annotation, factor)
marker_gene_annotation[] <- lapply(marker_gene_annotation, f)
marker_gene_abundances <- rel_abundances[[2]]
marker_gene_abundances$gene_call_id <- rownames(marker_gene_abundances)
marker_genes2 <- merge(marker_genes, marker_gene_annotation, by="annotation")
marker_genes2 <- merge(marker_genes2, marker_gene_abundances, by="gene_call_id")
#remove hgcB (probably non-specific HMM)
marker_genes2 <- marker_genes2[-which(marker_genes2$gene %in% "hgcB"),]
annotation_abundances <- aggregate(marker_genes2[,12:17], by=list(marker_genes2$gene), FUN=sum)
rownames(annotation_abundances) <- annotation_abundances$Group.1
annotation_abundances <- annotation_abundances[-1]
annotation_tax_table <- unique(marker_genes2[c(9:11)])
annotation_tax_table <- annotation_tax_table[order(annotation_tax_table$gene),]
rownames(annotation_tax_table) <- annotation_tax_table$gene
annotation_tax_table <- annotation_tax_table[-1]

resfam_genes <- genome_gene_lengths[[2]]
resfam_gene_annotation_raw <- read.csv("E:/hazen_metagenome/processed_files/resfams_annotations.csv", sep = "\t")
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


antismash_clusters <- read.table("E:/hazen_metagenome/processed_files/antismash_clusters.txt", sep="\t")[,1:2]
colnames(antismash_clusters) <- c("bin", "cluster")
antismash_clusters$bin <- gsub("\\sc_.*", "", antismash_clusters$bin)
antismash_clusters <- merge(antismash_clusters, rel_abundance_pathways, by="bin")
antismash_clusters2 <- aggregate(antismash_clusters[,3:8], by=list(antismash_clusters$cluster), FUN=sum)
rownames(antismash_clusters2) <- antismash_clusters2$Group.1
antismash_clusters2 <- antismash_clusters2[,-1]

RP_MAG_otu_table <- rel_abundance_pathways[1:6]
RP_MAG_tax_table <- as.matrix(new_RP_MAG_taxonomy[-1])

RP_MAG_phylo <- phyloseq(otu_table(RP_MAG_otu_table, taxa_are_rows=T), tax_table(RP_MAG_tax_table), sample_data(matam), phy_tree=as.phylo(new_RP_MAG_tree))
#RP_MAG_phylo <- transform_sample_counts(RP_MAG_phylo, function(x) 100 * x/sum(x))
kegg_phylo <- phyloseq(otu_table(kegg_pathways2,taxa_are_rows=T), sample_data(matam))
metacyc_phylo <- phyloseq(otu_table(metacyc_pathways2,taxa_are_rows=T), tax_table(as.matrix(rev(metacyc_tax_table))), sample_data(matam))
antismash_phylo <- phyloseq(otu_table(antismash_clusters2,taxa_are_rows=T), sample_data(matam))
marker_phylo <- phyloseq(otu_table(annotation_abundances,taxa_are_rows=T), tax_table(as.matrix(rev(annotation_tax_table))), sample_data(matam))
resfam_phylo <- phyloseq(otu_table(resfam_annotation_abundances,taxa_are_rows=T), tax_table(as.matrix(rev(resfam_annotation_tax_table))), sample_data(matam))

RP_MAG_distance <- DPCoA(RP_MAG_phylo)
ordu_RP_MAG <- ordinate(RP_MAG_phylo, "NMDS", distance = RP_MAG_distance$RaoDis)
continuous_sample_data <- sample_data(RP_MAG_phylo)[,3:10]
RP_MAG_env_distance <- vegdist(continuous_sample_data, method="euclidean")
RP_MAG_mantel <- mantel(RP_MAG_distance$RaoDis, RP_MAG_env_distance, method="pearson", permutations=10000)
ef_RP_MAG <- envfit(ordu_RP_MAG,sample_data(RP_MAG_phylo),permu=10000)

kegg_distance <- vegdist(wisconsin(sqrt(veganifyOTU(kegg_phylo))), distance = "bray")
ordu_kegg <- ordinate(kegg_phylo, "NMDS", distance = kegg_distance)
kegg_env_distance <- vegdist(continuous_sample_data, method="euclidean")
kegg_mantel <- mantel(kegg_distance, kegg_env_distance, method="pearson", permutations=10000)
ef_kegg <- envfit(ordu_kegg,sample_data(kegg_phylo),permu=10000)

metacyc_distance <- vegdist(wisconsin(sqrt(veganifyOTU(metacyc_phylo))), distance = "bray")
ordu_metacyc <- ordinate(metacyc_phylo, "NMDS", distance = metacyc_distance)
metacyc_env_distance <- vegdist(continuous_sample_data, method="euclidean")
metacyc_mantel <- mantel(metacyc_distance, metacyc_env_distance, method="pearson", permutations=10000)
ef_metacyc <- envfit(ordu_metacyc,sample_data(metacyc_phylo),permu=10000)

cluster_palette <- c("#03ff00",
                     "#d660cd",
                     "#00e1ff",
                     "#ffb100",
                     "#e65141",
                     "#9800ff")

dataset_names <- c("matam", "matam_func", "RP_MAG", "kegg", "metacyc")
rtsne_sample_data <- sample_data(matam)
rtsne_sample_names <- sample_names(matam)
for (i in 1:length(dataset_names)) {
set.seed(42)
if (dataset_names[i] %in% "RP_MAG") {
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
ggsave(file = paste0("E:/hazen_metagenome/results/rtsne_", dataset_names[i], ".svg"), plot=image, units="mm", width=225, height=150)
}


#phyla bars for the MAG data
phyla_RP_MAG <- RP_MAG_phylo
if(any(tax_table(phyla_RP_MAG)[,2] %in% "Proteobacteria")) {
  replacement_proteobacteria <- tax_table(phyla_RP_MAG)[tax_table(phyla_RP_MAG)[,2] %in% "Proteobacteria",3]
  replacement_proteobacteria[replacement_proteobacteria %in% NA] <- "Proteobacteria"
  tax_table(phyla_RP_MAG)[tax_table(phyla_RP_MAG)[,2] %in% "Proteobacteria",2] <- replacement_proteobacteria
}
if(any(is.na(tax_table(matam)[,2]))) {
  tax_table(phyla_RP_MAG)[which(is.na(tax_table(phyla_RP_MAG)[,2])),2] <- tax_table(phyla_RP_MAG)[which(is.na(tax_table(phyla_RP_MAG)[,2])),1]
}

tax_table(phyla_RP_MAG)[tax_table(phyla_RP_MAG)[,2] %in% NA,2] <- "Unknown"
phyla_RP_MAG <- tax_glom(phyla_RP_MAG, "Phylum")

phyla_RP_MAG <- transform_sample_counts(phyla_RP_MAG, function(x) 100 * x/sum(x))
bars_RP_MAG <- data.frame(t(otu_table(phyla_RP_MAG)))

colnames(bars_RP_MAG) <- c(tax_table(phyla_RP_MAG)[,2])
phyla_RP_MAG_levels <- c(phyla_levels_all[phyla_levels_all %in% colnames(bars_RP_MAG)], "Unknown")
bars_RP_MAG <- bars_RP_MAG[,phyla_RP_MAG_levels]
bars_RP_MAG <- reshape2::melt(as.matrix(bars_RP_MAG), varnames=c("Sample", "Phylum"), value.name="Abundance")
bars_RP_MAG$Phylum <- factor(bars_RP_MAG$Phylum, levels = phyla_RP_MAG_levels)
#match the colors of phyla in the tree figure
colors_matam <- phylum_to_color_map[which(phylum_to_color_map$Phylum %in% levels(bars_matam$Phylum)),]
colors_RP_MAG <- phylum_to_color_map[which(phylum_to_color_map$Phylum %in% levels(bars_RP_MAG$Phylum)),]
#add missing levels
colors_matam <- rbind(colors_matam, data.frame(Phylum=setdiff(unique(bars_matam$Phylum), unique(colors_matam$Phylum)), color=c("#4cab98","#7aa444","darkgrey")))
colors_RP_MAG <- rbind(colors_RP_MAG, data.frame(Phylum=c("Rhodothermaeota", "Unknown"), color=c("#ca5670","lightsteelblue4")))
matchcolors_matam <- colors_matam[match(levels(bars_matam$Phylum), colors_matam$Phylum),]
matchcolors_matam <- as.vector(matchcolors_matam$color)
matchcolors_RP_MAG <- colors_RP_MAG[match(phyla_RP_MAG_levels, colors_RP_MAG$Phylum),]
matchcolors_RP_MAG <- as.vector(matchcolors_RP_MAG$color)

#plot phyla for matam and RP_MAG
image <- ggplot(bars_matam, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_matam$Sample[nchar(as.character(bars_matam$Sample))!=1]) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=matchcolors_matam, name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("MATAM taxonomy")
ggsave(file = "E:/hazen_metagenome/results/matam_taxonomy_bars.svg", plot=image, units="mm", width=400, height=200)

image <- ggplot(bars_RP_MAG, aes(x=factor(Sample), y=Abundance, fill=factor(Phylum))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_RP_MAG$Sample[nchar(as.character(bars_RP_MAG$Sample))!=1]) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=matchcolors_RP_MAG, name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("Metagenome Assembled Genome taxonomy")
ggsave(file = "E:/hazen_metagenome/results/RP_MAG_taxonomy_bars.svg", plot=image, units="mm", width=400, height=200)

#plot bars for FAPROTAX and heatmaps for kegg, metacyc, antismash, and marker genes
#bars_func <- rbind(bars_func, data.frame(Sample = empty_levels, Function = levels(bars_func$Function)[1], Abundance = 0))
image <- ggplot(bars_func, aes(x=factor(Sample), y=Abundance, fill=factor(Function))) + scale_x_discrete(limits=rev(sample_order), breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_bar(stat="identity") + scale_fill_manual(values=matchcolors_matam, name = "Phylum", guide = guide_legend(reverse = T)) + 
  theme(axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + coord_flip() + ggtitle("FAPROTAX predictions")
ggsave(file = "E:/hazen_metagenome/results/FAPROTAX_function_bars.svg", plot=image, units="mm", width=400, height=200)

kegg_heatmap <- data.frame(otu_table(kegg_phylo))
colnames(kegg_heatmap) <- sample_order
kegg_order <- readLines("E:/hazen_metagenome/processed_files/kegg_pathways.order")
kegg_levels <- kegg_order[which(kegg_order %in% rownames(kegg_heatmap))]
kegg_heatmap <- kegg_heatmap[match(kegg_levels, rownames(kegg_heatmap)),]
kegg_heatmap <- reshape2::melt(log10(as.matrix(kegg_heatmap)), varnames=c("Pathway", "Sample"), value.name="Abundance")
kegg_heatmap$Pathway <- factor(kegg_heatmap$Pathway)

metacyc_heatmap <- data.frame(otu_table(metacyc_phylo))
colnames(metacyc_heatmap) <- sample_order

metacyc_levels <- metacyc_annotation$pathway[which(metacyc_annotation$pathway %in% rownames(metacyc_heatmap))]
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

image <- ggplot(kegg_heatmap, aes(x=factor(Sample), y=factor(Pathway))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.border = element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("KEGG pathways") +labs(fill="Log 10 abundance")
ggsave(file = "E:/hazen_metagenome/results/KEGG_pathways_heatmap.svg", plot=image, units="mm", width=400, height=600)

image <- ggplot(metacyc_heatmap, aes(x=factor(Sample), y=factor(Pathway))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 16), panel.border = element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("MetaCyc pathways") +labs(fill="Log 10 abundance")
ggsave(file = "E:/hazen_metagenome/results/metacyc_pathways_heatmap.svg", plot=image, units="mm", width=400, height=1200)

image <- ggplot(antismash_heatmap, aes(x=factor(Sample), y=factor(Cluster))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.border = element_blank(), axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("Antismash clusters") +labs(fill="Log 10 abundance")
ggsave(file = "E:/hazen_metagenome/results/antismash_clusters_heatmap.svg", plot=image, units="mm", width=400, height=600)

image <- ggplot(marker_heatmap, aes(x=factor(Sample), y=factor(Gene))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.border = element_blank(), axis.text.x = element_text(size=22, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("Marker gene abundances") +labs(fill="Log 10 abundance")
ggsave(file = "E:/hazen_metagenome/results/marker_abudance_heatmap.svg", plot=image, units="mm", width=400, height=600)

image <- ggplot(resfam_heatmap, aes(x=factor(Sample), y=factor(Resfam))) + scale_x_discrete(limits=sample_order, breaks = bars_func$Sample[nchar(as.character(bars_func$Sample))!=1]) + 
  geom_tile(aes(fill = Abundance), color = "white") + scale_fill_gradient(low = "white", high = "darkslategrey") +
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.border = element_blank(), axis.text.x = element_text(size=22, angle=45, vjust=1, hjust=1), axis.ticks.x = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ggtitle("Resfams abundances") +labs(fill="Log 10 abundance")
ggsave(file = "E:/hazen_metagenome/results/resfam_abudance_heatmap.svg", plot=image, units="mm", width=400, height=600)


#genome data here

genomes_qa <- readLines("E:/hazen_metagenome/processed_files/qa_trimmed.out")
genomes_qa <- genomes_qa[-c(1:3,(length(genomes_qa)))]
genomes_qa <- data.frame(matrix(unlist(sapply(genomes_qa, function(x) strsplit(x, "\\s+"))), nrow=length(genomes_qa), byrow=T))
genome_header <- c("Bin.Id", "Marker.lineage", "UID", "n.genomes", "n.markers", "n.marker.sets", "0", "1", "2", "3", "4", "5+", "Completeness", "Contamination", "Strain.heterogeneity")
colnames(genomes_qa) <- genome_header
genomes_qa_lineage <- read.table("E:/hazen_metagenome/processed_files/tree_qa_lineage.out", sep="\t", header=T, comment.char = "")
genomes_qa <- merge(genomes_qa, genomes_qa_lineage[,-c(12:ncol(genomes_qa_lineage))], by="Bin.Id")
genomes_qa$Bin.Id <- gsub("-contigs", "", genomes_qa$Bin.Id)
genomes_qa <- genomes_qa[order(genomes_qa$Bin.Id),]
genomes_qa_over50 <- genomes_qa[which(as.numeric(as.character(genomes_qa$Completeness)) > 50),]
genomes_qa <- genomes_qa[which(genomes_qa$Bin.Id %in% taxa_names(RP_MAG_phylo)),]
genomes <- as.list(as.data.frame(t(genomes_qa[,c(13:15,21:24)])))
names(genomes) <- genomes_qa$Bin.Id

genomes2 <- lapply(1:length(genomes),function(x) {
  df1 <- as.data.frame(t(as.numeric(as.character(genomes[[x]]))))
  colnames(df1) <- names(genomes[[x]])
  binname <- names(genomes)[x]
  y <- which(rownames(new_RP_MAG_taxonomy)==binname)
  df2 <- as.data.frame(kegg_pathways[which(kegg_pathways$bin==binname),3])
  colnames(df2) <- "kegg_pathways"
  df3 <- as.data.frame(metacyc_pathways[which(metacyc_pathways$bin==binname),c(2,4,5)])
  colnames(df3) <- c("metacyc_pathway", "process", "category")
  df4 <- as.data.frame(antismash_clusters[which(antismash_clusters$bin==binname),2])
  colnames(df4) <- "antismash_clusters"
  df5 <- as.data.frame(resfam_genes2[which(resfam_genes2$bin==binname),c(9,11)])
  colnames(df5) <- c("resfam_gene", "category")
  df6 <- as.data.frame(marker_genes2[which(marker_genes2$bin==binname),9:11])
  colnames(df6) <- c("marker_gene", "process", "category")
  list(stats=df1,taxonomy=new_RP_MAG_taxonomy[y,],kegg=df2,metacyc=df3,antismash=df4,resfams=df5,marker_genes=df6)
  })
                                                 
names(genomes2) <- genomes_qa$Bin.Id

#export data for DS-FDR
# dsfdr_out <- data.frame(otu_table(checkm_phylo))
# dsfdr_out <- make_biom(dsfdr_out)
# write_biom(dsfdr_out,"E:/Hazen_metagenome/processed_files/dsfdr_data.biom")

#random forest analysis and partial dependence

# #make tables of data available for models
# variables <- colnames(sample_data(matam))[-1]
# datasets <- c(rep("RP_MAG_phylo",length(variables)), rep("kegg_phylo",length(variables)), rep("metacyc_phylo",length(variables)), rep("antismash_phylo",length(variables)), rep("marker_phylo",length(variables)))
# modelstorun_frame <- data.frame(set=datasets, model=variables)
# continuous_variables <- variables[-1]


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
#   if (modelstorun_frame$set[j] %in% c("metacyc_phylo", "marker_phylo")) {
#     predictors_phylo <- tax_glom(eval(parse(text=paste0(modelstorun_frame$set[j]))), taxrank="process")
#     predictors <- t(otu_table(predictors_phylo))
#     colnames(predictors) <- data.frame(tax_table(predictors_phylo))$process
#   } else {
#     predictors <-
#       t(otu_table(eval(parse(
#         text = paste0(modelstorun_frame$set[j])
#       ))))
#   }
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
#     ggsave(file = paste0("E:/hazen_metagenome/results/", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
#            units="mm", width=1200, height=200)
#   } else if (length(unique(taxonomy_merged$variable)) >= 10) {
#     ggsave(file = paste0("E:/hazen_metagenome/results/", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
#            units="mm", width=600, height=200)
#   } else {
#     ggsave(file = paste0("E:/hazen_metagenome/results/", modelstorun_frame$set[k], "_pd_subplot_",modelstorun_frame$model[k],".svg"), plot=pd_ggplots[[k]],
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
# write.csv(models_out,"E:/Hazen_metagenome/results/final/results_table.csv",quote=F)


all_levels <- list(kegg_levels, metacyc_annotation[,1:2], antismash_levels, data.frame(resfam_gene_annotation$gene, resfam_gene_annotation$category), data.frame(marker_gene_annotation$gene, marker_gene_annotation$process))
genome_list_2 <- list(NULL)
genome_list_3 <- list(NULL)
for (k in 1:7) {
  if (k < 3) {
    genome_list_2[[k]] <- do.call(rbind,lapply(genomes2, '[[', k))
    genome_list_3[[k]] <- do.call(rbind,lapply(genomes2, '[[', k))
  } else if (any(k == c(3,5))) {
    common_elements <- lapply(genomes2, '[[', k)
    common_elements <- Reduce(intersect, do.call(c,common_elements))
    genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(genomes2, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements <- data.frame(table(all_elements[,1]))
    all_elements <- all_elements[match(all_levels[[k-2]], as.character(all_elements$Var1)),]
    colnames(all_elements) <- c("pathway", "count")
    all_elements <- all_elements[which(all_elements$count > 0),]
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    genome_list_3[[k]] <- all_elements
  } else if (k == 4) {
    common_elements <- lapply(genomes2, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("pathway", "processes")
    genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(genomes2, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$metacyc_pathway))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-2]]$pathway), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("pathway", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,2:4])
    all_elements2 <- unique(data.frame(table(all_elements2$process)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-2]]$process)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("process", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("pathway", "processes")
    genome_list_3[[k]] <- all_elements 
    } else if (k == 6) {
    common_elements <- lapply(genomes2, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("genes", "category")
    genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(genomes2, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$resfam_gene))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-2]]$resfam_gene_annotation.gene), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("gene", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,c("category","bin")])
    all_elements2 <- unique(data.frame(table(all_elements2$category)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-2]]$resfam_gene_annotation.category)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("category", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("genes", "category")
    genome_list_3[[k]] <- all_elements
  } else {
    common_elements <- lapply(genomes2, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("genes", "processes")
    genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(genomes2, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$marker_gene))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-2]]$marker_gene_annotation.gene), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("gene", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,2:4])
    all_elements2 <- unique(data.frame(table(all_elements2$process)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-2]]$marker_gene_annotation.process)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("process", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("genes", "processes")
    genome_list_3[[k]] <- all_elements
  }
}
names(genome_list_2) <- names(genomes2[[1]])
names(genome_list_3) <- names(genomes2[[1]])

#compare to NCBI genomes at phylum level
levels_to_keep <- levels(genome_list_2$taxonomy$Phylum)
levels_to_keep <- c(levels_to_keep, "Alphaproteobacteria", "Betaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria")

ncbi_to_keep <- new_ggtree[which(new_ggtree$Phylum_grp %in% levels_to_keep),]
ncbi_to_keep <- ncbi_to_keep[which(ncbi_to_keep$source %in% "repository"),]

ncbi_kegg_pathways <- read.table("E:/hazen_metagenome/processed_files/ncbi_genomes_kegg_pathways.txt",sep="\t")
CPR_kegg_pathways <- read.table("E:/hazen_metagenome/processed_files/CPR_genomes_kegg_pathways.txt",sep="\t")
ncbi_kegg_pathways <- rbind(ncbi_kegg_pathways, CPR_kegg_pathways)
colnames(ncbi_kegg_pathways) <- colnames_pathways
ncbi_kegg_pathways <- ncbi_kegg_pathways[which(ncbi_kegg_pathways$bin %in% ncbi_to_keep$label),]


#rownames(kegg_pathways2) <- kegg_pathways2$Group.1
#kegg_pathways2 <- kegg_pathways2[,-1]

ncbi_metacyc_pathways <- readLines("E:/hazen_metagenome/processed_files/ncbi_genomes_metacyc_pathways.txt")
ncbi_metacyc_pathways <- strsplit(ncbi_metacyc_pathways,"\t")
ncbi_metacyc_pathways <- do.call(rbind.data.frame, ncbi_metacyc_pathways)[1:3]
colnames(ncbi_metacyc_pathways) <- colnames_pathways
CPR_metacyc_pathways <- readLines("E:/hazen_metagenome/processed_files/CPR_genomes_metacyc_pathways.txt")
CPR_metacyc_pathways <- strsplit(CPR_metacyc_pathways,"\t")
CPR_metacyc_pathways <- do.call(rbind.data.frame, CPR_metacyc_pathways)[1:3]
colnames(CPR_metacyc_pathways) <- colnames_pathways
ncbi_metacyc_pathways <- rbind(ncbi_metacyc_pathways, CPR_metacyc_pathways)
ncbi_metacyc_pathways <- ncbi_metacyc_pathways[which(ncbi_metacyc_pathways$bin %in% ncbi_to_keep$label),]
ncbi_metacyc_pathways <- merge(ncbi_metacyc_pathways, metacyc_annotation, by="pathway")
#remove suspicious (archaeal) pathways
ncbi_metacyc_pathways <- ncbi_metacyc_pathways[-which(ncbi_metacyc_pathways$process %in% "Suspicious"),]

#rownames(metacyc_pathways2) <- metacyc_pathways2$Group.1
#metacyc_pathways2 <- metacyc_pathways2[,-1]

ncbi_marker_genes <- read.table("E:/hazen_metagenome/processed_files/ncbi_genomes_individual_genes.txt", sep = " ")
CPR_marker_genes <- read.table("E:/hazen_metagenome/processed_files/CPR_genomes_individual_genes.txt", sep = " ")
ncbi_marker_genes <- rbind(ncbi_marker_genes, CPR_marker_genes)[,-4]
colnames(ncbi_marker_genes) <- c("annotation", "bin", "identifier")

ncbi_marker_genes <- ncbi_marker_genes[which(ncbi_marker_genes$bin %in% ncbi_to_keep$label),]
ncbi_marker_genes <- merge(ncbi_marker_genes, marker_gene_annotation, by="annotation")
#remove hgcB (too non-specific)
ncbi_marker_genes <- ncbi_marker_genes[-which(ncbi_marker_genes$gene %in% "hgcB"),]

ncbi_resfams <- read.table("E:/hazen_metagenome/processed_files/ncbi_genomes_resfams.txt", sep = " ")
CPR_resfams <- read.table("E:/hazen_metagenome/processed_files/CPR_genomes_resfams.txt", sep = " ")
ncbi_resfams <- rbind(ncbi_resfams, CPR_resfams)[,-4]
colnames(ncbi_resfams) <- c("annotation", "bin", "identifier")

ncbi_resfams <- ncbi_resfams[which(ncbi_resfams$bin %in% ncbi_to_keep$label),]
ncbi_resfams <- merge(ncbi_resfams, resfam_gene_annotation, by="annotation")

ncbi_antismash_clusters <- read.table("E:/hazen_metagenome/processed_files/ncbi_genomes_antismash_clusters.txt", sep="\t")[,1:2]
CPR_antismash_clusters <- read.table("E:/hazen_metagenome/processed_files/CPR_genomes_antismash_clusters.txt", sep="\t")[,1:2]
ncbi_antismash_clusters <- rbind(ncbi_antismash_clusters, CPR_antismash_clusters)
colnames(ncbi_antismash_clusters) <- c("bin", "cluster")
ncbi_antismash_clusters$bin <- gsub("\\s.*$", "", ncbi_antismash_clusters$bin)
ncbi_antismash_clusters <- ncbi_antismash_clusters[which(ncbi_antismash_clusters$bin %in% ncbi_to_keep$label),]
new_antismash_levels <- factor(c(as.character(antismash_levels), setdiff(as.character(levels(ncbi_antismash_clusters$cluster)), as.character(antismash_levels))))

all_levels[[3]] <- new_antismash_levels

ncbi_genomes <- lapply(1:nrow(ncbi_to_keep), function(x) {
  binname <- ncbi_to_keep$label[x]
  y <- grep(binname, as.character(all_taxonomy$label))
  df2 <- as.data.frame(ncbi_kegg_pathways[which(ncbi_kegg_pathways$bin==binname),3])
  colnames(df2) <- "kegg_pathways"
  df3 <- as.data.frame(ncbi_metacyc_pathways[which(ncbi_metacyc_pathways$bin==binname),c(1,4,5)])
  colnames(df3) <- c("metacyc_pathway", "process", "category")
  df4 <- as.data.frame(ncbi_antismash_clusters[which(ncbi_antismash_clusters$bin==binname),2])
  colnames(df4) <- "antismash_clusters"
  df5 <- as.data.frame(ncbi_resfams[which(ncbi_resfams$bin==binname),c(4,6)])
  colnames(df5) <- c("resfam_gene", "category")
  df6 <- as.data.frame(ncbi_marker_genes[which(ncbi_marker_genes$bin==binname),c(4,5)])
  colnames(df6) <- c("marker_gene", "process")
  list(taxonomy=all_taxonomy[y,1:8],kegg=df2,metacyc=df3,antismash=df4,resfams=df5,marker_genes=df6)
})

names(ncbi_genomes) <- ncbi_to_keep$label


#continue from here ########### set new levels for NCBI genomes ###########
ncbi_genome_list_2 <- list(NULL)
ncbi_genome_list_3 <- list(NULL)
for (k in 1:6) {
  if (k < 2) {
    ncbi_genome_list_2[[k]] <- do.call(rbind,lapply(ncbi_genomes, '[[', k))
    ncbi_genome_list_3[[k]] <- do.call(rbind,lapply(ncbi_genomes, '[[', k))
  } else if (any(k == c(2,4))) {
    common_elements <- lapply(ncbi_genomes, '[[', k)
    common_elements <- Reduce(intersect, do.call(c,common_elements))
    ncbi_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(ncbi_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements <- data.frame(table(all_elements[,1]))
    all_elements <- all_elements[match(all_levels[[k-1]], as.character(all_elements$Var1)),]
    colnames(all_elements) <- c("pathway", "count")
    all_elements <- all_elements[which(all_elements$count > 0),]
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    ncbi_genome_list_3[[k]] <- all_elements
  } else if (k == 3) {
    common_elements <- lapply(ncbi_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("pathway", "processes")
    ncbi_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(ncbi_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$metacyc_pathway))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-1]]$pathway), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("pathway", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,2:4])
    all_elements2 <- unique(data.frame(table(all_elements2$process)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-1]]$process)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("process", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("pathway", "processes")
    ncbi_genome_list_3[[k]] <- all_elements 
  } else if (k == 5) {
    common_elements <- lapply(ncbi_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("genes", "category")
    ncbi_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(ncbi_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$resfam_gene))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-1]]$resfam_gene_annotation.gene), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("gene", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,c("category","bin")])
    all_elements2 <- unique(data.frame(table(all_elements2$category)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-1]]$resfam_gene_annotation.category)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("category", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("genes", "category")
    ncbi_genome_list_3[[k]] <- all_elements
  } else {
    common_elements <- lapply(ncbi_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("genes", "processes")
    ncbi_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(ncbi_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$marker_gene))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-1]]$marker_gene_annotation.gene), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("gene", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,2:3])
    all_elements2 <- unique(data.frame(table(all_elements2$process)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-1]]$marker_gene_annotation.process)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("process", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("genes", "processes")
    ncbi_genome_list_3[[k]] <- all_elements
  }
}
names(ncbi_genome_list_2) <- names(ncbi_genomes[[1]])
names(ncbi_genome_list_3) <- names(ncbi_genomes[[1]])

longbranch_genomes <- genomes2[which(names(genomes2) %in% over_the_line_MAG$label)]

longbranch_genome_list_2 <- list(NULL)
longbranch_genome_list_3 <- list(NULL)
for (k in 1:7) {
  if (k < 3) {
    longbranch_genome_list_2[[k]] <- do.call(rbind,lapply(longbranch_genomes, '[[', k))
    longbranch_genome_list_3[[k]] <- do.call(rbind,lapply(longbranch_genomes, '[[', k))
  } else if (any(k == c(3,5))) {
    common_elements <- lapply(longbranch_genomes, '[[', k)
    common_elements <- Reduce(intersect, do.call(c,common_elements))
    longbranch_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(longbranch_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements <- data.frame(table(all_elements[,1]))
    all_elements <- all_elements[match(all_levels[[k-2]], as.character(all_elements$Var1)),]
    colnames(all_elements) <- c("pathway", "count")
    all_elements <- all_elements[which(all_elements$count > 0),]
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    longbranch_genome_list_3[[k]] <- all_elements
  } else if (k == 4) {
    common_elements <- lapply(longbranch_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("pathway", "processes")
    longbranch_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(longbranch_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$metacyc_pathway))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-2]]$pathway), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("pathway", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,2:4])
    all_elements2 <- unique(data.frame(table(all_elements2$process)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-2]]$process)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("process", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("pathway", "processes")
    longbranch_genome_list_3[[k]] <- all_elements 
  } else if (k == 6) {
    common_elements <- lapply(longbranch_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("genes", "category")
    longbranch_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(longbranch_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$resfam_gene))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-2]]$resfam_gene_annotation.gene), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("gene", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,c("category","bin")])
    all_elements2 <- unique(data.frame(table(all_elements2$category)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-2]]$resfam_gene_annotation.category)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("category", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("genes", "category")
    longbranch_genome_list_3[[k]] <- all_elements
  } else {
    common_elements <- lapply(longbranch_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("genes", "processes")
    longbranch_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(longbranch_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$marker_gene))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-2]]$marker_gene_annotation.gene), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("gene", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,2:4])
    all_elements2 <- unique(data.frame(table(all_elements2$process)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-2]]$marker_gene_annotation.process)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("process", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("genes", "processes")
    longbranch_genome_list_3[[k]] <- all_elements
  }
}
names(longbranch_genome_list_2) <- names(longbranch_genomes[[1]])
names(longbranch_genome_list_3) <- names(longbranch_genomes[[1]])

longbranch_ncbi_genomes <- ncbi_genomes[which(names(ncbi_genomes) %in% over_the_line_ncbi$label)]

longbranch_ncbi_genome_list_2 <- list(NULL)
longbranch_ncbi_genome_list_3 <- list(NULL)
for (k in 1:6) {
  if (k < 2) {
    longbranch_ncbi_genome_list_2[[k]] <- do.call(rbind,lapply(longbranch_ncbi_genomes, '[[', k))
    longbranch_ncbi_genome_list_3[[k]] <- do.call(rbind,lapply(longbranch_ncbi_genomes, '[[', k))
  } else if (any(k == c(2,4))) {
    common_elements <- lapply(longbranch_ncbi_genomes, '[[', k)
    common_elements <- Reduce(intersect, do.call(c,common_elements))
    longbranch_ncbi_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(longbranch_ncbi_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements <- data.frame(table(all_elements[,1]))
    all_elements <- all_elements[match(all_levels[[k-1]], as.character(all_elements$Var1)),]
    colnames(all_elements) <- c("pathway", "count")
    all_elements <- all_elements[which(all_elements$count > 0),]
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    longbranch_ncbi_genome_list_3[[k]] <- all_elements
  } else if (k == 3) {
    common_elements <- lapply(longbranch_ncbi_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("pathway", "processes")
    longbranch_ncbi_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(longbranch_ncbi_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$metacyc_pathway))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-1]]$pathway), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("pathway", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,2:4])
    all_elements2 <- unique(data.frame(table(all_elements2$process)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-1]]$process)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("process", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("pathway", "processes")
    longbranch_ncbi_genome_list_3[[k]] <- all_elements 
  } else if (k == 5) {
    common_elements <- lapply(longbranch_ncbi_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("genes", "category")
    longbranch_ncbi_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(longbranch_ncbi_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$resfam_gene))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-1]]$resfam_gene_annotation.gene), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("gene", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,c("category","bin")])
    all_elements2 <- unique(data.frame(table(all_elements2$category)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-1]]$resfam_gene_annotation.category)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("category", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("genes", "category")
    longbranch_ncbi_genome_list_3[[k]] <- all_elements
  } else {
    common_elements <- lapply(longbranch_ncbi_genomes, '[[', k)
    common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
    common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
    common_elements <- list(common_elements1, common_elements2)
    names(common_elements) <- c("genes", "processes")
    longbranch_ncbi_genome_list_2[[k]] <- common_elements
    all_elements <- data.frame(do.call(rbind,lapply(longbranch_ncbi_genomes, '[[', k)))
    all_elements$bin <- sub("\\..*", "", row.names(all_elements))
    all_elements <- unique(all_elements)
    all_elements1 <- data.frame(table(all_elements$marker_gene))
    all_elements1 <- unique(all_elements1[match(as.character(all_levels[[k-1]]$marker_gene_annotation.gene), as.character(all_elements1$Var1)),])
    colnames(all_elements1) <- c("gene", "count")
    all_elements1 <- all_elements1[which(all_elements1$count > 0),]
    all_elements2 <- unique(all_elements[,2:3])
    all_elements2 <- unique(data.frame(table(all_elements2$process)))
    all_elements2 <- unique(all_elements2[match(unique(as.character(all_levels[[k-1]]$marker_gene_annotation.process)), as.character(all_elements2$Var1)),])
    colnames(all_elements2) <- c("process", "count")
    all_elements2 <- all_elements2[which(all_elements2$count > 0),]
    all_elements <- list(all_elements1, all_elements2)
    if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
    names(all_elements) <- c("genes", "processes")
    longbranch_ncbi_genome_list_3[[k]] <- all_elements
  }
}
names(longbranch_ncbi_genome_list_2) <- names(longbranch_ncbi_genomes[[1]])
names(longbranch_ncbi_genome_list_3) <- names(longbranch_ncbi_genomes[[1]])

metacyc_processes <- merge(ncbi_genome_list_3$metacyc$processes, genome_list_3$metacyc$processes, by="process", all=T)
metacyc_processes <- merge(metacyc_processes, longbranch_ncbi_genome_list_3$metacyc$processes, by="process", all=T)
colnames(metacyc_processes) <- c("process", "NCBI", "MAGs", "NCBI_long_branching")
metacyc_processes <- merge(metacyc_processes, longbranch_genome_list_3$metacyc$processes, by="process", all=T)
colnames(metacyc_processes)[5] <- "MAG_long_branching"
metacyc_processes <- metacyc_processes[match(as.character(ncbi_genome_list_3$metacyc$processes$process), metacyc_processes$process),]

marker_processes <- merge(ncbi_genome_list_3$marker_genes$processes, genome_list_3$marker_genes$processes, by="process", all=T)
marker_processes <- merge(marker_processes, longbranch_ncbi_genome_list_3$marker_genes$processes, by="process", all=T)
colnames(marker_processes) <- c("process", "NCBI", "MAGs", "NCBI_long_branching")
marker_processes <- merge(marker_processes, longbranch_genome_list_3$marker_genes$processes, by="process", all=T)
colnames(marker_processes)[5] <- "MAG_long_branching"
marker_processes <- marker_processes[match(as.character(ncbi_genome_list_3$marker_genes$processes$process), marker_processes$process),]

all_processes <- rbind(metacyc_processes, marker_processes)
all_processes$process <- gsub("_", " ", all_processes$process)
colnames(all_processes) <- gsub("_", " ", colnames(all_processes))

process_levels <- rev(c("Cell structure biosynthesis","Fatty acid and lipid biosynthesis","Amino acid biosynthesis","Aromatic compound biosynthesis",
  "Secondary metabolite biosynthesis","Antibiotic biosynthesis","Other biosynthesis","Nitrogen assimilation",
  "Selenate reduction","Amino acid degradation","Carbohydrate degradation","Aromatic compound degradation",
  "Secondary metabolite degradation","Chlorinated compound degradation","Glycerol degradation","DMSO degradation",
  "Phosphate degradation","Organic matter degradation","Hydrogen production","Phototrophy","Autotrophic CO2 assimilation",
  "Fermentation","Energy generation","Energy conservation","Methylotrophic metabolism","Metabolic pathway",
  "N fixation","Nitrification bacteria","Denitrification","DNRA","DNRA Polysulfide reduction",
  "Sulfur compound reduction","Sulfite reduction to sulfide reversible", "Sulfite reduction to sulfide",
  "Sulfur compound oxidation","Sulfide oxidation","Sulfite oxidation to sulfate",
  "Phosphorus uptake","Polyphosphate synthesis","Mercury resistance","Mercury methylation","Arsenic resistance",
  "Pb Zn resistance or homeostasis","Cd Pb resistance","Cd Co Cu Pb Zn resistance or homeostasis",
  "Cu resistance","Multimetal resistance regulation"))

all_processes$process <- factor(all_processes$process, levels=process_levels)
all_processes <- all_processes[order(all_processes$process),]
all_processes[,2:5] <- sweep(all_processes[,2:5], 2, c(2486, 55,53,18), FUN = '/')
all_processes[,2:5] <- round(all_processes[,2:5]*100, 0)
all_processes$category <- c(rep("Metal homeostasis", 8), rep("Biogeochemical cycles", 13), rep("Energy metabolism", 8), rep("Degradation", 9), rep("Biosynthesis", 9))
all_processes$category <- factor(all_processes$category, levels=rev(unique(all_processes$category)))

#chi-squared test and plot a map
all_processes_chi <- all_processes[2:5]
rownames(all_processes_chi) <- all_processes$process
all_processes_chi[is.na(all_processes_chi)] <- 0
cisq_test_all_processes <- chisq.test(all_processes_chi)
contrib_all_processes <- 100*cisq_test_all_processes$residuals^2/cisq_test_all_processes$statistic
#subset to processes that contribute more than equally (100 % / 47 processes) to the total chi-squared score
most_important_processes <- cisq_test_all_processes$residuals[rowSums(round(contrib_all_processes,3)) > 2.13,]
levels(rownames(most_important_processes)) <- process_levels
most_important_processes <- melt(most_important_processes)
image <- ggplot(data = most_important_processes, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color="white") + coord_flip() + 
  scale_fill_gradient2("Pearson residual", low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.y = element_text(size = 16), strip.text.y = element_text(angle=0, hjust=0), axis.text.x = element_text(angle=45, hjust=1), panel.border = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=28, hjust = 0.5))
ggsave(file = "E:/hazen_metagenome/results/chi_squared_heatmap.svg", plot=image, units="cm", width=10, height=12, scale=2.5)

all_processes_melt <- melt(all_processes)
all_processes_melt$variable <- factor(all_processes_melt$variable, levels=rev(colnames(all_processes)[2:5]))


image <- ggplot(all_processes_melt, aes(y=value, x=process, fill=variable)) + geom_col(position = position_dodge2(width = 0.75, preserve = "single"))  +#geom_point(aes(size=value), color="steelblue", shape=1, stroke=3) +
  coord_flip() + scale_fill_brewer(name="Collection",palette="Paired", guide = guide_legend(reverse = TRUE)) + 
  facet_grid(category~., scales="free", space="free") + ylab("Present in % of genomes") +
  theme(axis.text.y = element_text(size = 16), strip.text.y = element_text(angle=0, hjust=0), panel.border = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size=28, hjust = 0.5),
        legend.position = "bottom")
ggsave(file = "E:/hazen_metagenome/results/all_processes_bars.svg", plot=image, units="cm", width=10, height=12, scale=3)

#which MAGs take part in which processes
MAG_metacycs <- lapply(genomes2, '[[', 4)
MAG_markers <- lapply(genomes2, '[[', 7)
processes_of_interest <- c("sulfate reduction V", "sulfate activation for sulfonation", "sulfoacetaldehyde degradation I", "ammonia assimilation", "N_fixation", "Denitrification", "Nitrification_bacteria", 
                           "Sulfite_reduction_to_sulfide", "Sulfide_oxidation", "Sulfite_oxidation_to_sulfate", "Mercury_resistance",  "DNRA_Polysulfide_reduction")
process_tax_list <- list(NULL)
for (i in 1:4) {
  process_tax_list[[i]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl(processes_of_interest[i],x$metacyc_pathway)), MAG_metacycs))], '[[', 2))
}
for (i in 5:length(processes_of_interest)) {
  process_tax_list[[i]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl(processes_of_interest[i],x$process)), MAG_markers))], '[[', 2))
}
process_tax_list[[length(process_tax_list)+1]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl("nrfH",x$marker_gene)), MAG_markers))], '[[', 2))
process_tax_list[[length(process_tax_list)+1]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl("nosZ",x$marker_gene)), MAG_markers))], '[[', 2))
process_tax_list[[length(process_tax_list)+1]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl("nirK",x$marker_gene)), MAG_markers))], '[[', 2))
names(process_tax_list) <- c(processes_of_interest, "DNRA", "Nitrous oxide reduction", "Nitrite reduction to nitric oxide")
process_tax_list <- bind_rows(process_tax_list, .id="process")
process_tax_list <- merge(process_tax_list, data.frame(otu_table(RP_MAG_phylo)), by.x="label", by.y="row.names", all.x=T)
process_tax_list$mean <- rowMeans(as.matrix(process_tax_list[10:15]))
process_tax_list <- process_tax_list[-c(10:15)]
process_tax_list <- process_tax_list[order(process_tax_list$process, process_tax_list$mean, decreasing =T),]
write.csv(process_tax_list, "E:/hazen_metagenome/results/processes_of_interest_taxa.csv", quote=F, row.names = F)
process_phylum_list <- process_tax_list
process_phylum_list$Phylum <- as.character(process_phylum_list$Phylum)
process_phylum_list$Phylum[which(process_phylum_list$Phylum %in% "Proteobacteria")] <- as.character(process_phylum_list$Class[which(process_phylum_list$Phylum %in% "Proteobacteria")])
process_phylum_list <- aggregate(mean~process+Phylum, process_phylum_list, mean)
process_phylum_list <- process_phylum_list[order(process_phylum_list$process, process_phylum_list$mean, decreasing = T),]
write.csv(process_phylum_list, "E:/hazen_metagenome/results/processes_of_interest_phyla.csv", quote=F, row.names = F)

# checkm_model_shared <- list(NULL)
# checkm_model_present <- list(NULL)
# for (i in 1:length(variables)){
#   bins <- data.frame(bin=names(importance(models_2[[i]])), value=rep(NA, length(names(importance(models_2[[i]])))))
#   if (i == 1) {
#     bin_pd_data <- plot_pd_data[[i]][seq(1:(nrow(bins)*25)),]
#   } else {
#     bin_pd_data <- plot_pd_data[[i]]
#       }
#   for (j in 1:nrow(bins)){
#     if (i == 1) {
#       bin_dots_y <- bin_pd_data$prediction[which(as.character(bin_pd_data$variable) == as.character(bins[j,1]))]
#     } else {
#         bin_dots_y <- bin_pd_data$response[which(as.character(bin_pd_data$variable) == as.character(bins[j,1]))]
#         }
#     if (bin_dots_y[1] < bin_dots_y[length(bin_dots_y)]) {bins[j,2] <- "positive"} else {bins[j,2] <- "negative"}
#   }
#   positives <- as.character(bins$bin[which(bins$value %in% "positive")])
#   positives <- genomes2[positives]
#   positives2 <- list(NULL)
#   positives3 <- list(NULL)
#   for (k in 1:7) {
#     if (k < 3) {
#       positives2[[k]] <- do.call(rbind,lapply(positives, '[[', k))
#       positives3[[k]] <- do.call(rbind,lapply(positives, '[[', k))
#       }
#     else if (k < 6) {
#       common_elements <- lapply(positives, '[[', k)
#       common_elements <- Reduce(intersect, do.call(c,common_elements))
#       positives2[[k]] <- common_elements
#       all_elements <- data.frame(table(do.call(rbind,lapply(positives, '[[', k))))
#       all_elements <- all_elements[match(all_levels[[k-2]], as.character(all_elements$Var1)),]
#       colnames(all_elements) <- c("pathway", "count")
#       all_elements <- all_elements[which(all_elements$count > 0),]
#       if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
#       positives3[[k]] <- all_elements
#     } else if (k == 6) {
#       common_elements <- lapply(positives, '[[', k)
#       common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
#       common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
#       common_elements <- list(common_elements1, common_elements2)
#       names(common_elements) <- c("genes", "category")
#       positives2[[k]] <- common_elements
#       all_elements <- data.frame(do.call(rbind,lapply(positives, '[[', k)))
#       all_elements1 <- data.frame(table(all_elements$resfams))
#       all_elements1 <- all_elements1[match(as.character(all_levels[[k-2]]$resfam_gene_annotation.gene), as.character(all_elements1$Var1)),]
#       colnames(all_elements1) <- c("gene", "count")
#       all_elements1 <- all_elements1[which(all_elements1$count > 0),]
#       all_elements2 <- unique(data.frame(table(all_elements$category)))
#       all_elements2 <- all_elements2[match(unique(as.character(all_levels[[k-2]]$resfam_gene_annotation.category)), as.character(all_elements2$Var1)),]
#       colnames(all_elements2) <- c("category", "count")
#       all_elements2 <- all_elements2[which(all_elements2$count > 0),]
#       all_elements <- list(all_elements1, all_elements2)
#       if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
#       names(all_elements) <- c("genes", "category")
#       positives3[[k]] <- all_elements
#     } 
#     else {
#       common_elements <- lapply(positives, '[[', k)
#       common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
#       common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
#       common_elements <- list(common_elements1, common_elements2)
#       names(common_elements) <- c("genes", "processes")
#       positives2[[k]] <- common_elements
#       all_elements <- data.frame(do.call(rbind,lapply(positives, '[[', k)))
#       all_elements1 <- data.frame(table(all_elements$marker_gene))
#       all_elements1 <- all_elements1[match(as.character(all_levels[[k-2]]$marker_gene_annotation.gene), as.character(all_elements1$Var1)),]
#       colnames(all_elements1) <- c("gene", "count")
#       all_elements1 <- all_elements1[which(all_elements1$count > 0),]
#       all_elements2 <- unique(data.frame(table(all_elements$process)))
#       all_elements2 <- all_elements2[match(unique(as.character(all_levels[[k-2]]$marker_gene_annotation.process)), as.character(all_elements2$Var1)),]
#       colnames(all_elements2) <- c("process", "count")
#       all_elements2 <- all_elements2[which(all_elements2$count > 0),]
#       all_elements <- list(all_elements1, all_elements2)
#       if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
#       names(all_elements) <- c("genes", "processes")
#       positives3[[k]] <- all_elements
#     }
#   }
#   names(positives2) <- names(genomes2[[1]])
#   names(positives3) <- names(genomes2[[1]])
#   negatives <- as.character(bins$bin[which(bins$value %in% "negative")])
#   negatives <- genomes2[negatives]
#   negatives2 <- list(NULL)
#   negatives3 <- list(NULL)
#   for (k in 1:6) {
#     if (k < 3) {
#       negatives2[[k]] <- do.call(rbind,lapply(negatives, '[[', k))
#       negatives3[[k]] <- do.call(rbind,lapply(negatives, '[[', k))
#       }
#     else if (k < 6) {
#       common_elements <- lapply(negatives, '[[', k)
#       common_elements <- Reduce(intersect, do.call(c,common_elements))
#       negatives2[[k]] <- common_elements
#       all_elements <- data.frame(table(do.call(rbind,lapply(negatives, '[[', k))))
#       all_elements <- all_elements[match(all_levels[[k-2]], as.character(all_elements$Var1)),]
#       colnames(all_elements) <- c("pathway", "count")
#       all_elements <- all_elements[which(all_elements$count > 0),]
#       if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
#       negatives3[[k]] <- all_elements
#     } else {
#       common_elements <- lapply(negatives, '[[', k)
#       common_elements1 <- Reduce(intersect, lapply(common_elements, '[[', 1))
#       common_elements2 <- Reduce(intersect, lapply(common_elements, '[[', 2))
#       common_elements <- list(common_elements1, common_elements2)
#       names(common_elements) <- c("genes", "processes")
#       negatives2[[k]] <- common_elements
#       all_elements <- data.frame(do.call(rbind,lapply(negatives, '[[', k)))
#       all_elements1 <- data.frame(table(all_elements$marker_gene))
#       all_elements1 <- all_elements1[match(as.character(all_levels[[k-2]]$marker_gene_annotation.gene), as.character(all_elements1$Var1)),]
#       colnames(all_elements1) <- c("gene", "count")
#       all_elements1 <- all_elements1[which(all_elements1$count > 0),]
#       all_elements2 <- unique(data.frame(table(all_elements$process)))
#       all_elements2 <- all_elements2[match(unique(as.character(all_levels[[k-2]]$marker_gene_annotation.process)), as.character(all_elements2$Var1)),]
#       colnames(all_elements2) <- c("process", "count")
#       all_elements2 <- all_elements2[which(all_elements2$count > 0),]
#       all_elements <- list(all_elements1, all_elements2)
#       if (length(rownames(all_elements)) > 0)  {rownames(all_elements) <- c(1:nrow(all_elements))}
#       names(all_elements) <- c("genes", "processes")
#       negatives3[[k]] <- all_elements
#     }
#   }
#   names(negatives2) <- names(genomes2[[1]])
#   result_list1 <- list(positives2, negatives2)
#   names(result_list1) <- c("positives", "negatives")
#   checkm_model_shared[[i]] <- result_list1
#   names(checkm_model_shared)[i] <- as.character(modelstorun_frame$model[i])
#   names(negatives3) <- names(genomes2[[1]])
#   result_list2 <- list(positives3, negatives3)
#   names(result_list2) <- c("positives", "negatives")
#   checkm_model_present[[i]] <- result_list2
#   names(checkm_model_present)[i] <- as.character(modelstorun_frame$model[i])
# }

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
  top10_names[[i]] <- names(sort(taxa_sums(top10), decreasing = TRUE)[1:10])
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

ggplot(hgc_splits, aes(xmin = start_bp, xmax = end_bp, y = split_id, fill = gene, forward = direction)) + geom_gene_arrow() + 
  facet_wrap(~ split_id, scales="free") +
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
#   page <- readLines(paste0("E:/hazen_metagenome/", APD_categories[i], ".html"))
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
# candidateAP <- readFASTA("E:/hazen_metagenome/PROKKA_03022018_AP_cand2.faa")
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
