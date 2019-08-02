if (!requireNamespace("BiocManager")) {
  install.packages("BiocManager")
}

packages <- c("ggplot2", "ranger", "svglite", "ggthemes", "phyloseq", "vegan", "Rtsne", "dbscan", "mefa", "caret", "data.table", "plyr", "protr", 
              "data.table", "MLPUGS", "compiler", "DESeq2", "ape", "gggenes", "car", "edarf", "biomformat", "limma", "seqinr", "taxize", "ggtree", 
              "dplyr", "tidytree", "robust")


is.installed <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, ask=F, version= "3.8")
  }
  sapply(pkg, require, character.only = TRUE)
}
is.installed(packages)

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
too_short_seq_bins <- which(nchar(as.character(outfile)) -nchar( gsub("-", "", outfile, fixed=T)) > alignment_length*0.25)
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
new_ggtree <- as_tibble(phy_tree(concatenated_tree))
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
write.tree(new_ggtree2@phylo, file = "E:/hazen_metagenome/results/Full_tree_with_supports.nwk")

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

image <- ggtree(new_ggtree2, aes(color=group), ladderize = T, size=0.1) + geom_text2(aes(subset=!isTip, label=round(as.numeric(label),2)), hjust=-.3, size=0.3) + 
  geom_tiplab(aes(label=new_ggtree$fulltax), size=0.3, hjust=.1) + 
  geom_tippoint(aes(alpha=new_ggtree$source), color="red", shape=20, size=0.2) + 
  scale_alpha_manual(values=c("1", "0")) +
  geom_treescale(linesize=0.2, fontsize=1, offset=1) +
  scale_color_manual(values=c(treepalette2)) + ggplot2::ylim(0, 5997) + ggplot2::theme(plot.margin=unit(c(0,0,0,0),"mm"))

ggsave(file = "E:/hazen_metagenome/results/Figure_x.pdf", plot=image, units="cm", width=15, height=75, limitsize = F, scale=2.5, device=cairo_pdf)

image <- ggtree(new_ggtree2, layout="circular", aes(color=group), ladderize = T, size=0.1) + 
  #geom_nodelab(aes(subset=(!isTip & as.numeric(label) > 0.9)), label="\u25CF", size=0.2, hjust=0, color="black") + 
  geom_tippoint(aes(alpha=new_ggtree$source), color="red", shape=18, size=0.2) + 
  scale_alpha_manual(values=c("1", "0")) +
  geom_treescale(linesize=0.2, fontsize=1, offset=1, x=-2) +
  scale_color_manual(values=c(treepalette2))

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

#chi-squared test
longbranching_chisq <- data.frame(row.names=c("short", "long"), NCBI=c(3656,246), MAG=c(37,18))
longbranching_chisq_test <- chisq.test(longbranching_chisq, simulate.p.value = T, B = 1000)

#calculate pairwise distances of the best NCBI matchest for the 2 unknown bins and plot a tree
LH_MA_65_9_alignment <- read.alignment("E:/hazen_metagenome/processed_files/arb-silva.de_LH_MA_65_9_review.fasta", format="fasta")
LH_MA_65_9_distance <- dist.alignment(LH_MA_65_9_alignment, matrix="identity", gap=0)
LH_MA_65_9_distance2 <- 1-(as.matrix(LH_MA_65_9_distance)[,1]^2)
#edit names to retain only organismal identifiers
names(LH_MA_65_9_distance2) <- gsub("_16S_ribosomal_RNA_gene&g_partial_sequence$", "", names(LH_MA_65_9_distance2))
names(LH_MA_65_9_distance2) <- gsub("_partial_16S_rRNA_gene&g", "", names(LH_MA_65_9_distance2))
names(LH_MA_65_9_distance2) <- gsub("_16S_ribosomal_RNA&g_partial_sequence&g", "", names(LH_MA_65_9_distance2))
LH_MA_65_9_tree <- read.tree("E:/hazen_metagenome/processed_files/arb-silva.de_LH_MA_65_9_review.tree")
#root on Acetomicrobium hydrogeniformans ATCC BAA-1850 strain OS1
LH_MA_65_9_tree  <- root(LH_MA_65_9_tree, outgroup = 22, resolve.root=T)
LH_MA_65_9_ggtree <- as_tibble(phy_tree(LH_MA_65_9_tree))
LH_MA_65_9_ggtree$label <- gsub("_16S_ribosomal_RNA_gene&g_partial_sequence$", "", LH_MA_65_9_ggtree$label)
LH_MA_65_9_ggtree$label <- gsub("_partial_16S_rRNA_gene&g", "", LH_MA_65_9_ggtree$label)
LH_MA_65_9_ggtree$label <- gsub("_16S_ribosomal_RNA&g_partial_sequence", "", LH_MA_65_9_ggtree$label)
LH_MA_65_9_ggtree$group <- NA
LH_MA_65_9_ggtree$group[c(2,3,4,5,20,21,22,38)] <- "classified"
LH_MA_65_9_ggtree$group[c(45)] <- "MAG"
LH_MA_65_9_ggtree$group[is.na(LH_MA_65_9_ggtree$group)] <- "unknown"
LH_MA_65_9_ggtree$identity <- NA 
LH_MA_65_9_ggtree$identity[1:length(LH_MA_65_9_distance2)] <- LH_MA_65_9_distance2[!is.na(match(LH_MA_65_9_ggtree$label, names(LH_MA_65_9_distance2)))]
LH_MA_65_9_ggtree2 <- as.treedata(LH_MA_65_9_ggtree)
image <- ggtree(LH_MA_65_9_ggtree2, aes(color=group), ladderize = T, size=0.5) + geom_point2(aes(subset=(isTip & as.numeric(identity>0.9))), shape=18, size=3) + 
  geom_rootedge(0.03) + geom_nodelab(aes(subset=(!isTip & as.numeric(label) > 0.75)), label="\u25CF", size=4, color="black", hjust=1) + 
  geom_treescale(x=1, linesize=1, fontsize=5, offset=1) + geom_tiplab(hjust=-0.02) + scale_color_manual(values=c("darkgreen", "red", "black"))
ggsave(file = "E:/hazen_metagenome/results/LH_MA_65_9_tree_review.svg", plot=image, units="cm", width=10, height=10, limitsize = F, scale=3)

LH_MA_57_9_alignment <- read.alignment("E:/hazen_metagenome/processed_files/arb-silva.de_LH_MA_57_9_review.fasta", format="fasta")
LH_MA_57_9_distance <- dist.alignment(LH_MA_57_9_alignment, matrix="identity", gap=0)
LH_MA_57_9_distance2 <- 1-(as.matrix(LH_MA_57_9_distance)[,1]^2)
#edit names to retain only organismal identifiers
names(LH_MA_57_9_distance2) <- gsub("&a", ":", names(LH_MA_57_9_distance2))
names(LH_MA_57_9_distance2) <- gsub("_16S_ribosomal_RNA_gene&g_partial_sequence$", "", names(LH_MA_57_9_distance2))
names(LH_MA_57_9_distance2) <- gsub("_RNA_for_16S_rRNA&g_partial_sequence&g", "", names(LH_MA_57_9_distance2))
names(LH_MA_57_9_distance2) <- gsub("_gene_for_16S_ribosomal_RNA&g_partial_sequence&g", "", names(LH_MA_57_9_distance2))
names(LH_MA_57_9_distance2) <- gsub("_gene_for_16S_rRNA&g_partial_sequence&g", "", names(LH_MA_57_9_distance2))
names(LH_MA_57_9_distance2) <- gsub("_partial_16S_rRNA_gene&g", "", names(LH_MA_57_9_distance2))
names(LH_MA_57_9_distance2) <- gsub("_16S_ribosomal_RNA&g_partial_sequence", "", names(LH_MA_57_9_distance2))
LH_MA_57_9_tree <- read.tree("E:/hazen_metagenome/processed_files/arb-silva.de_LH_MA_57_9_review.tree")
#root on Acetomicrobium hydrogeniformans ATCC BAA-1850 strain OS1
LH_MA_57_9_tree  <- root(LH_MA_57_9_tree, outgroup = 28, resolve.root=T)
LH_MA_57_9_ggtree <- as_tibble(phy_tree(LH_MA_57_9_tree))
LH_MA_57_9_ggtree$label <- gsub("&a", ":", LH_MA_57_9_ggtree$label)
LH_MA_57_9_ggtree$label <- gsub("_16S_ribosomal_RNA_gene&g_partial_sequence$", "", LH_MA_57_9_ggtree$label)
LH_MA_57_9_ggtree$label <- gsub("_RNA_for_16S_rRNA&g_partial_sequence&g", "", LH_MA_57_9_ggtree$label)
LH_MA_57_9_ggtree$label <- gsub("_gene_for_16S_ribosomal_RNA&g_partial_sequence&g", "", LH_MA_57_9_ggtree$label)
LH_MA_57_9_ggtree$label <- gsub("_gene_for_16S_rRNA&g_partial_sequence&g", "", LH_MA_57_9_ggtree$label)
LH_MA_57_9_ggtree$label <- gsub("_partial_16S_rRNA_gene&g", "", LH_MA_57_9_ggtree$label)
LH_MA_57_9_ggtree$label <- gsub("_16S_ribosomal_RNA&g_partial_sequence", "", LH_MA_57_9_ggtree$label)
LH_MA_57_9_ggtree$group <- NA
LH_MA_57_9_ggtree$group[c(15,24,26,27,28)] <- "classified"
LH_MA_57_9_ggtree$group[c(29)] <- "MAG"
LH_MA_57_9_ggtree$group[is.na(LH_MA_57_9_ggtree$group)] <- "unknown"
LH_MA_57_9_ggtree$identity <- NA 
LH_MA_57_9_ggtree$identity[1:length(LH_MA_57_9_distance2)] <- LH_MA_57_9_distance2[!is.na(match(LH_MA_57_9_ggtree$label, names(LH_MA_57_9_distance2)))]
LH_MA_57_9_ggtree2 <- as.treedata(LH_MA_57_9_ggtree)
image <- ggtree(LH_MA_57_9_ggtree2, aes(color=group), ladderize = T, size=0.5) + geom_point2(aes(subset=(isTip & as.numeric(identity>0.9))), shape=18, size=3) +
  geom_rootedge(0.03) + geom_nodelab(aes(subset=(!isTip & as.numeric(label) > 0.75)), label="\u25CF", size=4, color="black", hjust=1) + 
  geom_treescale(x=0.6, linesize=1, fontsize=5, offset=1) + geom_tiplab(hjust=-0.02) + scale_color_manual(values=c("darkgreen", "red", "black"))
ggsave(file = "E:/hazen_metagenome/results/LH_MA_57_9_tree_review.svg", plot=image, units="cm", width=10, height=10, limitsize = F, scale=3)

# #read in the taxonomic matam data
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
chemistry_data_plot$value.x[is.na(chemistry_data_plot$value.x)] <- 0

image <- ggplot(chemistry_data_plot, aes(x=Depth.cm, y=value.x, group=Site)) + geom_point(aes(fill = variable, color = variable), shape = 18, size=7, alpha = 0.5) + 
  facet_grid(Site~., scales="free_y", space="free_y") + geom_path(aes(color = variable, group = variable)) +
  geom_errorbar(aes(ymin=sdmin, ymax=sdmax, color = variable), width = 0.5) + coord_flip() + scale_y_continuous(breaks = seq(0,450,100), limits = c(0,450)) +
  scale_color_manual(values = chemistry_palette) + ggtitle("Physicochemical variability") + xlab("Depth from sediment surface (cm)") + 
  theme(axis.text.y = element_text(hjust = 1, size = 22), panel.spacing = unit(2, "lines"), strip.text.y = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), 
        legend.position="bottom", plot.margin = unit(c(0.8,0.5,6.2,0.5), "lines"), plot.title = element_text(size=28, hjust = 0.5))
ggsave(file = "E:/hazen_metagenome/results/physicochemistry_new.svg", plot=image, units="cm", width=12, height=16, scale = 2)

#construct phyloseq objects
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

func <- read.csv("E:/hazen_metagenome/matam_func_table.txt",sep="\t",skip=1,row.names=1)
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

physchem_vectors <- data.frame(model = names(ef_matam$vectors$r))
physchem_vectors <- cbind(physchem_vectors, scores(ef_matam, display="vectors"), pvalues=ef_matam$vectors$pvals)
sig_physchem_vectors <- physchem_vectors[physchem_vectors$pvalues < 0.05,]
nonsig_physchem_vectors <- physchem_vectors[physchem_vectors$pvalues > 0.05,]
plot(ordu_matam)
multiplier_matam <- ordiArrowMul(ef_matam)
image <- ggplot(data = NMDS_data_matam, aes(NMDS1, NMDS2)) + geom_point(aes(color=factor(Sample), shape=factor(Site)), size = 5) + coord_fixed() + labs(color = "Sample") +
  geom_segment(data=sig_physchem_vectors, aes(x=0, xend=NMDS1*multiplier_matam, y=0, yend=NMDS2*multiplier_matam), arrow = arrow(length = unit(0.5, "cm")), colour="black", size=1) +
  geom_segment(data=nonsig_physchem_vectors, aes(x=0, xend=NMDS1*multiplier_matam, y=0, yend=NMDS2*multiplier_matam), arrow = arrow(length = unit(0.5, "cm")), colour="black", linetype="twodash", size=1) +
  geom_text(data=physchem_vectors, aes(x=NMDS1*multiplier_matam+multiplier_matam*0.1, y=NMDS2*multiplier_matam+multiplier_matam*0.1, label=model),size=5)
ggsave("E:/Hazen_metagenome/results/matam_ordi.svg", plot=image, units="mm", width=300, height=300)


matam_func_data <- data.frame(otu_table(func))
matam_func_distance <- vegdist(t(matam_func_data), distance = "bray")
set.seed(42)
ordu_matam_func <- ordinate(func, "NMDS", distance = matam_func_distance)
ef_matam_func <- envfit(ordu_matam_func,sample_data(func),permu=10000)
NMDS_data_matam_func <- data.frame(sample_data(func))
NMDS_data_matam_func$NMDS1 <- ordu_matam_func$points[ ,1]
NMDS_data_matam_func$NMDS2 <- ordu_matam_func$points[ ,2]


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

new_RP_MAG_taxonomy <- read.table("E:/hazen_metagenome/processed_files/RP_MAGs_taxonomy.txt", sep="\t", header=T, na.strings=c("","NA"))
rownames(new_RP_MAG_taxonomy) <- new_RP_MAG_taxonomy$label
new_RP_MAG_tree <- as_tibble(treeio::drop.tip(new_ggtree2, tip=new_ggtree$label[which(new_ggtree$source %in% "repository")]))
new_RP_MAG_tree <- full_join(new_RP_MAG_tree %>% select(-(5:ncol(new_RP_MAG_tree))), new_RP_MAG_taxonomy, by="label")

#scale to binned reads in each sample
rel_abundance_pathways <- rel_abundances[[1]]
rel_abundance_pathways$bin <- rownames(rel_abundances[[1]])
rel_abundance_pathways <- rel_abundance_pathways[which(rel_abundance_pathways$bin %in% new_RP_MAG_tree$label),]
for (i in 1:6) {
  rel_abundance_pathways[,i] <-  sapply(rel_abundance_pathways[,i], function(x) x <- x/sum(rel_abundance_pathways[,i])*100)
}

sequencing <- read.csv("E:/hazen_metagenome/sample_sequencing_stats.csv",sep="\t")
sequencing <- sequencing[-which(sequencing$sample %in% "total"),]
sequencing <- sequencing[,-c(3,4)]
colnames(sequencing) <- c("sample", "Number of reads", "Total base pairs")
sequencing$sample <- sample_order
# sequencing$binned <- colSums(rel_abundances[[1]])*0.01*sequencing$n.reads
# sequencing$quality.bins <- colSums(rel_abundance_pathways[,-7])*0.01*sequencing$n.reads
sequencing <- reshape2::melt(sequencing)

image <- ggplot(sequencing, aes(x=factor(sample), y=value)) + geom_bar(stat="identity", color="black") + scale_x_discrete(limits=(sample_order)) + facet_wrap(~variable, scales="free") + 
  theme(axis.title.x = element_blank(), axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position="bottom", plot.title = element_text(size=28, hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Abundance (%)") + ggtitle("Sequencing statistics") + scale_y_continuous(labels = scales::comma)
ggsave(file = "E:/hazen_metagenome/results/sequencing_statistics.svg", plot=image, units="mm", width=300, height=200)

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


#set new levels for NCBI genomes
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
colnames(metacyc_processes) <- c("process", "Reference genomes", "MAGs", "Long-branching reference genomes")
metacyc_processes <- merge(metacyc_processes, longbranch_genome_list_3$metacyc$processes, by="process", all=T)
colnames(metacyc_processes)[5] <- "Long-branching MAGs"
metacyc_processes <- metacyc_processes[match(as.character(ncbi_genome_list_3$metacyc$processes$process), metacyc_processes$process),]

marker_processes <- merge(ncbi_genome_list_3$marker_genes$processes, genome_list_3$marker_genes$processes, by="process", all=T)
marker_processes <- merge(marker_processes, longbranch_ncbi_genome_list_3$marker_genes$processes, by="process", all=T)
colnames(marker_processes) <- c("process", "Reference genomes", "MAGs", "Long-branching reference genomes")
marker_processes <- merge(marker_processes, longbranch_genome_list_3$marker_genes$processes, by="process", all=T)
colnames(marker_processes)[5] <- "Long-branching MAGs"
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
  "Sulfite reduction to sulfide reversible", "Sulfite reduction to sulfide",
  "Sulfur compound oxidation","Sulfur oxidation","Sulfite oxidation to sulfate",
  "Phosphorus uptake", "Polyphosphate processing","Mercury resistance","Mercury methylation","Arsenic resistance",
  "Pb Zn resistance or homeostasis","Cd Pb resistance","Cd Co Cu Pb Zn resistance or homeostasis",
  "Cu resistance","Multimetal resistance regulation"))

all_processes$process <- factor(all_processes$process, levels=process_levels)
all_processes <- all_processes[order(all_processes$process),]
all_processes$category <- c(rep("Metal homeostasis", 8), rep("Biogeochemical cycles", 13), rep("Energy metabolism", 7), rep("Degradation", 9), rep("Biosynthesis", 9))
all_processes$category <- factor(all_processes$category, levels=rev(unique(all_processes$category)))

#separate chi-squared tests for MAGs and LBMs
MAG_processes_chi <- all_processes[c(2,3)]
rownames(MAG_processes_chi) <- all_processes$process
MAG_processes_chi[is.na(MAG_processes_chi)] <- 0
MAG_processes_chi <- MAG_processes_chi[rowSums(MAG_processes_chi) > 0,]
cisq_test_MAG_processes <- chisq.test(MAG_processes_chi, simulate.p.value = T, B = 10000)
contrib_MAG_processes <- 100*cisq_test_MAG_processes$residuals^2/cisq_test_MAG_processes$statistic
#subset to processes that contribute more than equally (100 % / 47 processes) to the total chi-squared score
most_important_MAG_processes <- cisq_test_MAG_processes$residuals[rowSums(round(contrib_MAG_processes,3)) > (100/47),]
levels(rownames(most_important_MAG_processes)) <- process_levels
most_important_MAG_processes <- melt(most_important_MAG_processes)

LBM_processes_chi <- all_processes[c(2,5)]
rownames(LBM_processes_chi) <- all_processes$process
LBM_processes_chi[is.na(LBM_processes_chi)] <- 0
LBM_processes_chi <- LBM_processes_chi[rowSums(LBM_processes_chi) > 0,]
cisq_test_LBM_processes <- chisq.test(LBM_processes_chi, simulate.p.value = T, B = 1000)
contrib_LBM_processes <- 100*cisq_test_LBM_processes$residuals^2/cisq_test_LBM_processes$statistic
#subset to processes that contribute more than equally (100 % / 47 processes) to the total chi-squared score
most_important_LBM_processes <- cisq_test_LBM_processes$residuals[rowSums(round(contrib_LBM_processes,3)) > (100/47),]
levels(rownames(most_important_LBM_processes)) <- process_levels
most_important_LBM_processes <- melt(most_important_LBM_processes)

image <- ggplot(data = most_important_MAG_processes, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color="white") + coord_flip() + 
  scale_fill_gradient2("Pearson residual", low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.y = element_text(size = 16), strip.text.y = element_text(angle=0, hjust=0), axis.text.x = element_text(angle=45, hjust=1), panel.border = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=28, hjust = 0.5))
ggsave(file = "E:/hazen_metagenome/results/chi_squared_heatmap_MAGs.svg", plot=image, units="cm", width=10, height=12, scale=2.5)

image <- ggplot(data = most_important_LBM_processes, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color="white") + coord_flip() + 
  scale_fill_gradient2("Pearson residual", low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.y = element_text(size = 16), strip.text.y = element_text(angle=0, hjust=0), axis.text.x = element_text(angle=45, hjust=1), panel.border = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=28, hjust = 0.5))
ggsave(file = "E:/hazen_metagenome/results/chi_squared_heatmap_LBMs.svg", plot=image, units="cm", width=10, height=12, scale=2.5)

compare_processes_chi <- all_processes[c(3,5)]
rownames(compare_processes_chi) <- all_processes$process
compare_processes_chi[is.na(compare_processes_chi)] <- 0
compare_processes_chi <- compare_processes_chi[rowSums(compare_processes_chi) > 0,]
cisq_test_compare_processes <- chisq.test(compare_processes_chi, simulate.p.value = T, B = 1000)

all_processes[,2:5] <- sweep(all_processes[,2:5], 2, c(2486,55,53,18), FUN = '/')
all_processes[,2:5] <- round(all_processes[,2:5]*100, 0)
all_processes[is.na(all_processes)] <- 0
all_processes_melt <- melt(all_processes)
all_processes_melt$variable <- factor(all_processes_melt$variable, levels=rev(colnames(all_processes)[2:5]))
#remove LBM groups from the figure because of changes to the MS
all_processes_melt <- all_processes_melt[which(all_processes_melt$variable %in% c("MAGs", "Reference genomes")),]

image <- ggplot(all_processes_melt, aes(y=value, x=process, fill=variable)) + geom_col(position = position_dodge2(width = 0.75, preserve = "single"))  +#geom_point(aes(size=value), color="steelblue", shape=1, stroke=3) +
  coord_flip() + scale_fill_manual(name="Collection", values=c("#b2df8a", "#33a02c")) + 
  facet_grid(category~., scales="free", space="free") + ylab("Present in % of genomes") + guides(fill = guide_legend(nrow = 2,byrow = F, reverse = T)) +
  theme(axis.text.y = element_text(size = 16), strip.text.y = element_text(angle=0, hjust=0), panel.border = element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size=28, hjust = 0.5),
        legend.position = "bottom")
ggsave(file = "E:/hazen_metagenome/results/all_processes_bars.svg", plot=image, units="cm", width=10, height=12, scale=3)

#summarize stats of genomes
MAG_summary_table <- genome_list_2[1]$stats
MAG_summary_table$est.Genome.Size..Mbp. <- MAG_summary_table$Genome.size..Mbp. / (MAG_summary_table$Completeness/100)
MAG_summary_table$LBM <- ifelse(rownames(MAG_summary_table) %in% rownames(longbranch_genome_list_2[1]$stats), "yes", "no")
MAG_summary_table <- rbind(MAG_summary_table, c(colMeans(MAG_summary_table[-ncol(MAG_summary_table)]), "NA"))
rownames(MAG_summary_table) <- c(as.character(genome_list_2[2]$taxonomy$Species), "Means")
write.csv(MAG_summary_table, "E:/hazen_metagenome/results/MAG_summary_table.csv", quote=F, row.names = T)
MAG_summary_table2 <- data.frame(Species = rownames(MAG_summary_table)[-nrow(MAG_summary_table)], LBM = MAG_summary_table[-nrow(MAG_summary_table),ncol(MAG_summary_table)])

#which MAGs take part in which processes
MAG_metacycs <- lapply(genomes2, '[[', 4)
MAG_markers <- lapply(genomes2, '[[', 7)
processes_of_interest <- c("sulfate reduction V", "sulfate activation for sulfonation", "sulfoacetaldehyde degradation I", "ammonia assimilation", "N_fixation", "Denitrification", "Nitrification_bacteria", 
                           "Sulfite_reduction_to_sulfide", "Sulfide_oxidation", "Mercury_resistance",  "DNRA_Polysulfide_reduction")
process_tax_list <- list(NULL)
for (i in 1:4) {
  process_tax_list[[i]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl(processes_of_interest[i],x$metacyc_pathway)), MAG_metacycs))], '[[', 2))
}
for (i in 5:length(processes_of_interest)) {
  process_tax_list[[i]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl(processes_of_interest[i],x$process)), MAG_markers))], '[[', 2))
}
process_tax_list[[length(process_tax_list)+1]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl("nrfH|nir[B-D]",x$marker_gene)), MAG_markers))], '[[', 2))
process_tax_list[[length(process_tax_list)+1]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl("nosZ",x$marker_gene)), MAG_markers))], '[[', 2))
process_tax_list[[length(process_tax_list)+1]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl("nirK",x$marker_gene)), MAG_markers))], '[[', 2))
process_tax_list[[length(process_tax_list)+1]] <- do.call(rbind, lapply(genomes2[which(mapply(function(x) any(grepl("nar[G-J]",x$marker_gene)), MAG_markers))], '[[', 2))
names(process_tax_list) <- c(processes_of_interest, "DNRA", "Nitrous oxide reduction", "Nitrite reduction to nitric oxide", "Nitrate reduction to nitrite")
process_tax_list <- bind_rows(process_tax_list, .id="process")
process_tax_list <- merge(process_tax_list, data.frame(otu_table(RP_MAG_phylo)), by.x="label", by.y="row.names", all.x=T)
process_tax_list$mean <- rowMeans(as.matrix(process_tax_list[10:15]))
process_tax_list <- process_tax_list[-c(10:15)]
process_tax_list <- process_tax_list[order(process_tax_list$process, process_tax_list$mean, decreasing =T),]
process_tax_list <- join(process_tax_list, MAG_summary_table2, by = "Species")[-1]
write.csv(process_tax_list, "E:/hazen_metagenome/results/processes_of_interest_taxa.csv", quote=F, row.names = F)
process_phylum_list <- process_tax_list
process_phylum_list$Phylum <- as.character(process_phylum_list$Phylum)
process_phylum_list$Phylum[which(process_phylum_list$Phylum %in% "Proteobacteria")] <- as.character(process_phylum_list$Class[which(process_phylum_list$Phylum %in% "Proteobacteria")])
process_phylum_list <- aggregate(mean~process+Phylum, process_phylum_list, sum)
process_phylum_list <- process_phylum_list[order(process_phylum_list$process, process_phylum_list$mean, decreasing = T),]
write.csv(process_phylum_list, "E:/hazen_metagenome/results/processes_of_interest_phyla.csv", quote=F, row.names = F)
process_LBM_list <- process_tax_list
process_LBM_list <- aggregate(mean~process+LBM, process_LBM_list, sum)
process_LBM_list <- process_LBM_list[order(process_LBM_list$process, process_LBM_list$mean, decreasing = T),]
write.csv(process_LBM_list, "E:/hazen_metagenome/results/processes_of_interest_LBM.csv", quote=F, row.names = F)

#review checks
#1. add N50 to Table S4
new_S4 <- do.call(rbind,lapply(genomes2, '[[', 1))
new_S4 <- merge(new_S4, genome_gene_lengths[[1]][,c(1,4,5)], by.x="row.names", by.y="bins")
new_S4 <- new_S4[order(new_S4$Completeness, decreasing = T),]
write.csv(new_S4, "E:/hazen_metagenome/results/new_Table_S4.csv", quote=F, row.names = F) #added as a column to the Table S4 in SI

#2. Check how many proteins are missing in the LBM compared to others and is there a significant difference

raw_files2 <- unlist(raw_files)
raw_files2 <- raw_files2[which(grepl("Bin",raw_files2))]
LBM_n_RP <- list(NULL)
for (i in 1:length(names(longbranch_genomes))) {
  binname <- names(longbranch_genomes)[i]
  LBM_n_RP[[i]] <- data.frame(Bin=binname, n_RP=length(grep(paste0(binname,"$"), raw_files2)), LBM="yes")
}

non_LBM <- names(genomes2[!(names(genomes2) %in% names(longbranch_genomes))])
non_LBM_n_RP <- list(NULL)
for (i in 1:length(non_LBM)) {
  binname <- non_LBM[i]
  non_LBM_n_RP[[i]] <- data.frame(Bin=binname, n_RP=length(grep(paste0(binname,"$"), raw_files2)), LBM="no")
}

RP_comparison <- rbind(do.call(rbind, LBM_n_RP), do.call(rbind, non_LBM_n_RP))
ggplot(RP_comparison, aes(x=LBM, y=n_RP)) + geom_boxplot()
summary(aov(n_RP ~ LBM, data=RP_comparison))

#3. Check if there is a difference between the lengths of the alignments (how many gaps)
sequences <- data.frame(Bin=gsub(">", "", outfile[which(grepl("Bin", outfile))]), seq=nchar(gsub("-","",outfile[which(grepl("Bin", outfile))+1])))
sequences <- merge(sequences, RP_comparison, by="Bin")

ggplot(sequences, aes(x=LBM, y=seq)) + geom_boxplot()
summary(aov(seq ~ LBM, data=sequences))

#4. Check branch length against n missing bases
sequences$ngap <- nchar(gsub("[ATGC]","",outfile[which(grepl("Bin", outfile))+1]))
sequences <- merge(sequences, branchlength_comparison[which(branchlength_comparison$source %in% "metagenome"),c("label", "branch.length")], by.x="Bin", by.y="label")
ggplot(sequences, aes(x=branch.length, y=ngap/6099*100)) + geom_point(aes(color=factor(LBM)), size=3) + geom_vline(xintercept=repository_quantile) + geom_hline(yintercept = 25) +
  xlab("Terminal branch length") + ylab("% gaps in alignment") + geom_smooth(method="lm") + scale_color_manual(values=c("red", "black"), name="LBM?")

robust_model <- lmRob(log10(ngap) ~ branch.length, data=sequences)
sequences2 <- sequences$LBM[which(sequences$ngap/6099*100 < 25)]
table(sequences2)

#5. Check the fragmentation of the RP in the bins

#on how many elements are the RP in each bin?

#following bash code run on CAC
#for genome in $(tail -n +2 processed_files/RP_MAGs_taxonomy.txt | awk {'print $1'}):
#  do
#for RP in $(awk {'print $1'} processed_files/ribosomal_proteins/RP_annotations.txt):
#  do
#grep ${RP} SAMPLES-SUMMARY-CONCOCT2/bin_by_bin/${genome}/${genome}-gene_calls.txt | awk -v x=${genome} -v y=${RP} -v OFS="\t" '{print x, y, $2}' >> processed_files/ribosomal_proteins/RP_contigs_in_bins.txt
#done
#done

RP_frag <- read.csv("E:/hazen_metagenome/processed_files/ribosomal_proteins/RP_contigs_in_bins.txt", header=F, sep="\t")
RP_unique_contigs <- data.frame(setDT(RP_frag)[, .(count = uniqueN(V3)), by = V1])
RP_unique_contigs <- RP_unique_contigs[order(as.numeric(RP_unique_contigs$count)),]
RP_unique_contigs$V1 <- factor(RP_unique_contigs$V1, levels=RP_unique_contigs[order(as.numeric(RP_unique_contigs$count)),1])
RP_unique_RP <- data.frame(setDT(RP_frag)[, .(count = uniqueN(V2)), by = V1])
RP_frag2 <- merge(RP_unique_contigs, RP_unique_RP, by="V1")
RP_frag2$LBM <- ifelse(RP_frag2$V1 %in% do.call(rbind, LBM_n_RP)$Bin,"yes","no")

#how robust are the trees for each protein family:
#is there a significant connection between the top 10 closest relatives and number of unique RP contigs

fragmented_bins <- as.character(RP_frag2$V1)

new_trees <- list(NULL)
i = 1
for (uniqueRP in RP_annotations$Pfam) {
  treefile <- paste0("E:/hazen_metagenome/SAMPLES-SUMMARY-CONCOCT2/ribosomal_proteins/", uniqueRP, "_trimal.tre")
  new_trees[[i]] <- ape::read.tree(treefile)
  i = i+1
}

outlists <- list(NULL)
for (i in 1:length(fragmented_bins)) {
  outlist <- list(NULL)
  fragmented_bin <- fragmented_bins[i]
  for (j in 1:16) {
    binid <- paste0(RP_annotations$Pfam[j],"_",fragmented_bin)
    if (any(new_trees[[j]]$tip.label %in% binid)) {
      distancestobin <- cophenetic.phylo(new_trees[[j]])[binid,]
      distancestobin <- distancestobin[order(as.numeric(distancestobin))]
      names(distancestobin) <- gsub("^(.*?_).*?","",names(distancestobin))
      binitself <- names(distancestobin[grep(fragmented_bin, names(distancestobin))])
      matchnames <- names(head(distancestobin[-grep("Bin", names(distancestobin))],10))
      outlist[[j]] <- data.frame(append(binitself, matchnames))
    }
  }
  outlist <- tail(sort(table(do.call(rbind, outlist))))
  outlists[[i]] <- outlist
}

names(outlists) <- fragmented_bins

outlists2 <- list(NULL)
for (i in 1:length(outlists)) {
  percent_top10 <- outlists[[i]]
  total_RP <- percent_top10[which(names(percent_top10) %in% names(outlists)[i])]
  percent_top10 <- percent_top10[-grep("Bin", names(percent_top10))]
  percent_top10 <- round(max(percent_top10)/as.numeric(total_RP)*100,2)
  if (length(percent_top10) == 0) {percent_top10 <- 0}
  outlists2[[i]] <- percent_top10
}

RP_frag2 <- cbind(RP_frag2,do.call(rbind,outlists2))
colnames(RP_frag2)[5] <- "perc.same.top10"

ggplot(RP_frag2, aes(x=V1, y=count.x, color=factor(count.y))) + geom_point(aes(shape=factor(LBM)), size=5) + xlab("MAG") + ylab("Unique RP contigs (n)") +
  scale_color_brewer(palette="Set2", name="Ribosomal\nproteins (n)") + scale_shape_manual(values=c(16,17), name="LBM?") +  theme(axis.text.x = element_text(angle = 90, hjust=0))
RP_fragmentation_model <- lm(log10(count.x) ~ as.numeric(perc.same.top10), data=RP_frag2)
summary(RP_fragmentation_model)

#6. T-tests of differential abundance of individual marker genes
sample_data(marker_phylo)$O2yesno <- ifelse(sample_data(marker_phylo)$O2.mgL>0,"1","0")
mt(marker_phylo, "Site", method="fdr", test="t")
mt(marker_phylo, "O2yesno", method="fdr", test="t")

#7. check if phyla proportions are similar in matam data
# get % proportions as fractions
matam_long <- reshape2::melt(as.matrix(otu_table(phyla_matam)/100))
colnames(matam_long) <- c("phylum_otu", "sample", "value")
