##Assignment 4
#### Packages -----
library(tidyverse)
library(Biostrings)
library(phylotools)
library(ape)
library(ggtree)
library(picante)
library(phytools)
library(vegan)
library(viridis)

#### Retrieve public data via BOLD API ------
# retrieved on November 26, 2025
# Carabidae dataset using tropical countries where Carabidae data is available according to BOLD Taxonomy Browser
#dfCarabidae <- read_tsv("https://www.boldsystems.org/index.php/API_Public/combined?taxon=Carabidae&geo=Costa%20Rica|Ghana|Philippines|South%20Africa|Peru&format=tsv")
#write_tsv(dfCarabidae, "../data/df_Carabidae.tsv")

df_Carabidae <- read_tsv("../data/df_Carabidae.tsv")

#### Prepare dataset for phylogenetic analyses -----
Carabidae <- df_Carabidae %>%
  filter(!is.na(bin_uri) & !is.na(elev) & !is.na(lat) & markercode == "COI-5P") %>%
  filter(lat >= -23.5 & lat <= 23.5) %>% #tropical latitudes to account for countries that span multiple zones
  select("processid", "sampleid", "bin_uri", "elev", "nucleotides","subfamily_name", "genus_name", "species_name")

# remove sequences below 500bp
Carabidae <- Carabidae %>%
  mutate(seq_length = str_count(nucleotides, "[ACGTacgt]")) %>%
  filter(seq_length >= 500)

# modify BIN to convert ":" to "_" for Newick file compatibility
Carabidae <- Carabidae %>%
  mutate(bin = gsub(":", "_", bin_uri))

#### Generate elevation bins -----
# set minimum elevation as ground-level (0)
min_elev <- 0
# determine maximum from dataset
max_elev <- max(Carabidae$elev)

# construct series of intervals increasing by 100m
elev_bins <- seq(min_elev, max_elev, by = 100)

# associate elev_bins into Carabidae dataset
Carabidae <- Carabidae %>%
  mutate(
    elev_bin = cut(
      elev,
      breaks = elev_bins,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 5 # number of allowable digits
    )
  ) 

# compute midpoints for each bin
midpoints <- (head(elev_bins, -1) + tail(elev_bins, -1)) / 2

# attach midpoints
Carabidae <- Carabidae %>%
  mutate(elev_midpoint = midpoints[as.numeric(elev_bin)])

rm(df_Carabidae, elev_bins, midpoints)

#### Select BIN representative -----
# need to select a single representative high-quality sequence for phylogenetic tree building and analyses
# use Barcode Index Number as a proxy for species

# count unique BINs in dataset
Carabidae %>%
  summarise(distinct_bins = n_distinct(bin_uri))
#398 distinct BINs

# summary of sequence lengths
summary(nchar(Carabidae$nucleotides))

# convert 'nucleotides' string to sequence compatible with Biostrings
Carabidae <- as.data.frame(Carabidae)
seq <- DNAStringSet(Carabidae$nucleotides)

# generate sequence data including sequence length and ambiguous base count
seq_info <- data.frame(
  processid = Carabidae$processid,
  bin = Carabidae$bin,
  length = width(seq),
  N = letterFrequency(seq, "N"),
  seq = seq, 
  stringsAsFactors = FALSE
)

# select highest quality sequence as longest read with least ambiguous bases "N"
high_quality_seq <- seq_info %>%
  group_by(bin) %>%
  arrange(desc(length), N) %>%
  dplyr::slice(1) %>%
  ungroup()

# generate BIN representative sequences dataset
bin_rep <- high_quality_seq %>%
  mutate(header = paste(bin, processid, sep = "|")) %>%
  select(seq.name = header, seq.text = seq)

dat2fasta(dat = bin_rep, outfile = "../data/bin_rep.fasta")

# separate pipe delimited header to inner join bin_rep with Carabidae using processid (unique identifer)
bin_rep <- bin_rep %>%
  separate_wider_delim(
    cols = seq.name,
    delim = "|",
    names = c("bin", "processid")) %>% 
  dplyr::rename(nucleotides = seq.text)

# BIN representative data frame
df_binrep <- inner_join(Carabidae, bin_rep, by = "processid") %>%
  select(-c(bin.y, nucleotides.y)) %>%
  dplyr::rename(bin = bin.x, nucleotides = nucleotides.x)

rm(seq, seq_info, high_quality_seq, bin_rep)

#### Construct Community Matrix -----
# extract elevation bin and BIN from full dataset and count BIN frequency per elevation bin
community_matrix <- Carabidae %>%
  group_by(elev_bin, bin) %>%
  summarise(freq = n()) %>%
  pivot_wider(id_cols = elev_bin, names_from = bin, values_from = freq, values_fill = 0) %>%
  column_to_rownames(var="elev_bin") %>%
  as.data.frame()

# convert community_matrix to logical presence/absence matrix
community_matrix[community_matrix > 0] <- 1
#### Visualize Maximum Likelihood BIN Phylogeny -----
# 1. Align and trim the BIN representative FASTA file
# 2. Run model selection in MEGA12
#    • Best-fit model: GTR + Gamma + Invariant sites (GTR+G+I)
# 3. Reconstruct phylogenetic tree using IQ-TREE2
#    • 1000 bootstrap replicates

# BIN removed after initial tree inspection:
# - BIN BOLD:AEE8199 was identified as contaminated via BOLD v5 ID Engine
# - Removed from aligned FASTA
# - Tree was reconstructed using the same parameters

Carabidae <- Carabidae %>%
  filter(bin_uri != "BOLD:AEE8199")
df_binrep <- df_binrep %>%
  filter(bin_uri != "BOLD:AEE8199")
community_matrix <- community_matrix %>%
  select(-BOLD_AEE8199)

# save a copy of dataframes
write_tsv(Carabidae, "../data/Carabidae.tsv")
write_tsv(df_binrep, "../data/df_binrep.tsv")

# read in tree with original labels BIN|processid
# keep identifier (processid) for record-tracking / misidentification or contamination checks
tree_original <- read.tree("../data/bin_rep_aln.fasta.treefile")

# create a copy of the tree with BIN label only, 
tree <- tree_original
tree$tip.label <- sapply(strsplit(tree$tip.label, "\\|"), `[`, 1)

# plot phylogenetic tree
ggtree(tree) + 
  geom_tiplab(size=0)

ggsave("../figs/MLtree.png", plot = last_plot() , width = 10, height = 6, dpi = 100)

rm(tree_original)

#### Compare Community Matrix and Phylogenetic Tree -----
# compare community matrix & phylogenetic tree
clean_tree <- match.phylo.comm(phy = tree, comm = community_matrix)$phy
clean_comm <- match.phylo.comm(phy = tree, comm = community_matrix)$comm

#### Calculate Phylogenetic alpha (𝛼) Diversity (PD) across elevation bins -----
# check if tree is rooted
is.rooted(clean_tree)
# FALSE

# root tree using midpoint
clean_tree_rooted <- midpoint_root(clean_tree)
is.rooted(clean_tree_rooted)
# TRUE

rm(clean_tree)

# calculate Faith's phylogenetic diversity (PD)
PD <- pd(samp = clean_comm, tree = clean_tree_rooted, include.root = TRUE)

# calculate standard effect size (ses) of PD for each site (elev_bin)
ses.pd <- ses.pd(samp = clean_comm, 
                 tree = clean_tree_rooted, 
                 null.model = "taxa.labels", 
                 runs = 1000)

# include elevation midpoints to plot elev vs PD or SR
# extract elev_bins present in community matrix
comm_elev_bins <- rownames(community_matrix)

# split elev_bin string to upper and lower bounds
comm_elev_bin_split <- str_split_fixed(comm_elev_bins, ",", 2)

# remove bracket character and convert lower and upper bounds to numeric
lower <- as.numeric(str_remove_all(comm_elev_bin_split[,1], "\\["))
upper <- as.numeric(str_remove_all(comm_elev_bin_split[,2], "\\)|\\]"))

# calculate midpoints
comm_elev_midpoints <- (lower + upper) / 2

# attach midpoints
PD$elev_midpoint <- comm_elev_midpoints

rm(comm_elev_bin_split, lower, upper)

# model PD against elevation
lm_elev_PD <- lm(PD ~ elev_midpoint, data = PD)
summary(lm_elev_PD)

# plot elevation vs PD
ggplot(data = PD, aes(x = elev_midpoint, y = PD)) + 
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ x, se = TRUE, show.legend = FALSE) +
  labs(title = "Phylogenetic Diversity across Tropical Elevations",
       x = "Elevation (m)",
       y = "Phylogenetic Diversity") + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../figs/elev_vs_PD.png", plot = last_plot() , width = 10, height = 6, dpi = 100)

# model SR against elevation
lm_elev_SR <- lm(SR ~ elev_midpoint, data = PD)
summary(lm_elev_SR)

# plot elevation vs SR
ggplot(data = PD, aes(x = elev_midpoint, y = SR)) + 
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ x, se = TRUE, show.legend = FALSE) +
  labs(title = "BIN Richness across Tropical Elevations",
       x = "Elevation (m)",
       y = "BIN Richness") + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../figs/elev_vs_SR.png", plot = last_plot() , width = 10, height = 6, dpi = 100)

# correlation test between PD & SR
cor.test(PD$PD, PD$SR)
# p-value = 3.946e-16

# model PD against SR
lm_PD_SR <- lm(PD ~ SR, data = PD)

#plot PD vs SR
ggplot(data = PD, aes(x = SR, y = PD)) + 
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ x, se = TRUE, show.legend = FALSE) +
  labs(title = "Phylogenetic Diversity vs. BIN Richness",
       x = "Phylogenetic Diversity",
       y = "BIN Richness") + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../figs/SR_vs_PD.png", plot = last_plot() , width = 10, height = 6, dpi = 100)

#extract residuals
PD$residuals <- resid(lm_PD_SR)

#plot elevation vs residuals of PD & SR
ggplot(data = PD, aes(x = elev_midpoint, y = residuals)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "blue") + # straight line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "Residual Phylogenetic Diversity Across Elevation",
       x = "Elevation (m)",
       y = "Residuals of PD & SR") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../figs/elev_vs_resid.png", plot = last_plot() , width = 10, height = 6, dpi = 100)

rm(ses.pd, lm_PD_SR)

#### Calculate Phylogenetic beta (𝛽) Diversity (PBD) across elevation bins -----
# calculate dissimilarity between communities using Phylosor
sor <- phylosor(samp = clean_comm,
                tree = clean_tree_rooted)

# calculate phylogenetic distance between communities using Unifrac
uni <- unifrac(comm = clean_comm, tree = clean_tree_rooted)

# create list of distance matrics for Principal Coordinate Analysis (PCoA)
dist_list <- list(
  PhyloSor = sor,
  UniFrac  = uni
)

# loop over dist_list to run PCoA and plot
lapply(names(dist_list), function(metric) {
  
  # run PCoA
  res <- pcoa(dist_list[[metric]], correction = "cailliez")
  
  # retrieve first 2 principal axes
  coords <- as.data.frame(res$vectors[, 1:2])
  colnames(coords) <- c("Axis1", "Axis2")
  
  # add elevation factor
  coords$Elevation <- factor(rownames(coords), levels = rownames(coords))
  
  # calculate variance explained
  eig <- res$values$Eigenvalues
  var_exp <- round(eig / sum(eig) * 100, 2)
  
  # plot
  p <- ggplot(coords, aes(x = Axis1, y = Axis2, color = Elevation)) +
    geom_point(size = 2) +
    scale_color_viridis_d(option = "C") +
    labs(
      title = paste0("Phylogenetic Beta Diversity (", metric, ") Across Elevation"),
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)"),
      color = "Elevation Bin"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # save plot
  ggsave(
    filename = paste0("../figs/PCoA_", metric, ".png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 100
  )
  
  # view plots
  return(p)
})

rm(dist_list)

# construct dataframe to group elevation bins into groups (low=<700m, mid=700m-1500m, high=>1500m) 
group <- data.frame(
  elev_group = cut(
    comm_elev_midpoints,
    breaks = c(min_elev, 700, 1500, (max_elev + 1)),
    labels = c("Low", "Mid", "High"),
    include.lowest = TRUE)
)

# conduct Permutational Multivariate Analysis of Variance (PERMANOVA) on Unifrac data
permanova_uni <- adonis2(uni ~ elev_group, 
                         data = group, 
                         permutations = 999)

permanova_uni

rm(group, min_elev, max_elev)

#### Calculate Net Relatedness Index (NRI) -----
# create a phylogenetic distance matrix
cophenDist <- cophenetic.phylo(clean_tree_rooted)

# estimate standardized effect size of the mean phylogenetic distances (MPD) on the cophenetic matrix
ses.mpd <- ses.mpd(clean_comm, 
                   cophenDist, 
                   null.model = "taxa.labels", 
                   abundance.weighted = FALSE, 
                   runs = 100)

# calculate NRI
NRI <- as.matrix(-1 * ((ses.mpd[,2] - ses.mpd[,3]) / ses.mpd[,4]))
rownames(NRI) <- row.names(ses.mpd)
colnames(NRI) <- "NRI"

rm(ses.mpd)

#### Calculate Nearest Taxon Index (NTI) -----
# estimate standardized effect size of the mean nearest neighbour phylogenetic distances (MNTD) on the cophenetic matrix
ses.mntd <- ses.mntd(clean_comm, 
                     cophenDist, 
                     null.model = "taxa.labels", 
                     abundance.weighted = FALSE, 
                     runs = 100) 

# calculate NTI
NTI <- as.matrix(-1 * ((ses.mntd[,2] - ses.mntd[,3]) / ses.mntd[,4]))
rownames(NTI) <- row.names(ses.mntd)
colnames(NTI) <- "NTI"

# convert NTI to dataframe to plot
NTI <- as.data.frame(NTI)

# attach elevation midpoints to plot elev vs NTI
NTI$elev_midpoint <- comm_elev_midpoints

# attach p-values to colourize significant elevations
NTI$p.value <- ses.mntd$mntd.obs.p

# add column for significance
NTI$Significance <- ifelse(NTI$p.value < 0.05, "Significant", "Not Significant")
NTI$Significance[is.na(NTI$Significance)] <- "Not Significant"

# convert significance to factor
NTI$Significance <- factor(NTI$Significance,
                          levels = c("Not Significant", "Significant"))

# plot elevation vs NTI
ggplot(data = NTI, aes(x = elev_midpoint, y = NTI)) + 
  geom_point(aes(colour = Significance)) +
  geom_smooth(method = "glm", formula = y ~ x, se = TRUE, show.legend = FALSE) +
  scale_color_manual(values = c("Not Significant" = "black", "Significant" = "red")) +
  labs(title = "Nearest Taxon Index across Tropical Elevations",
       x = "Elevation (m)",
       y = "NTI") + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))
  
ggsave("../figs/elev_vs_NTI.png", plot = last_plot() , width = 10, height = 6, dpi = 100)

rm(ses.mntd, cophenDist, comm_elev_bins, comm_elev_midpoints)