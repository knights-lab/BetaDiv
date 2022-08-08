# Figures for "Microbiome beta diversity analysis and visualization"
# Suzie Hoops and Dan Knights
# May 2022

# FEMS Microbiology Reviews publishing requirements
#    - Color figures are encouraged and free of charge.
#    - All lines should be 1.5 point (0.5 mm) wide (size=0.5+), broken
#      lines styles allowed to differentiate multiple plots.
#    - Letters and numbers must be non-serif font (Helvetica allowed), 16pt font.
#    - Symbols in figure must be 3mm in diameter. (size=3+)
#    - Lines drawn should not go through hollow symbols.
#    - Numbers in axis labels must have preceding zero if less than unity.
#    - Larger composite figures can occupy two columns.
#    - Print size 
#    - Minimum resolution for final publication (submit low res for peer review):
#      400 dpi for color images, 600 dpi for line drawings
#    - Consider also making a "feature image"



##### Set Up #####
# Loading libraries/packages
library(here)
library(vegan)
library(lsa)
library(ggplot2)
library(GGally)
library(gridExtra)
library(grid)
library(Rtsne)
library(ROptSpace)
set.seed(25)
# Working directory should be where this script is located (assuming same directory as 'data' and 'results')
curr_dir <- here()
if (grepl('Google', curr_dir) != T) {
  curr_dir <- rstudioapi::getActiveDocumentContext()$path
}
curr_dir <- gsub('/Volumes', '~', curr_dir)
curr_dir <- gsub('GoogleDrive', 'Google Drive', curr_dir)
curr_dir <- gsub('MyDrive', 'My Drive', curr_dir)
curr_dir <- gsub('beta_review_figs.r', '', curr_dir) # remove this file name
curr_dir <- gsub('bin/', '', curr_dir) # move out of the bin folder to parent folder
print(paste0('Setting working directory to script location: ', curr_dir))
setwd(curr_dir)



##### Helper functions #####
# Custom colors and shapes to match the IMP paper
g_col <- c("#e38cbd", "#7e7e7e", "#d65ca0", "#4ea99b")
names(g_col) <- c("Hmong2nd", "Control", "HmongThai", "KarenThai")
g_shp <- c(17, 17, 16, 16)
names(g_shp) <- c("Hmong2nd", "Control", "HmongThai", "KarenThai")

# Custom colors and shapes for soil
soil_fill <- c("black","#0933d9","#2f56f9","#748efb","#adbdfc","white")
names(soil_fill) <- c("<4", "4-5", "5-6", "6-7", "7-8", ">8")
soil_shp <- c(22,24,25,21,23,23)
names(soil_shp) <- c("<4", "4-5", "5-6", "6-7", "7-8", ">8")
soil_col <- c("white","white","white","white","white","black")
names(soil_col) <- c("<4", "4-5", "5-6", "6-7", "7-8", ">8")

# PCoA plot for distance comparison
## expects a distance object
pcoa_plot <- function (dist_obj, k=2, meta_df, title_text="", legendpos="none", adjx=1, adjy=1) {
  tmp_pc <- cmdscale(dist_obj, k=k, eig=T)
  pc <- as.data.frame(tmp_pc$points)
  colnames(pc) <- paste0("PC", 1:k)
  pc$imp_group <- meta_df$Sample.Group
  pc$SampleID <- meta_df$SampleID
  var <- (eigenvals(tmp_pc)/sum(eigenvals(tmp_pc)))[1:k]
  var <- round(var*100, 1)
  pc[,1] <- pc[,1]*adjx
  pc[,2] <- pc[,2]*adjy
  p <- ggplot(pc, aes(x=PC1, y=PC2)) +
    geom_point(size=4, alpha=0.7, aes(shape=imp_group, col=imp_group)) +
    scale_color_manual(values = g_col[unique(pc$imp_group)]) +
    scale_shape_manual(values = g_shp[unique(pc$imp_group)]) +
    labs(title=title_text, x=paste0("PC 1 [", var[1], "%]"), y=paste0("PC 2 [", var[2], "%]"), col="Group", shape="Group") +
    theme_bw() + 
    theme(legend.position = legendpos,
          legend.box.background = element_rect(colour = "black"),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(0,0,0,0,'cm'),
          panel.grid = element_blank(),
          text = element_text(family="Helvetica", size=16))
  return(p)
}

# 2D plots for ordination comparison
## expects a dataframe with ordination dimensions, ordered the same as the metadata
## expects IMP metadata
plot_2d <- function (dims, meta_df, title_text="", eig=NULL, legendpos="none") {
  df_2d <- as.data.frame(dims)
  colnames(df_2d) <- paste0("Dim", 1:ncol(df_2d))
  df_2d$imp_group <- meta_df$Sample.Group
  df_2d$sampleid <- rownames(df_2d)
  p <- ggplot(df_2d, aes(x=Dim1, y=Dim2)) +
    geom_point(size=4, alpha=0.7, aes(shape=imp_group, col=imp_group)) +
    scale_color_manual(values = g_col[unique(df_2d$imp_group)]) +
    scale_shape_manual(values = g_shp[unique(df_2d$imp_group)]) +
    theme_bw() +
    theme(legend.position = legendpos,
          legend.box.background = element_rect(colour = "black"),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(0,0,0,0,'cm'),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(family="Helvetica", size=16))
  if (is.null(eig)) {
    p <- p + labs(title=title_text, x="Dimension 1", y="Dimension 2", col="", shape="")
  } else {
    eig <- eig/sum(eig)
    var <- round(eig*100, 1)
    p <- p + labs(title=title_text, x=paste0("Dimension 1 [",var[1],"%]"), y=paste0("Dimension 2 [",var[2],"%]"), col="", shape="")
  }
  return(p)
}

# Legend Only Plot
plot_imp_legend <- function(meta_df) {
  p_df <- as.data.frame(matrix(rep(0, nrow(meta_df) * 2), nrow=nrow(meta_df)))
  colnames(p_df) <- c("Dim1", "Dim2")
  p_df$imp_group <- meta_df$Sample.Group
  p <- ggplot(p_df, aes(x=Dim1, y=Dim2)) +
    geom_point(size=4, alpha=1, aes(shape=imp_group, col=imp_group)) +
    scale_color_manual(values = g_col[unique(p_df$imp_group)]) +
    scale_shape_manual(values = g_shp[unique(p_df$imp_group)]) +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          title = element_blank(),
          legend.position = c(0.5,0.5),
          legend.box.background = element_rect(colour = "black", size=1.2),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(0.25,0.25,0.25,0.25,'cm'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          text = element_text(family="Helvetica", size=16))
  p
}

# Procrustes plot
proc_plot <- function (prot_obj, meta_df, mytitle="Procrustes") {
  stacked <- data.frame(rbind(prot_obj$X[,1:2], prot_obj$Yrot[,1:2]))
  colnames(stacked) <- c("xdir", "ydir");
  stacked$imp_group <- rep(meta_df$Sample.Group, 2)
  segm <- data.frame(cbind(prot_obj$X[,1:2], prot_obj$Yrot[,1:2]))
  colnames(segm) <- c("xstart", "ystart", "xend", "yend");
  segm$group <- meta_df$Sample.Group
  p <- ggplot(stacked, aes(x=xdir, y=-1*ydir)) +
    geom_point(size=3, alpha=0.7, aes(col=imp_group, shape=imp_group)) +
    scale_color_manual(values = g_col[unique(stacked$imp_group)]) +
    scale_shape_manual(values = g_shp[unique(stacked$imp_group)]) +
    labs(title=mytitle, x="", y="", col="", shape="") +
    geom_segment(aes(x=xstart, y=-1*ystart, xend=xend, yend=-1*yend, col=group),
                 data=segm, arrow=arrow(length=unit(0.015, "npc"))) +
    annotate("label", x=-0.01, y=-0.07, size=6, label = paste0("protest signif: ", prot_obj$signif)) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          text = element_text(family="Helvetica", size=16))
  p
}



##### Load example data set: Vangay et al. U.S. Immigration westernizes human gut #####
# loading datafiles (obtained from GitHub)
meta <- read.table("data/imp/map.txt", sep="\t", header=T)
rownames(meta) <- meta$SampleID
rare_otu <- read.delim("data/imp/final_otu.txt", header=T, row=1) # OTU table rarefied to ~10mil counts
# refine to no first gen samples (these were all over the place)
meta <- meta[!(grepl("1st", meta$Sample.Group)),]    # this omits 137 Hmong1st and 236 Karen1st (278 total left)
meta <- meta[meta$SampleID %in% colnames(rare_otu),] # ignore 9 samples not matching otu table (269 total left)
rare_otu <- rare_otu[,meta$SampleID]
rare_otu <- rare_otu[rowSums(rare_otu) > 0,]         # restrict to non-empty OTU rows (6803 remian)



##### Load example data set: Lauber et al. 88 Soils #####
# loading datafiles (obtained from Qiita)
meta_soil <- read.table("data/soil/clean_map.txt", sep="\t", header=T)
meta_soil <- meta_soil[,c("SampleID","instrument_name", "library_construction_protocol",
                          "target_gene","collection_date","country","depth","elevation",
                          "latitude","longitude","ph","silt_clay","soil_moisture_deficit",
                          "soil_type","specific_location","texture")]
rownames(meta_soil) <- meta_soil$SampleID
soil_otu <- read.delim("data/soil/44766_clean_otus.txt", header=T, row=1) # OTU table normalized
identical(colnames(soil_otu), rownames(meta_soil))
# refine to samples without low seq counts
drop <- c("BB1", "GB5", "HI1", "LQ1", "MD3", "MD4") # drop low counts (6 samples)
meta_soil <- meta_soil[!(rownames(meta_soil) %in% drop),]
soil_otu <- soil_otu[,!(colnames(soil_otu) %in% drop)]
soil_otu <- soil_otu[rowSums(soil_otu) > 0,]
# create ph group variable
meta_soil$ph_group <- cut(meta_soil$ph, breaks=c(0,4,5,6,7,8,14), right=F)
levels(meta_soil$ph_group) <- c("<4", "4-5", "5-6", "6-7", "7-8", ">8")



##### Load example data set: Yatsunenko et al. Human gut across age & geography #####
# metadata
meta_gg <- read.table("data/global_gut/clean_map.txt", sep="\t", header=T)
meta_gg <- meta_gg[,c("Sample_ID","body_site","diet","sample_type","sex","weight","age", "geo_loc_name")]
rownames(meta_gg) <- meta_gg$Sample_ID
meta_gg <- meta_gg[!(is.na(meta_gg$age)),]
# unweighted unifrac distances
gg_uunifrac <- as.dist(read.delim("data/global_gut/clean_uw_unifrac_dist.txt", header=T, row=1))



##### Run Figure library scripts to generate figures #####
source("lib/figure1_litreview.R")
source("lib/figure2_pairwisedists.R")
source("lib/supplemental_outliers.R")
source("lib/figure3_dist_pcoa.R")
source("lib/figure4_visualization.R")
source("lib/figure5_constrained.R")
source("lib/figure6_artifacts.R")
print('All scripts complete.')

# ##### Graphical abstract ???
# mini_dist_df <- data.frame(BrayCurtis=c(bray), uUniFrac=c(u_unifrac),
#                            Aitchison=c(aitch), Color_Within_Group=c(group_dist))
# mini_pair_plt <- ggpairs(mini_dist_df, columns=1:3, ggplot2::aes(col=Color_Within_Group, fill=Color_Within_Group, alpha=0.4),
#                         lower = list(continuous = wrap("points", alpha=0.1, pch=16)),
#                         upper = list(continuous = custom_corr) ) + 
#   scale_color_manual(values=c("#9d02d7", "#ffb14e")) +
#   scale_fill_manual(values=c("#9d02d7", "#ffb14e")) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         text = element_text(family = "Helvetica", size=14))
# png("Figures_Tables/full_res_graphs/graphical_abstract.png", width=3, height=3, res=700, units="in")
# mini_pair_plt
# dev.off()

