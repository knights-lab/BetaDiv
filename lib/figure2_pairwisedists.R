##### FIG 2 DISTANCE METRICS : Pairwise Plot #####
# Depends on being run within driver script 'beta_review_figs.R'

# Compute distances:
## Bray-Curtis
bray <- vegdist(t(rare_otu), method="bray")
## Kulczynski Distance
kulcz <- vegdist(t(rare_otu), method="kulczynski")
## Manhattan Distance
manh <- vegdist(t(rare_otu), method="manhattan")
## Cosine dissimilarity
cos <- cosine(as.matrix(rare_otu)) # measures cosine similarity
cos <- as.dist(1-cos)              # get the dissimilarity and convert to dist object
## Jaccard
jacc <- vegdist(t(rare_otu), method="jaccard", binary = T)
## Unweighted UniFrac (load from IMP analyses GitHub)
u_unifrac <- read.delim('data/imp/unweighted_unifrac_dm.txt', sep='\t', row=1)
u_unifrac <- as.dist(u_unifrac[colnames(rare_otu),colnames(rare_otu)])
## Weighted UniFrac (load from IMP analyses GitHub)
w_unifrac <- read.delim('data/imp/weighted_unifrac_dm.txt', sep='\t', row=1)
w_unifrac <- as.dist(w_unifrac[colnames(rare_otu),colnames(rare_otu)])
## Euclidean Distance
euclid <- vegdist(t(rare_otu), method="euclidean")
## Aitchison Distance (Euclidean distance in CLR transformed data, requires non-zero data)
aitch <- vegdist(t(rare_otu), method="aitchison", pseudocount=1)
## Robust Aitchison
deicode <- read.table("results/deicode_out_imp/distance-matrix.tsv", sep="\t", header=T, row=1)
robust_a <- as.dist(deicode)
## Chi-square Distance
chisq <- vegdist(t(rare_otu), method="chisq")


# Pairwise Scatterplots Visual
## color by within group or outside group distance
group_dist <- dist(as.numeric(as.factor(meta$Sample.Group)))
group_dist <- as.factor(c(group_dist) < 1) # T/F list of w/in group or not
levels(group_dist) <- c("out", "in")
## create data frame of distances for pairwise graphing
dist_df <- data.frame(BrayCurtis=c(bray), Kulczynski=c(kulcz), Manhattan=c(manh), Cosine=c(cos), Jaccard=c(jacc), 
                      uUniFrac=c(u_unifrac), wUniFrac=c(w_unifrac), Euclidean=c(euclid),
                      Aitchison=c(aitch), rAitchison=c(robust_a), ChiSq=c(chisq),
                      Color_Within_Group=c(group_dist))
## for trying the same plot but with rank (this actually results in worse clouds)
# rank_df <- as.data.frame(apply(dist_df[,-ncol(dist_df)], 2, rank))
# rank_df$Color_Within_Group <- c(group_dist)

## CUSTOM FUNCTION for upper triangle in ggpairs
custom_corr <- function (data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  my_corr <- cor(x, y, method="pearson")
  adjx <- mean(x)
  adjy <- mean(y)
  
  p <- ggplot(data) +
    geom_blank(mapping) +
    geom_smooth(mapping, method="lm", linetype="dashed") +
    annotation_custom(grid::textGrob(paste0("corr: ", round(my_corr, 2)), gp=gpar(fontsize=16)),
                      #xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  p
}

## ggplot for publication
pairwise_plt <- ggpairs(dist_df, columns=1:11, ggplot2::aes(col=Color_Within_Group, fill=Color_Within_Group, alpha=0.4),
                        lower = list(continuous = wrap("points", alpha=0.1, pch=16)),
                        upper = list(continuous = custom_corr) ) + 
  scale_color_manual(values=c("#9d02d7", "#ffb14e")) +
  scale_fill_manual(values=c("#9d02d7", "#ffb14e")) +
  theme_bw() +
  theme(title = element_text(face="bold"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(family="Helvetica", size=16))
png("Figures_Tables/full_res_graphs/dist_pairwise.png", width=12, height=12, res=700, units="in")
pairwise_plt
dev.off()

print('Done with figure 2.')


