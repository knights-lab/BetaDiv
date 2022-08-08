##### SUPP FIG 1 & 2 EXPERIMENT : How are outliers affecting the Euclidean distance cloud? ######
# Depends on being run within driver script 'beta_review_figs.R'


# Test: color by between the Chi-sq outliers or not:
tmp <- meta[,c("SampleID","Sample.Group")]
rownames(tmp) <- tmp$SampleID
tmp$Outlying <- 0
tmp[c("CS.253", "CS.378", "CS.283", "T.CS.149", "CS.339"),"Outlying"] <- 1:5
outl_dist <- dist(as.numeric(tmp$Outlying))
outl_dist <- as.factor(c(outl_dist) < 1)
levels(outl_dist) <- c("outliers", "typical")
## create data frame of distances for pairwise graphing
dist_df2 <- data.frame(BrayCurtis=c(bray), UnweightedUniFrac=c(u_unifrac), Euclidean=c(euclid),
                       ChiSquared=c(chisq), Color_Within_Group=c(outl_dist))
## ggplot
png("Figures_Tables/full_res_graphs/supplemental_chisq_pairwise.png", width=7, height=7, units="in", res=700)
ggpairs(dist_df2, columns=1:4, ggplot2::aes(col=Color_Within_Group, fill=Color_Within_Group),
        lower=list(continuous = wrap("points", alpha=0.1, pch=16))) +
  scale_color_manual(values=c("red", "grey")) +
  scale_fill_manual(values=c("red", "grey")) +
  theme_bw() +
  theme(title = element_text(face="bold"),
        panel.grid = element_blank())
dev.off()

# Seems that only some of the outliers are the problem for Euclidean - which ones?
##    Let's try using a different subset of outliers (choosing just the worst one)
tmp$Outlying <- 0
tmp[c("T.CS.149"),"Outlying"] <- 1
outl_dist <- dist(as.numeric(tmp$Outlying))
outl_dist <- as.factor(c(outl_dist) < 1)
levels(outl_dist) <- c("outlier", "typical")
dist_df3 <- data.frame(BrayCurtis=c(bray), UnweightedUniFrac=c(u_unifrac), Euclidean=c(euclid),
                       ChiSquared=c(chisq), Color_Within_Group=c(outl_dist))
png("Figures_Tables/full_res_graphs/supplemental_euclid_pairwise.png", width=7, height=7, units="in", res=700)
ggpairs(dist_df3, columns=1:4, ggplot2::aes(col=Color_Within_Group, fill=Color_Within_Group),
        lower=list(continuous = wrap("points", alpha=0.1, pch=16))) +
  scale_color_manual(values=c("red", "grey")) +
  scale_fill_manual(values=c("red", "grey")) +
  theme_bw() +
  theme(title = element_text(face="bold"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90))
dev.off()
## seems the outlier for Euclidean is a single sample: T.CS.149



##### SUPP FIG 1 & 2 EXPERIMENT: Why the Chi-Sq outliers? Why the one Euclidean outlier?? #####

# Removing these outliers makes a more normal looking plot
outliers <- c("CS.253", "CS.378", "CS.283", "T.CS.149", "CS.339")
idx_out <- which(colnames(rare_otu) %in% outliers)
chisq2 <- vegdist(t(rare_otu[,-idx_out]), method="chisq")
w_plt <- pcoa_plot(chisq, k=2, meta, "Chi-square")
wo_plt <- pcoa_plot(chisq2, k=2, meta[-idx_out,], "Chi-Square w/o outliers")  # more reasonable w/o outliers
png('Figures_Tables/full_res_graphs/supplemental_chisq_compare.png', width=7, height=4, res=700, units="in")
grid.arrange(w_plt, wo_plt, nrow=1)
dev.off()

# What is different about these ChiSq outliers?
## --- empty OTUs : outliers appear to be in the top quantile... but not by much
colMeans(as.matrix(rare_otu)[,idx_out] == 0)
quantile(colMeans(as.matrix(rare_otu)[,-idx_out] == 0))
tmp <- data.frame(value=colMeans(as.matrix(rare_otu) == 0), col_outl=c("typical","outlier")[(colnames(rare_otu) %in% outliers)+1])
png("Figures_Tables/full_res_graphs/supplemental_sparsity_outliers.png", width=5, height=5, res=700, units="in")
ggplot(tmp, aes(x=col_outl, y=value, col=col_outl)) + 
  geom_boxplot() +
  geom_jitter(pch=16, size=3) +
  scale_color_manual(values = c("red", "grey")) + 
  labs(x="", y="Percentage of Missing OTUs (value of 0)", title="Sparsity of Samples") +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(face="bold"),
        panel.grid = element_blank())
dev.off()
## --- large number of rare OTUs in other samples?
## ----> find the predominant species in the outliers (choosing top 20 for now)
dom_otu <- sort(rowSums(rare_otu[,outliers]), decreasing=T)[1:20]
dom_otu <- rare_otu[names(dom_otu),]
dom_otu$OTU <- rownames(dom_otu)
library(reshape)
dom_otu <- melt(dom_otu[1:5,], id=c("OTU"))
dom_otu$group <- c("typical", "outlier")[(dom_otu$variable %in% outliers)+1]
dom_otu$group_OTU <- paste0("otu ", dom_otu$OTU, " ", dom_otu$group)
dom_ord <- paste0("otu ", rep(unique(dom_otu$OTU), each=2), rep(c(" typical", " outlier"), 5))
png("Figures_Tables/full_res_graphs/supplemental_domOTUs_outliers.png", width=5, height=5, res=700, units="in")
ggplot(dom_otu, aes(x=factor(group_OTU, levels=dom_ord), y=value, col=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(pch=16, width = 0.1, size=4, alpha=0.8) +
  scale_color_manual(values=c("red", "grey")) +
  labs(x="OTUs", y="Value per sample", title="Dominant OTUs in Outliers") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, vjust=1.2, hjust=1.2),
        panel.grid = element_blank(),
        title = element_text(face="bold"))
dev.off()

## ----> compare to plot of gap between first and second most abundant samples within each species
# imp_mat <- as.matrix(rare_otu)
# gap_otu <- rep(0, nrow(imp_mat))
# names(gap_otu) <- rownames(imp_mat)
# for (s in rownames(imp_mat)) {
#   top <- sort(imp_mat[s,], decreasing=T)[1:2]
#   diff <- unname(abs(top[1] - top[2]))
#   gap_otu[s] <- diff
# }
# gap_otu <- sort(gap_otu, decreasing=T)
# gap_otu <- rare_otu[names(gap_otu),]
# gap_otu$OTU <- rownames(gap_otu)
# gap_otu <- melt(gap_otu[1:10,], id=c("OTU"))
# gap_otu$group <- c("typical", "outlier")[(gap_otu$variable %in% outliers)+1]
# gap_otu$group_OTU <- paste0("otu ", gap_otu$OTU, " ", gap_otu$group)
# gap_ord <- paste0("otu ", rep(unique(gap_otu$OTU), each=2), rep(c(" typical", " outlier"), 5))
# ggplot(gap_otu, aes(x=factor(group_OTU, levels=gap_ord), y=value, col=group)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(pch=16, width = 0.1, size=3) +
#   scale_color_manual(values=c("red", "grey")) +
#   labs(x="OTUs", y="Value per sample", title="Outlier samples within OTUs") +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle=45, vjust=1.2, hjust=1.2),
#         panel.grid = element_blank(),
#         title = element_text(face="bold"))
## ----> still not the answer... looks like similar gaps happen but aren't outliers
## since these only appear in chisq and euclidean, the problem lies in some unshared OTUs
# num_unshared <- rep(0, ncol(rare_otu)) # pairwise number of unshared OTUs
# for (i in 1:ncol(rare_otu)) {
#   for (j in 1:ncol(rare_otu)) {
#     unsh <- sum((rare_otu[,i] > 0) & (rare_otu[,j] == 0))
#     totalij <- sum(colSums(rare_otu[,c(i,j)]) > 0)
#     num_unshared[i] <- num_unshared[i] + (unsh / totalij)
#   }
# }
# names(num_unshared) <- colnames(rare_otu)
# num_unshared <- sort(num_unshared, decreasing=T)
# head(num_unshared, 10)
# ## ----> why is it that Euclidean doesn't have a problem in PCoA?
# tmp <- cmdscale(chisq, k=2)
# plot(tmp)
# tmp2 <- cmdscale(euclid, k=2)
# plot(tmp2)
# cors <- apply(rare_otu, 1, FUN = function (x) { cor(x, tmp[,2], method="pearson") })
# cors <- cors[order(abs(cors), decreasing = T)]
## ----> top corrs w/ PC2: OTU 4481, 5476, 6156, 4453, 6824, 515, 4287, 1431, 1314
## ----> but not all of these are disproportionate for the outliers, only 4 of these 10.
## ----> chisq is based on probabilities, how could this be the heart of our problem?

## ----> try chisq with NMDS
# chi_nmds <- metaMDS(chisq, k=3, try=20, trymax=50, maxit=500)
# plot_2d(chi_nmds$points, meta, title_text="Chi-Sq with NMDS") # could not converge! Looks bad too
# chi_nmds <- metaMDS(chisq2, k=3, try=20, trymax=50, maxit=500)
# plot_2d(chi_nmds$points, meta[-idx_out,], title_text="Chi-Sq with NMDS (no outliers)") # still bad!
## --- Evenness : seems to be in lower 2 quartiles of pielou/bulla evenness
# microbiome::evenness(rare_otu[,idx_out], 'bulla')
# quantile(microbiome::evenness(rare_otu[,-idx_out], 'bulla')[,1])

print('Done with supplemental figures.')


