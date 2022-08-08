##### FIG 4: ORDINATION METHOD FIGURES #####
# Depends on being run within driver script 'beta_review_figs.R'
# NOTE: Flipping Axes to match original paper Fig 2A (https://doi.org/10.1016/j.cell.2018.10.029)

# FIG 4 PART A : Compare Vis Methods visually : Unweighted UniFrac with all methods
## PCoA - use pcuuf from above
pc <- cmdscale(u_unifrac, k=3, eig=T)
tmp <- pc$points
tmp[,1] <- tmp[,1]*-1
tmp[,2] <- tmp[,2]*-1
pc_plt <- plot_2d(tmp, meta, "PCoA", pc$eig)
## NMDS
nmds <- metaMDS(u_unifrac, k=3, try=20, trymax=50, maxit=500)
tmp <- nmds$points
tmp[,1] <- tmp[,1]*-1
tmp[,2] <- tmp[,2]*-1
nmds_plt <- plot_2d(tmp, meta, "NMDS")
## t-SNE
tsne <- Rtsne(u_unifrac, is_distance=T, dims=3)
tmp <- tsne$Y
tmp[,2] <- tmp[,2]*-1
tsne_plt <- plot_2d(tmp, meta, "t-SNE", legendpos = c(0.21,0.8))
## save a png with this figure
png("Figures_Tables/full_res_graphs/ordination_compare.png", width=8, height=8, res=700, units="in")
grid.arrange(pc_plt, nmds_plt, tsne_plt, nrow=2)
dev.off()


# FIG 4 PART B : Procrustes analysis (2D plots, but using as many dimensions as possible)
## Procrustes analysis (using protest to get 'significance' of Procrustes statistic, Peres-Neto & Jackson (2001))
## how protest works:
##       permutes this function: procr <- function(X, Y) sum(svd(crossprod(X, Y), nv=0, nu=0)$d)
##            which essentially finds orthogonal space between our scaled matrices and computes SVD singular values (in d) from this space
##       compare sum of singular vlaues to the sum of suares from procrustes : perm >= ( sqrt(1 - sol$ss) ) - EPS
##            if the sum of singular values is larger than sqrt(1-sum of squares), then significance score increases
##       final signif score is the percentage of permutations which "pass" : i.e. svd singular values of orthogonal matrix exceed S.o.S. of original matrices
pcoa_nmds <- protest(pc$points, nmds$points, scale=T, scores = "sites")
pcoa_tsne <- protest(pc$points[,1:3], tsne$Y, scale=T, scores = "sites")
nmds_tsne <- protest(nmds$points[,1:3], tsne$Y, scale=T, scores = "sites")

## Plots for Procrustes
## --- PCoA and NMDS
p_n <- proc_plot(pcoa_nmds, meta, "Procrustes: PCoA & NMDS")
## --- PCoA and t-SNE
p_t <- proc_plot(pcoa_tsne, meta, "Procrustes: PCoA & t-SNE")
## --- PCoA and CA
n_t <- proc_plot(nmds_tsne, meta, "Procrustes: NMDS & t-SNE")
## --- Joined plot
png("Figures_Tables/full_res_graphs/ordination_procrustes.png", width=12, height=4, res=700, units="in")
grid.arrange(p_n, p_t, n_t, nrow=1)
dev.off()

print('Done with figure 4.')


