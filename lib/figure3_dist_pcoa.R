##### FIG 3 DISTANCE PCoA PLOTS:  PCoA Comparison #####
# Depends on being run within driver script 'beta_review_figs.R'

# Generate PCoA plots
pcbray <- pcoa_plot(bray, k=2, meta, "Bray-Curtis", adjx=-1, adjy=-1)
pckulc <- pcoa_plot(kulcz, k=2, meta, "Kulczynski", adjx=-1, adjy=-1)
pcmanh <- pcoa_plot(manh, k=2, meta, "Manhattan", adjx=-1, adjy=-1)
pccos <- pcoa_plot(cos, k=2, meta, "Cosine", adjx=-1, adjy=-1)
pcjacc <- pcoa_plot(jacc, k=2, meta, "Jaccard", adjx=-1, adjy=-1)
pcuuf <- pcoa_plot(u_unifrac, k=2, meta, "Unweighted UniFrac", adjx=-1, adjy=-1)
pcwuf <- pcoa_plot(w_unifrac, k=2, meta, "Weighted UniFrac", adjx=-1)
pceuc <- pcoa_plot(euclid, k=2, meta, "Euclidean", adjx=-1, adjy=-1)
pcait1 <- pcoa_plot(aitch, k=2, meta, "Aitchison", adjx=-1, adjy=-1)
pcait2 <- pcoa_plot(robust_a, k=2, meta, "Robust Aitchison", adjx=-1, adjy=-1)
pcchi <- pcoa_plot(chisq, k=2, meta, "Chi-Square", adjx=-1)
pclegend <- plot_imp_legend(meta)
# Print out for publication
png("Figures_Tables/full_res_graphs/dist_ordination.png", width=13, height=9, res=700, units="in")
grid.arrange(pcbray, pckulc, pcmanh, pccos, pcjacc, pcuuf,
             pcwuf, pceuc, pcait1, pcait2, pcchi, pclegend, nrow=3)
dev.off()

print('Done with figure 3.')

