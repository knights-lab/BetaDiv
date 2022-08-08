##### FIG 5: CONSTRAINED ORDINATION #####
# Depends on being run within driver script 'beta_review_figs.R'

# Canonical Correspondence Analysis: Ethnicity, Age, BMI in IMP study
## compute CCA - compare with and without outliers causing issues in CA
cca_imp <- cca(t(rare_otu) ~ Age + BMI + Ethnicity,  data=meta)
plot(cca_imp)     # with outliers
cca_imp2 <- cca(t(rare_otu[,-idx_out]) ~ Age + BMI + Ethnicity,  data=meta[-idx_out,])
plot(cca_imp2)    # without outliers
## my plot version
png("Figures_Tables/full_res_graphs/constrained_cca.png", width=5, height=5, res=700, units="in")
plot(cca_imp2, type="n", main="Constrained CA", cex.main=1.5)
points(cca_imp2, pch=g_shp[meta$Sample.Group[-idx_out]], col=alpha(g_col[meta$Sample.Group[-idx_out]], 0.4), cex=2)
text(cca_imp2, dis="cn", cex=1)
dev.off()

# Redundancy Analysis (RDA): same variables for IMP
## compute rda (try using all samples though)
rda_imp <- rda(t(rare_otu) ~ Age + BMI + Ethnicity, data=meta)
plot(rda_imp)
## my plot version
png("Figures_Tables/full_res_graphs/constrained_rda.png", width=5, height=5, res=700, units="in")
plot(rda_imp, type="n", main="Redundancy Analysis", cex.main=1.5)
points(rda_imp, pch=g_shp[meta$Sample.Group], col=alpha(g_col[meta$Sample.Group], 0.4), cex=2)
text(rda_imp, dis="cn", cex=1)
dev.off()

print('Done with figure 5.')


