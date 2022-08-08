##### FIG 6: ARTIFACTS #####
# Depends on being run within driver script 'beta_review_figs.R'

# PART A: Sampling Frame (side by side)
## Make a random 2000 x 200 matrix
x <- matrix(rnorm(12000),nr=600)
## Increase values in the first 10 features in the
##     first 200 rows so that they form a separate cluster
x[1:200,1:10] <- x[1:200,1:10] + rnorm(2000,1.5,.25)
## Increase values in the last 10 features in the
##     middle 200 rows so that they form yet another separate cluster
x[201:400,11:20] <- x[201:400,11:20] + rnorm(2000,1.5,.25)
## Create IDs of the three clusters
group.id <- rep(1:3,times=c(200,200,200))
sf_cols <- alpha(c('#530070','#00948b','#faeb00')[group.id], 0.7) # using my color scheme (viridis)
sf_pch <- c(21,22,24)[group.id]
## Visualize all points (saved as png)
png('Figures_Tables/full_res_graphs/sampling_frame_example.png', width=9, height=5, res=700, units="in")
par(mfrow=c(1,2))
ix <- 1:600
sf_pc <- cmdscale(dist(x[ix,]))
plot(sf_pc[,1], sf_pc[,2], bg=sf_cols[ix], col='#00000077', pch=sf_pch[ix], main='200 points per group', xlab='PC 1', ylab='PC 2')
legend('topleft',c('Cluster 1','Cluster 2','Cluster 3'), pch=c(21,22,24), pt.bg=c('#530070','#00948b','#faeb00'),col='#00000077')
## Visualize with one group diminished
nkeep <- 4
ix <- (200-nkeep+1):600
sfr_pc <- cmdscale(dist(x[ix,]))
plot(sfr_pc[,1], sfr_pc[,2], bg=sf_cols[ix], col='#00000077', pch=sf_pch[ix], main='Same data, 4 Cluster 1 points', xlab='PC 1', ylab='PC 2')
dev.off()
par(mfrow=c(1,1))


# PART B: Arch in Soil Data
## Colors & shapes from original paper
soil_c <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "black") # colors (borders)
soil_bg <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white") # fills (background)
soil_s <- c(22, 25, 24, 21, 23, 23) # shapes
ord <- cut(meta_soil$ph, breaks = c(0,4:8,14)) # pH groups (<4, 4-5, 5-6, 6-7, 7-8, >8)
soil_c <- soil_c[ord]; soil_bg <- soil_bg[ord]; soil_s <- soil_s[ord];
## Compute PCoA (note we are using Jaccard instead of unweighted UniFrac because no tree available)
dist_soil  <- vegdist(t(soil_otu), method="jaccard")
pc_soil <- cmdscale(dist_soil, k=2, eig=T)
pc_soil$points[,1] <- -1 * pc_soil$points[,1]
pc_soil$points[,2] <- -1 * pc_soil$points[,2]
## Plotting
calc.perc.var <- function (eigen, dimension) {
  percents <- round((eigen/sum(eigen))*100, 1) # rounds percentages to one decimal place
  return(percents[dimension])
}
png("Figures_Tables/full_res_graphs/arch_soil_pcoa.png", width=6, height=6, res=700, units="in")
par(mgp=c(0.8, 0, 0))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pc_soil$points, xaxt='n', yaxt='n',
     xlim=range(pc_soil$points)+c(-0.1, 0.1), ylim=range(pc_soil$points)+c(-0.1, 0.1),
     col=soil_c, bg=soil_bg, pch=soil_s, cex=1.2, lwd=3,
     xlab=paste0("PC1 [", calc.perc.var(pc_soil$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(pc_soil$eig,2), "%]"), cex.lab=2.4)
legend("topright", legend=c("<4", "4-5", "5-6", "6-7", "7-8", ">8"), pch=soil_s,
       pt.bg=soil_bg, col=soil_c, inset=c(-0.3,0), title="pH")
dev.off()
## PCoA of soil with unweighted UniFrac distances (as used in the paper)


# # PART B (?): Arch in global gut study (?)
# ## Uneven sampling of the age gradient... downsampling under age 2 & ages 14-16
# hist(meta_gg$age[meta_gg$geo_loc_name == "USA"], breaks=40, main="Histogram of USA age distribution", xlab="age")
# u2_keep <- sample(rownames(meta_gg)[meta_gg$age < 2 & meta_gg$geo_loc_name == "USA"], size=15)
# teen_keep <- sample(rownames(meta_gg)[meta_gg$age >= 14 & meta_gg$age <= 16 & meta_gg$geo_loc_name == "USA"], size=30)
# other_keep <- rownames(meta_gg)[meta_gg$age > 2 & meta_gg$age < 14 & meta_gg$geo_loc_name == "USA"]
# other_keep <- c(other_keep, rownames(meta_gg)[meta_gg$age > 16 & meta_gg$age <= 30 & meta_gg$geo_loc_name == "USA"])
# ## Restrict the data sets
# gg_keep <- c(u2_keep, teen_keep, other_keep)
# meta_gg <- meta_gg[gg_keep,]
# temp_mat <- as.matrix(gg_uunifrac)
# temp_mat <- temp_mat[gg_keep, gg_keep]
# gg_uunifrac <- as.dist(temp_mat)
# ## PCoA of global gut, set up dataframe
# gg_pc <- cmdscale(gg_uunifrac, k=2, eig=T)
# gg_var <- (gg_pc$eig/sum(gg_pc$eig))[1:2]
# gg_var <- round(gg_var*100, 1)
# gg_pc <- as.data.frame(gg_pc$points)
# colnames(gg_pc) <- paste0("PC", 1:2)
# identical(meta_gg$Sample_ID, rownames(gg_pc))
# gg_pc$age <- meta_gg$age
# ## ggplot
# ggplot(gg_pc, aes(x=PC1, y=PC2, col=age)) +
#   geom_point(size=4, alpha=0.7) +
#   scale_color_gradient(low="#ffd700", high="#0000ff") +
#   labs(title="", x=paste0("PC 1 [", gg_var[1], "%]"),
#        y=paste0("PC 2 [", gg_var[2], "%]"), col="age") +
#   theme_bw() + 
#   theme(panel.grid = element_blank(),
#         text = element_text(family="Helvetica", size=16))

print('Done with figure 6.')


