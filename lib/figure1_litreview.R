##### FIG 1 LITERATURE REVIEW (final formatting done in Adobe Illustrator) #####
# Depends on being run within driver script 'beta_review_figs.R'

# Loading data
zotero <- read.csv("data/zotero/reviewed_subset_zotero_07-20-2022.csv", header=T)
zotero <- zotero[,c("Title","Author","Publication.Year","Date",
                    "Publication.Title","DOI","Item.Type","Manual.Tags")]
colnames(zotero) <- c("Title", "Author","Publication.Year","Date","Journal","DOI","Item.Type","Tags")

# Simplify Tags
tags <- strsplit(zotero$Tags, "; ") # creates a list of the tags split
zotero$Sequencing <- sapply(tags, function (x) { paste0(x[grepl("sequencing", x)], collapse="; ") })
zotero$Biome <- sapply(tags, function (x) { paste0(gsub("biome-", "", x[grepl("biome-", x)]), collapse="; ") })
zotero$Variables <- sapply(tags, function (x) { paste0(gsub("var-", "", x[grepl("var-", x)]), collapse="; ") })
zotero$AlphaAnalysis <- sapply(tags, function (x) { paste0(gsub("alpha-analysis-", "", x[grepl("alpha-analysis-", x)]), collapse="; ") })
zotero$AlphaDiv <- sapply(tags, function (x) { paste0(gsub("alpha-[^a]", "", x[grepl("alpha-", x)]), collapse="; ") })
zotero$BetaDist <- sapply(tags, function (x) { paste0(gsub("dist-", "", x[grepl("dist-", x)]), collapse="; ") })
zotero$BetaOrd <- sapply(tags, function (x) { paste0(gsub("ordin-", "", x[grepl("ordin-", x)]), collapse="; ") })
zotero$BetaAnalysis <- sapply(tags, function (x) { paste0(gsub("beta-analysis-", "", x[grepl("beta-analysis-", x)]), collapse="; ") })
zotero$FurtherAnalysis <- sapply(tags, function (x) { paste0(gsub("further-analysis-", "", x[grepl("further-analysis-", x)]), collapse="; ") })
zotero$Arch <- sapply(tags, function (x) { paste0(gsub("arch-", "", x[grepl("arch-", x)]), collapse="; ") })
zotero$Study <- sapply(tags, function (x) { paste0(gsub("study-", "", x[grepl("study-", x)]), collapse="; ") })

# Biomes Distribution
table(zotero$Biome)

# Distance/Dissimilarity Measures Distribution
table(zotero$BetaDist)

# Ordination Distribution
table(zotero$BetaOrd)

# Statistical tests Distribution
table(zotero$BetaAnalysis)
simple_ba <- strsplit(zotero$BetaAnalysis, "; ")
unique(unlist(simple_ba))
## simplify these further since we have so many tests
simple_ba <- lapply(simple_ba, function (x) { ifelse(x %in% c("PERMANOVA", "ANOSIM"), "PERMANOVA/ANOSIM", x) })
simple_ba <- lapply(simple_ba, function (x) { ifelse(x %in% c("Kruskal-Wallis test","Student t-test","Wilcoxon test",
                                                              "Welch's t-test","Spearman correlation","linear regression model",
                                                              "logistic regression", "LASSO conditional logistic regression","ANCOM",
                                                              "ANOVA","ANCOVA","Bartlett's test","GLM","LDM","Tukey test",
                                                              "Fisher's exact test","LSD"), "per taxon tests", x) })
simple_ba <- lapply(simple_ba, function (x) { ifelse(x %in% c("LDA", "PLS-DA"), "discriminant analysis", x) })
simple_ba <- lapply(simple_ba, function (x) { ifelse(x %in% c("UPGMA","hierarchical clustering","PAM","Mantel test"), "clustering", x) })
simple_ba <- lapply(simple_ba, function (x) { ifelse(x %in% c("LCBD","db-RDA","SIMPROF","distance decay",
                                                              "Intraclass correlation coefficient (ICC)","beta MNTD",
                                                              "Random Forest","F/B ratio","pRDA","Procrustes","MiRKAT",
                                                              "OTU overlap"), "other", x) })
zotero$BetaAnalysisSimple <- sapply(simple_ba, function (x) { paste0(unique(x), collapse="; ") })
sort(table(unlist(simple_ba)), decreasing=T)
tmp <- lapply(simple_ba, function (x) { ifelse(x %in% c("clustering", "dispersion tests"), "other", x) })
tmp <- sapply(tmp, function (x) { paste0(sort(unique(x)), collapse="; ") })
table(tmp)[order(names(table(tmp)))]

# Arch Artifacts Distribution
table(zotero$Arch)

print('Done with figure 1.')


