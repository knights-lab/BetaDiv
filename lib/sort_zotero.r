# Sorting Lit Review Articles
# for "Microbiome beta diversity analysis and visualization"
# Suzie Hoops and Dan Knights
# July 2022

##### Set Up #####
library(here)
set.seed(25)
# set up the current directory
curr_dir <- here()
if (grepl('Google', curr_dir) != T) {
  curr_dir <- rstudioapi::getActiveDocumentContext()$path
}
curr_dir <- gsub('/Volumes', '~', curr_dir)
curr_dir <- gsub('GoogleDrive', 'Google Drive', curr_dir)
curr_dir <- gsub('MyDrive', 'My Drive', curr_dir)
curr_dir <- gsub('sort_zotero.r', '', curr_dir)
print(paste0('Setting working directory to script location: ', curr_dir))
setwd(curr_dir)


##### Read in CSV of all literature found in search (over 500) #####
original <- read.csv("zotero_litsearch_06-27-2022.csv", header=T)
allz <- original[,c("Title","Author","Publication.Year","Date",
                    "Publication.Title","DOI","Manual.Tags")]
colnames(allz) <- c("Title", "Author","Publication.Year","Date","Journal","DOI","Tags")


##### Refine (remove out-of-scope) #####
# drop the previously-deemed "out-of-scope"
allz <- allz[!grepl("remove", allz$Tags),] # drops 19 (597 down to 578)
length(unique(allz$Journal))               # from 217 unique journals
# Create a list of lists for the tags
tags <- strsplit(allz$Tags, "; ")


##### Reorganize to get the columns we want #####
allz$Sequencing <- sapply(tags, function (x) { paste0(x[grepl("sequencing", x)], collapse="; ") })
allz$Biome <- sapply(tags, function (x) { paste0(gsub("biome-", "", x[grepl("biome-", x)]), collapse="; ") })
allz$Variables <- sapply(tags, function (x) { paste0(gsub("var-", "", x[grepl("var-", x)]), collapse="; ") })
allz$AlphaDiv <- sapply(tags, function (x) { paste0(gsub("alpha-", "", x[grepl("alpha-", x)]), collapse="; ") })
allz$BetaDist <- sapply(tags, function (x) { paste0(gsub("dist-", "", x[grepl("dist-", x)]), collapse="; ") })
allz$BetaOrd <- sapply(tags, function (x) { paste0(gsub("ordin-", "", x[grepl("ordin-", x)]), collapse="; ") })
allz$Analysis <- sapply(tags, function (x) { paste0(gsub("analysis-", "", x[grepl("analysis-", x)]), collapse="; ") })
allz$Tags <- NULL


##### Refine to random subset of 150 representative studies #####
# Select a random set
sub <- allz[sample(nrow(allz), 150),]
sub <- sub[order(as.numeric(rownames(sub))),]
subset <- sort(as.numeric(rownames(sub)), decreasing = F)
# Check how many have biomes and what they are
table(sub$Biome)             # decent spread (huamn gut 42, soil 9, water 9)
table(sub$Sequencing)        # decent spread (16S: 138, 16S & shotgun: 4, shotgun: 8)
length(unique(sub$Journal))  # 81 unique journals
# Write this list out so we can use Zotero to further analyze
subz <- original[subset,]
write.csv(subz, "subset_zotero_07-20-2022.csv", row.names = F)
write.csv(sub, "abbrev_subset_zotero_07-20-2022.csv", row.names = F)



