## Libraries ----
library(data.table)
library(psych)
library(cowplot)
library(tidyverse)
library(dendextend)
library(corrplot)
library(arm)
library(dplyr)
library(vegan)
library(factoextra)
# custom functions
source("mm2in.R")
# install.packages("colorspace")
library(colorspace)

allchems <- read.table('clipboard',header = T)
chems <- allchems %>%
  # select chemical variables I understand (!)
  dplyr::select(XLogP:nAcid)

# generate centred/scaled data
chemsCS <- chems %>%
  mutate(across(everything(), arm::rescale))

## Clustering ----
# get correlations/distances
chemCor <- cor(chemsCS, method = "spearman")
chemDist <- as.dist(1 - chemCor) 
chemClust <- hclust(chemDist, "ward.D")
# find clusters
par(mfrow = c(1, 2))
plot(find_k(chemClust))
# build dendrogram
dend <- chemClust %>%
  as.dendrogram %>%
  color_branches(k = 6) %>%
  color_labels(k = 6)

## Do full PCA ----
# do PCA
fullPCA <- principal(chemsCS, nfactor = 3)
fullPCA
# get loadings
var_names <- rownames(fullPCA$loadings)
match_pca <- match(labels(dend), var_names)
fullLoads <- data.frame(
  variable = var_names[match_pca],
  group = labels_colors(dend),
  PC1 = fullPCA$loadings[match_pca, 1],
  PC2 = fullPCA$loadings[match_pca, 2],
  PC3 = fullPCA$loadings[match_pca, 3]
)

## Subset data to key traits ----
subChems <- chems %>%
  dplyr::select(nAcid, MW, HybRatio, nHBAcc, FMF, ALogP)
subChemsCS <- subChems %>%
  mutate(across(everything(), arm::rescale))

## Redo PCA with selected traits ----
subPCA <- principal(subChemsCS, nfactor = 3)
subPCA


## Compare PCAs ----
# procrustes rotation
vegan::protest(X = fullPCA, Y = subPCA)
# compare with full PCA using pairwise correlations of PC scores
cor(fullPCA$scores, subPCA$scores, method = "spearman")

## Various synthetic scores ----
scores <- allchems %>%
  dplyr::select(my_class) %>%
  # add PC scores from full PCA
  mutate(fullPC1 = as.vector(fullPCA$scores[, 1])) %>%
  mutate(fullPC2 = as.vector(fullPCA$scores[, 2])) %>%
  mutate(fullPC3 = as.vector(fullPCA$scores[, 3])) %>%
  # add PC scores from subset PCA
  mutate(subPC1 = as.vector(subPCA$scores[, 1])) %>%
  mutate(subPC2 = as.vector(subPCA$scores[, 2])) %>%
  mutate(subPC3 = as.vector(subPCA$scores[, 3]))

## Different correlation combinations ----
STcor <- cor(dplyr::select(scores, subPC1:subPC3), dplyr::select(chems, XLogP:nAcid), method = "spearman")
FScor <- cor(dplyr::select(scores, subPC1:subPC3), dplyr::select(scores, fullPC1:fullPC3), method = "spearman")


## Reorder those matched to traits ----
STcorSort <- STcor[, match(labels(dend), colnames(STcor))]

## Tests for full and subset PC axes ----
cor.test(x = scores$subPC1, y = scores$fullPC1, method = "spearman")
cor.test(x = scores$subPC2, y = scores$fullPC3, method = "spearman")
cor.test(x = scores$subPC3, y = scores$fullPC2, method = "spearman")



## Dendrogram ----
postscript(
  file = "chem_traits_dendro.eps",
  width = mm2in(150), height = mm2in(90)
)
plot(dend)
dev.off()

## Trait-by-trait correlation plot ----
postscript(
  file = "chem_traits_corr.eps",
  width = mm2in(150), height = mm2in(150)
)
corrplot::corrplot(
  chemCor, 
  order = "hclust", 
  hclust = "ward.D", 
  col = c("white", "black"), 
  outline = T, 
  cl.pos = "n", 
  addrect = 6, 
  rect.col = colorspace::rainbow_hcl(n = 6)
)
dev.off()

## Sub PC score correlation plot ----
postscript(
  file = "chem_traits_st_corr.eps",
  width = mm2in(150), height = mm2in(60)
)
corrplot::corrplot(
  STcorSort, 
  col = c("white", "black"), 
  outline = T, 
  cl.pos = "n", 
  addrect = 6, 
  rect.col = colorspace::rainbow_hcl(n = 6)
)
dev.off()

## Full and Sub PC scores ----
postscript(
  file = "chem_traits_fs_corr.eps",
  width = mm2in(60), height = mm2in(60)
)
corrplot::corrplot(
  FScor, 
  col = c("white", "black"), 
  outline = T, 
  cl.pos = "n"
)
dev.off()



## Define metabolite families of interest ----
compTypes <- c(
  "Benzenoids",
  "Lipids",
  "Nucleosides_and_nucleotides", 
  "Organic_acids",
  "Organic_oxygen",
  "Organoheterocyclic",
  "Phenylpropanoids_and_polyketides",
  "others"
)

## All chemical properties ----
traits <- allchems %>%
  bind_cols(., chemsCS) %>%
  # make class a factor with levels in bespoke order
  mutate(my_class = factor(my_class, levels = compTypes))

## Selected chemical properties ----
chosen <- scores %>%
  select(my_class, PC1 = subPC1, PC2 = subPC2, PC3 = subPC3) %>%
  bind_cols(., subChems)



## Calculate metabolite class centroids ----
centroids <- chosen %>%
  group_by(my_class) %>%
  summarise(across(PC1:PC3, mean)) %>%
  left_join(
    chosen, ., 
    by = "my_class", 
    suffix = c(".value", ".centre")
  ) %>%
  mutate(my_class = factor(my_class, levels = compTypes))
onlyCentroids <- centroids %>%
  select(my_class, PC1.centre:PC3.centre) %>%
  unique



## Biplots ----
# base plot
baseGroupPlot <- ggplot(centroids) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(col = "none") + 
  scale_colour_brewer(palette = "Dark2") +
  geom_segment() +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)
# PC1 vs PC2
pc1pc2 <- baseGroupPlot +
  aes(
    x = PC1.value, y = PC2.value, 
    xend = PC1.centre, yend = PC2.centre, 
    col = my_class
  ) +
  scale_x_continuous(limits = c(-2, 4.5)) +
  scale_y_continuous(limits = c(NA, NA)) +
  geom_point(
    aes(x = PC1.centre, y = PC2.centre),
    data = onlyCentroids,
    shape = 21, 
    fill = "white", 
    size = 2.5
  ) +
  xlab("PC1 (34%) - size") +
  ylab("PC2 (26%) - Complexity")
pc1pc2
# PC1 vs PC3
pc1pc3 <- baseGroupPlot +
  aes(
    x = PC1.value, y = PC3.value, 
    xend = PC1.centre, yend = PC3.centre, 
    col = my_class
  ) +
  scale_x_continuous(limits = c(-2, 4.5)) +
  scale_y_continuous(limits = c(NA, 4)) +
  geom_point(
    aes(x = PC1.centre, y = PC3.centre), 
    data = onlyCentroids,
    shape = 21, 
    fill = "white", 
    size = 2.5
  ) +
  xlab("PC1 (34%) - size") +
  ylab("PC3 (23%) - A-polarity")
pc1pc3

# do boxplots of scores
baseBox <- ggplot(centroids) +
  theme_bw() +
  guides(fill = "none") +
  scale_fill_brewer(palette = "Dark2") +
  geom_boxplot(outlier.size = 0.5)
# fill them
bpPC1 <- baseBox +
  theme(panel.grid = element_blank(), axis.text.y = element_blank()) +
  aes(
    x = my_class, y = PC1.value, 
    fill = my_class
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(-2, 4.5)) +
  xlab("") +
  ylab("PC1 (34%) - size")
bpPC2 <- baseBox +
  theme(panel.grid = element_blank(), axis.text.x = element_blank()) +
  aes(
    x = my_class, y = PC2.value, 
    fill = my_class
  ) +
  scale_y_continuous(limits = c(NA, NA)) +
  xlab("") +
  ylab("PC2 (26%) - Complexity")
bpPC3 <- baseBox +
  theme(panel.grid = element_blank(), axis.text.x = element_blank()) +
  aes(
    x = my_class, y = PC3.value, 
    fill = my_class
  ) +
  scale_y_continuous(limits = c(NA, 4)) +
  xlab("") +
  ylab("PC3 (23%) - A-polarity")
# combine
postscript(
  file = "chem_pca_classes.eps",
  width = mm2in(110), height = mm2in(190)
)
plot_grid(bpPC1, NA, pc1pc2, bpPC2,pc1pc3, bpPC3,
          nrow = 3, rel_heights = c(1, 2, 2), rel_widths = c(2, 1),
          axis = "tlbr", align = "hv")
dev.off()



#3D GPSNP
library(plot3D)
library(rgl)

# Define colors for the seven types
my_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
# Ensure my_class is a factor
chem_data$my_class <- factor(chem_data$my_class)

# Map classes to colors explicitly.  This handles cases where the levels aren't in a "natural" order.
class_colors <- my_colors[as.numeric(chem_data$my_class)] # Directly use numeric values as indices



# Create the 3D scatter plot
plot3d(
  x = chem_data$PC1, 
  y = chem_data$PC2, 
  z = chem_data$PC3, 
  type = "s",            # "s" for spheres
  size = 0.5,           # Adjust sphere size
  col = class_colors,   # Use defined colors
  box = TRUE,           # Show the bounding box
  xlab = "PC1",         # X-axis label
  ylab = "PC2",         # Y-axis label
  zlab = "PC3",         # Z-axis label
  #xlim = c(-1.5, 1.5),  # Set x-axis limits - Adjust these based on your data's range
  #ylim = c(-1.5, 1.0),  # Set y-axis limits
  #zlim = c(-1.5, 1.0),  # Set z-axis limits
  axes = TRUE# Turn on the axes
)
# Add a floor to the plot with points projected onto it
floor_level <- min(chem_data$PC3) - 0.2 # Slightly below the minimum Z
points3d(chem_data$PC1, chem_data$PC2, rep(floor_level, nrow(chem_data)), 
         size = 4,
         col = class_colors)

floor_level <- min(chem_data$PC2) - 0.2  # Slightly below the minimum Y
points3d(chem_data$PC1, rep(floor_level, nrow(chem_data)), chem_data$PC3,
         size = 4,
         col = class_colors)
floor_level <- max(chem_data$PC1) + 0.2  # Slightly below the minimum X
points3d(rep(floor_level, nrow(chem_data)), chem_data$PC2, chem_data$PC3,
         size = 4,
         col = class_colors)

rgl.postscript("myplot123-raw.pdf", fmt = "pdf", drawText = T)


# Create the 3D scatter plot
plot3d(
  x = chem_data$PC2, 
  y = chem_data$PC3, 
  z = chem_data$PC4, 
  type = "s",            # "s" for spheres
  size = 0.5,           # Adjust sphere size
  col = class_colors,   # Use defined colors
  box = TRUE,           # Show the bounding box
  xlab = "PC2",         # X-axis label
  ylab = "PC3",         # Y-axis label
  zlab = "PC4",         # Z-axis label
  #xlim = c(-1.5, 1.5),  # Set x-axis limits - Adjust these based on your data's range
  #ylim = c(-1.5, 1.0),  # Set y-axis limits
  #zlim = c(-1.5, 1.0),  # Set z-axis limits
  axes = TRUE# Turn on the axes
)
# Add a floor to the plot with points projected onto it

floor_level <- min(chem_data$PC3) - 0.2  # Slightly below the minimum Y
points3d(chem_data$PC2, rep(floor_level, nrow(chem_data)), chem_data$PC4,
         size = 4,
         col = class_colors)
floor_level <- max(chem_data$PC2) + 0.2  # Slightly below the minimum X
points3d(rep(floor_level, nrow(chem_data)), chem_data$PC3, chem_data$PC4,
         size = 4,
         col = class_colors)

rgl.postscript("myplot234.pdf", fmt = "pdf", drawText = T)