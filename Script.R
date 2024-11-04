pkgs_needed <- c("tidyverse", "factoextra", "ggpubr")
for (pkg in pkgs_needed) {
  if(!require(pkg, character.only = T)) library(pkg, character.only = T)
}
data("iris")
dim(iris)
head(iris)
str(iris)
pA <- iris |> 
  pivot_longer(!Species) |> 
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  theme_pubr(base_size = 9) +
  theme(axis.title.x = element_blank()) + 
  rotate_x_text(angle = 45) +
  ggtitle("Raw data")
pB <- 
  as.data.frame(scale(iris[,1:4], center = T, scale = F)) %>%
  pivot_longer(names(iris)[1:4]) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  theme_pubr(base_size = 9) +
  theme(axis.title.x = element_blank()) +
  rotate_x_text(angle = 45) +
  ggtitle("Centred data")
pC <- 
  as.data.frame(scale(iris[,1:4], center = T, scale = T)) %>%
  pivot_longer(names(iris)[1:4]) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  theme_pubr(base_size = 9) +
  theme(axis.title.x = element_blank()) +
  rotate_x_text(angle = 45) +
  ggtitle("Scaled & Centred data")
ggarrange(pA, pB, pC, ncol = 3, labels = "AUTO")

iris |> 
  select(!Species) |> 
  scale(center = T, scale = T) |> 
  cor()

iris |> 
  ggplot(aes(x = scale(Petal.Length), y = scale(Petal.Width))) +
  geom_point(aes(color = Species)) +
  geom_smooth(method = "lm") + # adds line of best fit defined by linear regression
  theme_classic()

iris.svd <- iris |> 
  select(!Species) |> 
  prcomp(scale = F)

summary(iris.svd)

PCs_sd <- apply(iris.svd$x,2,sd)
PCs_sd

PCs_sd^2/sum(PCs_sd^2)

cumsum(PCs_sd^2/sum(PCs_sd^2))

raw_sds <- apply(iris[,-5], 2, sd)
raw_sds_prop <- raw_sds^2/sum(raw_sds^2)
cumsum(raw_sds_prop[order(-raw_sds_prop)])

plot(cumsum(PCs_sd^2/sum(PCs_sd^2)), type = "o", lty = 1, ylim = c(0.5,1),
     ylab = "Cumulative Proportion of Variance", xlab = "Dimensions")
points(cumsum(raw_sds_prop[order(-raw_sds_prop)]), col = "red", pch = "*")
lines(cumsum(raw_sds_prop[order(-raw_sds_prop)]), col = "red", lty = 2)
legend(3,0.7, legend = c("PCA", "raw"), col = c("black", "red"), lty = c(1,2))

barplot(PCs_sd^2/sum(PCs_sd^2), ylab = "Proportion of Variance")

fviz_eig(iris.svd, addlabels = T) +
  theme_pubr(base_size = 9)

str(iris.svd)

fviz_pca_var(iris.svd) + 
  xlim(c(-0.5,2)) + ylim(c(-0.5,0.2))

fviz_pca_var(iris.svd, select.var = list(cos2 = 2)) + xlim(c(-0.5,2)) + ylim(c(-0.5,0.2))

dim(iris.svd$x)

PC_scores <- scale(iris[,-5], center = T, scale = F) %*% iris.svd$rotation

head(PC_scores)
table(PC_scores == iris.svd$x)

fviz_pca_ind(iris.svd, habillage = iris$Species, label = "none", addEllipses = T, invisible = "quali") +
  theme_pubr(base_size = 9)

fviz_pca_biplot(iris.svd, habillage = iris$Species, label = "var", addEllipses = T, invisible = "quali") +
  theme_pubr(base_size = 9) +
  expand_limits(x = 6)

###GSE48035
#load data
load("DS1.Rdata")
#run PCA
DS1.svd <- DS1_lcpm.filtered.norm |> 
  t() |> #important to transpose data frame e.g. samples must be rows
  prcomp(scale = F)
#screeplot
fviz_eig(DS1.svd, labels = T)
#visualise data in 2d
fviz_pca_ind(DS1.svd, repel = F)
#visualise sample meta data
fviz_pca_ind(DS1.svd, habillage = DS1_sample.metaData$Sample, label = "none", addEllipses = T, invisible = "quali") +
  theme_pubr(base_size = 9)

fviz_pca_ind(DS1.svd, habillage = DS1_sample.metaData$Tissue.type, label = "none", addEllipses = T, invisible = "quali") +
  theme_pubr(base_size = 9)

fviz_pca_ind(DS1.svd, habillage = DS1_sample.metaData$Protocol, label = "none", addEllipses = T, invisible = "quali") +
  theme_pubr(base_size = 9)

#change row names from gene_ID to row names
gene.names.idx <- match(rownames(DS1_lcpm.filtered.norm), DS1_gene.metaData$gene_id)

gene.names <- DS1_gene.metaData$gene_name[gene.names.idx]

gene.names.unique <- make.unique(gene.names)

rownames(DS1.svd$rotation) <- gene.names.unique

#plot feature coordinates
fviz_pca_var (DS1.svd, select.var = list(cos2 = 20), repel = T)

#find top 10 genes with most variance
as.data.frame(DS1.svd$rotation[,1:2]) |> 
  slice_max(n=10, order_by = PC1)

#create a biplot
fviz_pca_biplot(DS1.svd, select.var = list(name = c("GFAP", "MKI67")),
                habillage = DS1_sample.metaData$Tissue.type,
                label = "var" , 
                invisible = "quali",
                addEllipses = T) + ggtitle("Tissue Type")
