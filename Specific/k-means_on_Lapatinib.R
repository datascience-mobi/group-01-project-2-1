# Performing a k-means on Lapatinib
library(cluster)
library(dplyr)


LapatinibFold = select(Fold_Change, contains("Lapa"))

#Determining the number of clusters
topVarFold = apply(LapatinibFold, 1, var)
summary(topVarFold)

# Using the most variable, thus informative genes
topVarFold75 = LapatinibFold[topVarFold > quantile(topVarFold, probs = 0.75), ]
dim(topVarFold75)

km = kmeans(x = t(topVarFold75), centers = 2, nstart = 10)
km$tot.withinss

km = kmeans(x = t(topVarFold75), centers = 3, nstart = 10)
km$tot.withinss

#running a loop for the best n (searching for "ellbow")
wss = as.data.frame(sapply(2:7, function(k) {
  kmeans(x =t(topVarFold75), centers = k)$tot.withinss})  )            
plot(2:7, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares", main = "Determining the amount of clusters from Foldchange")

#no real ellbow :(


# Using the silhouett method
D = dist(t(topVarFold75))
km = kmeans(x = t(topVarFold75), centers = 3, nstart = 10)
s = silhouette(km$cluster, D)
plot(s)

#plot kmeans
library(factoextra)
library(fpc)
clus <- kmeans(topVarFold75, centers=3)
clusplot(topVarFold75, clus$cluster, color=TRUE, shade=TRUE, labels=2, lines=0, col.clus = c("purple", "pink", "lightseagreen"), main = "K-means Plot", col.p = c("purple", "pink", "lightseagreen") )
plotcluster(topVarFold75, clus$cluster)

#dendrogram
cor.mat = cor(topVarFold75, method = "spearman")
cor.dist = as.dist(1 - cor.mat)
cor.hc = hclust(cor.dist, method = "ward.D2")
cor.hc = as.dendrogram(cor.hc)
library(dendextend)
cor.hc %>% set("labels_col", value = c("purple", "pink", "lightseagreen"), k=3) %>%  plot(main = "Color labels \n by cluster", horiz = TRUE)


#PCA
pca = prcomp(topVarFold75, center = T, scale. = T)  
print(pca)
plot(pca, type = "l") #First two componets explain most of the variability in the data
plot(pca$rotation[, 1], pca$rotation[, 2], col = c("maroon1", "turquoise3"), pch = 19, xlab = "PC1", ylab = "PC2")

