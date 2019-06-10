# Performing a k-means on Lapatinib
library(cluster)

Fold_Change = data.frame(scale(Fold_Change))
LapatinibFold = select(Fold_Change, contains("Lapa"))

#Determining the number of clusters
topVarFold = apply(LapatinibFold, 1, var)
summary(topVarTreated)

# Using the most variable, thus informative genes
topVarFold75 = LapatinibFold[topVarFold > quantile(topVarFold, probs = 0.75), ]
dim(topVarFold75)

km = kmeans(x = t(topVarFold75), centers = 2, nstart = 10)
km$tot.withinss

km = kmeans(x = t(topVarFold75), centers = 3, nstart = 10)
km$tot.withinss

#running a loop for the best n (searching for "ellbow")
wss = sapply(2:7, function(k) {
  kmeans(x =t(topVarFold75), centers = k)$tot.withinss})            
plot(2:7, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares", main = "Determining the amount of clusters from Foldchange")

#no real ellbow :(


# Using the silhouett method
D = dist(t(topVarFold75))
km = kmeans(x = t(topVarFold75), centers = 3, nstart = 10)
s = silhouette(km$cluster, D)
plot(s)


