# Getting the Data
library(rstudioapi)
library(cluster) 



#Boxplots,check for normalization ; scale the Data',

boxplot(Treated)
boxplot(Untreated)
boxplot(Fold_Change)
list = list(Treated,Untreated,Fold_Change)
nlist = lapply(list,scale)
Treated = as.data.frame(nlist[[1]])
Untreated = as.data.frame(nlist[[2]])
Fold_Change = as.data.frame(nlist[[3]])
boxplot(Treated, col= df, ylab = "Gene expression profile", main = "Teated genexpressionprofiles",xaxt = "n")
boxplot(Untreated, ylab = "Gene expression profile", main = "Untreated genexpressionprofiles",xaxt = "n")
boxplot(Fold_Change, ylab = "Gene expression profile", main = "Fold_Change genexpressionprofiles",xaxt = "n")
boxplot(Treated[,1:10], ylab = "Gene expression profile", main = "First 10 reated genexpressionprofiles")

#couloured boxplots, creating a susbet of Metadata based on the amount of tested genes in Treated
Treated1 = readRDS(paste0(wd, "/Data/NCI_TPW_gep_treated.rds"))
df = data.frame(t(Treated1))
df.data <- data.frame(sample = rownames(df))
adjustedMeda = subset(Metadata, sample %in% intersect(Metadata$sample, df.data$sample))
rm(df,df.data, Treated1)

palette(rainbow(15))
boxplot(Treated, border=adjustedMeda$drug, ylab = "Gene expression profile", main = "Teated genexpressionprofiles",xaxt ="n")

#Densityplot, abline show the 3 quantiles (      25%       50%       75%    )

plot(density(NCI_TPW_gep_treated), "Densityplot Treated vs Untreated")
lines(density(NCI_TPW_gep_untreated), col = "indianred2")
legend("topright", legend = c("treated", "untreated"), col = c("black", "indianred2"), pch = 20)
abline(v = quantile(NCI_TPW_gep_treated)[2:4], col = c("lightblue", "blue",  "orange"), lty = 2)

# Performing a k-means on Treated
#Determining the number of clusters
topVarTreated = apply(Treated, 1, var)
summary(topVarTreated)

# Using the most variable, thus informative genes
topVarTreated75 = Treated[topVarTreated > quantile(topVarTreated, probs = 0.75), ]
dim(topVarTreated75)

km = kmeans(x = t(topVarTreated75), centers = 2, nstart = 10)
km$tot.withinss

km = kmeans(x = t(topVarTreated75), centers = 3, nstart = 10)
km$tot.withinss

#running a loop for the best n (searching for "ellbow")
wss = sapply(2:7, function(k) {
kmeans(x = t(topVarTreated75), centers = k)$tot.withinss})            
plot(2:7, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares", main = "Determining the amount of clusters from Treated")

# Using the silhouett method
D = dist(t(topVarTreated75))
km = kmeans(x = t(topVarTreated75), centers = 10, nstart = 10)
s = silhouette(km$cluster, D)
plot(s)


