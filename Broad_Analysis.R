# Getting the Data
library(rstudioapi)
library(clusters) #meinst du vielleicht "library(cluster)"? das package "clusters" gibt es nicht ;)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

NCI_TPW_gep_treated = readRDS(paste0(wd, "/Data/NCI_TPW_gep_treated.rds"))
NCI_TPW_gep_untreated = readRDS(paste0(wd, "/Data/NCI_TPW_gep_untreated.rds"))
Metadata = read.delim(paste0(wd, "/Data/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t", stringsAsFactors = TRUE)
Cellline_Annotation = read.delim(paste0(wd, "/Data/cellline_annotation.tsv"), header = TRUE, sep = "\t", stringsAsFactors = TRUE)
Drug_Annotation = read.delim(paste0(wd, "/Data/drug_annotation.tsv"), header = TRUE, sep = "\t", stringsAsFactors = TRUE)
CCLE_mutations = readRDS(paste0(wd, "/Data/CCLE_mutations.rds"))
CCLE_copynumber = readRDS(paste0(wd, "/Data/CCLE_copynumber.rds"))
CCLE_basalexpression = readRDS(paste0(wd, "/Data/CCLE_basalexpression.rds"))
NegLogGI50 = as.data.frame(readRDS(paste0(wd, "/Data/NegLogGI50.rds")))
Fold_Change = NCI_TPW_gep_treated - NCI_TPW_gep_untreated
Treated = data.frame(NCI_TPW_gep_treated)
Untreated = data.frame(NCI_TPW_gep_untreated)
Fold_Change <- data.frame(Fold_Change)

# Sort out for Lapatinib

#dat = NCI_TPW_gep_treated[, 421:474]
#LapatimibTreat <- data.frame(dat)
#rm(dat)
#dat = NCI_TPW_gep_untreated[, 421:474]
#LapatimibUnTreat <- data.frame(dat)
#rm(dat)

#Boxplots,check for normalization ; scale the Data', FC takes long

boxplot(Treated)
boxplot(Untreated)
boxplot(Fold_Change)
Treated = scale(Treated)
Untreated = scale(Untreated)
Fold_Change = scale(Fold_Change)
boxplot(Treated, ylab = "Gene expression profile", main = "Teated genexpressionprofiles",xaxt = "n")
boxplot(Untreated, ylab = "Gene expression profile", main = "Untreated genexpressionprofiles",xaxt = "n")
boxplot(Fold_Change, ylab = "Gene expression profile", main = "Fold_Change genexpressionprofiles",xaxt = "n")
boxplot(Treated[,1:10], ylab = "Gene expression profile", main = "First 10 reated genexpressionprofiles")


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

#perfoming PCA
pca = prcomp(topVarTreat75, center = F, scale. = F)
plot(pca, type = "l")
plot(pca$rotation[, 1], pca$rotation[, 2], col = c("blue","red"), xlab = "PC1", 
ylab = "PC2")
