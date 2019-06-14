#PCA: May 23rd 2019
#applying fold change to continous variables 
library(ggplot2)
library(gridExtra)


pca <- prcomp(t(Fold_Change), scale = TRUE)
#transform df by function t() because prcomp expects rows to be samples, which is not the case for Fold_Change
#goal: show relation between related samples

plot(pca$x[,1], pca$x[,2]) 

# get name of sample with highest pc1 value
rownames(pca$x)[order(pca$x[,1], decreasing = TRUE)[1]]


pca.var <- pca$sdev^2  # sdev calculates variation each PC accounts for
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
# since percentages make more sense then normal variation values
# calculate % or variation, which is much more interesing

barplot(pca.var.per, main = "Scree plot", xlab = "Principal Components", ylab = "% variation")
plot(pca.var.per[1:10], type = "l", xlab = "Principal Components", ylab = "% variation")
plot(pca.var.per[1:8], type = "l", xlab = "Principal Components", ylab = "% variation")
plot(cumsum(pca.var.per[1:15]), type = "l", xlab = "Principal Components", ylab = "% variation") 



pca.data <- data.frame(sample = rownames(pca$x),
                       x = pca$x[,1],
                       y = pca$x[,2]
                       )

View(pca.data)


## get names of top 10 genes that contribute most to pc1
loading_scores_1 <- pca$rotation[,1]

loading_scores <- pca$rotation[,1:6]
ranked_pca <- sort(loading_scores[,1], decreasing = TRUE)
View(ranked_pca)

gene_score <- abs(loading_scores_1) ## sort magnitude
mean(gene_score)
max(gene_score)
min(gene_score)

gene_score_ranked <- sort(gene_score, decreasing = TRUE)

head(gene_score_ranked)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes # show names of top 10 genes



pca$rotation[top_10_genes,1] ##show scores and +/- sign
barplot(pca$rotation[, 1], horiz = TRUE, main = "PC1")


heatmap(var(Fold_Change))



### cleaned up coloring 
metad.cl <- subset(Metadata, sample %in% intersect(Metadata$sample, pca.data$sample)) ## adjust row length of metadata to pca.data


metad.cl$mechanism <- Drug_Annotation$Mechanism[match(metad.cl$drug, Drug_Annotation$Drug)]
metad.cl$msi <- Cellline_Annotation$Microsatellite_instability_status[match(metad.cl$cell, Cellline_Annotation$Cell_Line_Name)]


plot(pca$x[,1], pca$x[,2], col = metad.cl$drug, xlab = "PC1", ylab = "PC2", main =  "Colored by drug")
plot(pca$x[,1], pca$x[,2], col = metad.cl$dose, xlab = "PC1", ylab = "PC2", main = "Colored by dose")
plot(pca$x[,1], pca$x[,2], col = metad.cl$tissue, xlab = "PC1", ylab = "PC2", main = "Colored by tissue")

plot(pca$x[,1], pca$x[,2], col = metad.cl$msi, xlab = "PC1", ylab = "PC2", main = "Colored by Microsatellite Instability Status")
plot(pca$x[,1], pca$x[,2], col = metad.cl$mechanism, xlab = "PC1", ylab = "PC2", main = "Colored by Mechanism of action")


## plotting all important PCs, try for loop
scores <- data.frame(colnames(Fold_Change), pca$x[,1:6])
pc1.2 <- qplot(x=PC1, y=PC2, data = scores, colour =metad.cl$drug) +
  theme(legend.position="none")
pc1.3 <- qplot(x=PC1, y=PC3, data = scores, colour =metad.cl$drug) +
  theme(legend.position="none")
pc2.3 <- qplot(x=PC2, y=PC3, data = scores, colour =metad.cl$drug) +
  theme(legend.position="none")

pc1.2 
pc1.3 
pc2.3
grid.arrange(pc1.2, pc1.3, pc2.3)




j<-1
while(j < 7){
  i <- 2
  while(i < 7){
    sapply(1:length(i), function(x){
      sapply(1:length(i), function(y){
        plot(x=pca$x[,j], y=pca$x[,i],  col = metad.cl$tissue, xlab = j, ylab = i, main = "Colored by tissue")
        })
      })
    i <- i+1}
  j <- j+1}
