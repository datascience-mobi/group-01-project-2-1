# PCA specific for Lapatinib: June 6th 2019
library(ggplot2)
library(gridExtra)

Fold_Change = NCI_TPW_gep_treated - NCI_TPW_gep_untreated
Treated = data.frame(NCI_TPW_gep_treated)
Untreated = data.frame(NCI_TPW_gep_untreated)


L_fc <- Fold_Change[,421:474]

# PCA
pca <- prcomp(t(L_fc), scale = TRUE)
plot(pca$x[,1], pca$x[,2]) #overview

rownames(pca$x)[order(pca$x[,1], decreasing = TRUE)[1]] # smaple with highest pc1 value

pca.var <- pca$sdev^2  # sdev calculates variation each PC accounts for
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 


barplot(pca.var.per, main = "Scree plot", xlab = "Principal Components", ylab = "% variation")
plot(pca.var.per[1:10], type = "l", xlab = "Principal Components", ylab = "% variation")
plot(pca.var.per[1:8], type = "l", xlab = "Principal Components", ylab = "% variation") # either 2 or 4 PCs
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
tail(gene_score_ranked)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes # show names of top 10 genes


metad.cl <- subset(Metadata, sample %in% intersect(Metadata$sample, pca.data$sample))  

Cell_Anno <- subset(Cellline_Annotation, Cell_Line_Name %in% intersect(metad.cl$cell, Cellline_Annotation$Cell_Line_Name))
Cell_Anno <- Cell_Anno[ order(Cell_Anno$Cell_Line_Name), ] # sort alphabetically
Cell_Anno <- cbind(Cell_Anno, sample = metad.cl$sample)

plot(pca$x[,1], pca$x[,2], col = Cell_Anno$Doubling_Time, xlab = "PC1", ylab = "PC2", main = "Colored by tissue")


## plotting all important PCs, try for loop
scores <- data.frame(colnames(L_fc), pca$x[,1:6])
pc1.2 <- qplot(x=PC1, y=PC2, data = scores, colour =metad.cl$tissue) +
  theme(legend.position="none")
pc1.3 <- qplot(x=PC1, y=PC3, data = scores, colour =metad.cl$tissue) +
  theme(legend.position="none")
pc2.3 <- qplot(x=PC2, y=PC3, data = scores, colour =metad.cl$tissue) +
  theme(legend.position="none")

pc1.2 
pc1.3 
pc2.3
grid.arrange(pc1.2, pc1.3, pc2.3)


pcs <- sapply(1:4, function(i){
  pca$x[,i]
}
)


plot(x=pcs, y=pcs,  col = metad.cl$drug, xlab = "PC"+i, ylab = "PC"+i+1, main = "Colored by drug")

i <- 1
while(i < 6){
  plot(x=pca$x[,i], y=pca$x[,i+1],  col = metad.cl$drug, xlab = i, ylab = i+1, main = "Colored by drug")
  i <- i+1}
