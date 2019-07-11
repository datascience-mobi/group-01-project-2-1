# PCA specific for Lapatinib: June 6th 2019
library(ggplot2)
library(gridExtra)
library(dplyr)
library(RColorBrewer)


L_fc <- select(Fold_Change, contains("Lapa"))


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

pca.data <- data.frame(pca$x)
rownames(pca.data) <- gsub(x = rownames(pca.data), pattern = "X786", replacement = "786")
pca.data <- cbind(sample =rownames(pca.data), pca.data)


View(pca.data)

## get names of top 10 genes that contribute most to pc1
loading_scores_1 <- pca$rotation[,1]

loading_scores <- pca$rotation[,1:6]
ranked_pca <- sort(loading_scores[,1], decreasing = TRUE)
View(ranked_pca)

gene_score <- abs(loading_scores_1) ## sort magnitude
gene_score_ranked <- sort(gene_score, decreasing = TRUE)


top_10_genes <- names(gene_score_ranked[1:10])
top_10_genes # show names of top 10 genes


### Metadata matrix for coloring 
Metadata$sample <- gsub(x = Metadata$sample, pattern = "-", replacement = ".")
metad.cl <- subset(Metadata, sample %in% intersect(Metadata$sample, pca.data$sample))  
metad.cl$msi <- Cellline_Annotation$Microsatellite_instability_status[match(metad.cl$cell, Cellline_Annotation$Cell_Line_Name)]
metad.cl$inoculation_d <- Cellline_Annotation$Inoculation_Density[match(metad.cl$cell, Cellline_Annotation$Cell_Line_Name)]
metad.cl$doubling_time <- Cellline_Annotation$Doubling_Time[match(metad.cl$cell, Cellline_Annotation$Cell_Line_Name)]
metad.cl$cancer_type <- Cellline_Annotation$Cancer_type[match(metad.cl$cell, Cellline_Annotation$Cell_Line_Name)]




##overview of most important PCs
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




#color vectors for coloring by msi and tissue
colormsi <- brewer.pal(3, "Set1")
color_msi = colormsi[metad.cl$msi]
msi <- levels(metad.cl$msi)

magma <- magma(9)
color_tissue = magma[metad.cl$tissue]
tissue <- levels(metad.cl$tissue)


## colored by msi
#plot PC1 and PC2
plot(pca$x[,1], 
     pca$x[,2], 
     col = color_msi,
     pch = 19, 
     xlab = paste("PC1 (",pca.var.per[1],"%)"), 
     ylab = paste("PC2 (",pca.var.per[2],"%)"))
#create legend
legend("bottomright", 
       legend = msi, 
       col = colormsi, 
       pch = 19, 
       xpd = "TRUE",
       bty = "n"
)
#create title
mtext("PCA of Fold Change colored by MSI", 
      side = 3, 
      line = -2.5,
      cex = 1.2,
      font = 2, 
      outer = TRUE)

#plot PC2 and PC3
plot(pca$x[,2], 
     pca$x[,3], 
     col = color_msi,
     pch = 19, 
     xlab = paste("PC2 (",pca.var.per[2],"%)"), 
     ylab = paste("PC3 (",pca.var.per[3],"%)"))
#create legend
legend("topleft", 
       legend = msi, 
       col = colormsi, 
       pch = 19, 
       xpd = "TRUE",
       bty = "n"
)
#create title
mtext("PCA of Fold Change colored by MSI", 
      side = 3, 
      line = -2.5,
      cex = 1.2,
      font = 2, 
      outer = TRUE)



##colored by tissue
#plot PC1 and PC2
plot(pca$x[,1], 
     pca$x[,2], 
     col = color_tissue,
     pch = 19, 
     xlab = paste("PC1 (",pca.var.per[1],"%)"), 
     ylab = paste("PC2 (",pca.var.per[2],"%)"))
#create legend
legend("bottomright", 
       legend = tissue, 
       col = magma, 
       pch = 19, 
       xpd = "TRUE",
       bty = "n"
)
#create title
mtext("PCA of Fold Change colored by tissue", 
      side = 3, 
      line = -2.5,
      cex = 1.2,
      font = 2, 
      outer = TRUE)

#plot PC2 and PC3
plot(pca$x[,2], 
     pca$x[,3], 
     col = color_tissue,
     pch = 19, 
     xlab = paste("PC2 (",pca.var.per[2],"%)"), 
     ylab = paste("PC3 (",pca.var.per[3],"%)"))
#create legend
legend("bottomleft", 
       legend = tissue, 
       col = magma, 
       pch = 19, 
       xpd = "TRUE",
       bty = "n"
)
#create title
mtext("PCA of Fold Change colored by tissue", 
      side = 3, 
      line = -2.5,
      cex = 1.2,
      font = 2, 
      outer = TRUE)

