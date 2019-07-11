#PCA: May 23rd 2019
#applying fold change to continous variables 
library(ggplot2)
library(gridExtra)
library(viridis)


pca <- prcomp(t(Fold_Change), scale = TRUE)
#transform df by function t() because prcomp expects rows to be samples, which is not the case for Fold_Change
#goal: show relation between related samples

plot(pca$rotation[,1], pca$rotation[,2]) 

# get name of sample with highest pc1 value
rownames(pca$x)[order(pca$x[,1], decreasing = TRUE)[1]]



# sdev calculates variation each PC accounts for
pca.var <- pca$sdev^2  
# since percentages make more sense then normal variation values
# calculate % or variation, which is much more interesing
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 


barplot(pca.var.per, main = "Scree plot", xlab = "Principal Components", ylab = "% variation")
plot(pca.var.per[1:10], type = "l", xlab = "Principal Components", ylab = "% variation")
plot(pca.var.per[1:8], type = "l", xlab = "Principal Components", ylab = "% variation")
plot(cumsum(pca.var.per[1:15]), type = "l", xlab = "Principal Components", ylab = "% variation") 

#creating data frame with all pcs 
#cleaning up sample names as they differed between matrices 
pca.data <- data.frame(pca$x)
rownames(pca.data) <- gsub(x = rownames(pca.data), pattern = "X786", replacement = "786")
pca.data <- cbind(sample =rownames(pca.data), pca.data)

## get names of top 10 genes that contribute most to pc1
loading_scores_1 <- pca$rotation[,1]

loading_scores <- pca$rotation[,1:6]
ranked_pca <- sort(loading_scores[,1], decreasing = TRUE)
View(ranked_pca)

gene_score <- abs(loading_scores_1) ## sort magnitude
gene_score_ranked <- sort(gene_score, decreasing = TRUE)


top_10_genes <- names(gene_score_ranked[1:10])
top_10_genes # show names of top 10 genes



pca$rotation[top_10_genes,1] ##show scores and +/- sign



### Metadata color matrix for coloring 
Metadata$sample <- gsub(x = Metadata$sample, pattern = "-", replacement = ".")

metad.cl <- subset(Metadata, Metadata$sample %in% pca.data$sample) ## adjust row length of metadata to pca.data


metad.cl$mechanism <- Drug_Annotation$Mechanism[match(metad.cl$drug, Drug_Annotation$Drug)]
metad.cl$msi <- Cellline_Annotation$Microsatellite_instability_status[match(metad.cl$cell, Cellline_Annotation$Cell_Line_Name)]
View(metad.cl)



#overview of clustering by different PCs
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


#color vectors for coloring by drug and tissue
viridis <- viridis(9)
color_tissue = viridis[metad.cl$tissue]
tissue <- levels(metad.cl$tissue)

magma <- magma(15)
color_drug = magma[metad.cl$drug]
drug <- levels(metad.cl$drug)


## colored by drug
#plot PC1 and PC2
plot(pca$x[,1], 
     pca$x[,2], 
     col = color_drug,
     pch = 19, 
     xlab = paste("PC1 (",pca.var.per[1],"%)"), 
     ylab = paste("PC2 (",pca.var.per[2],"%)"))
#create legend
legend("topleft", 
       legend = drug, 
       col = magma, 
       pch = 19, 
       xpd = "TRUE",
       bty = "n"
)
#create title
mtext("PCA of Fold Change  colored by drug", 
      side = 3, 
      line = -2.5,
      cex = 1.2,
      font = 2, 
      outer = TRUE)



#plot PC2 and PC3
plot(pca$x[,2], 
     pca$x[,3], 
     col = color_drug,
     pch = 19, 
     xlab = paste("PC2 (",pca.var.per[2],"%)"), 
     ylab = paste("PC3 (",pca.var.per[3],"%)"))
#create legend
legend("right", 
       legend = drug, 
       col = magma, 
       pch = 19, 
       xpd = "TRUE",
       bty = "n",
       inset = c(0, 2)
)
#create title
mtext("PCA of Fold Change  colored by drug", 
      side = 3, 
      line = -2.5,
      cex = 1.2,
      font = 2, 
      outer = TRUE)


## colored by tissue
#plot PC1 and PC2
plot(pca$x[,1], 
     pca$x[,2], 
     col = color_tissue,
     pch = 19, 
     xlab = paste("PC1 (",pca.var.per[1],"%)"), 
     ylab = paste("PC2 (",pca.var.per[2],"%)"))
#create legend
legend("topleft", 
       legend = tissue, 
       col = viridis, 
       pch = 19, 
       xpd = "TRUE",
       bty = "n"
)
#create title
mtext("PCA of Fold Change  colored by tissue", 
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
legend("right", 
       legend = tissue, 
       col = viridis, 
       pch = 19, 
       xpd = "TRUE",
       bty = "n",
       inset = c(0, 2)
)
#create title
mtext("PCA of Fold Change  colored by tissue", 
      side = 3, 
      line = -2.5,
      cex = 1.2,
      font = 2, 
      outer = TRUE)




