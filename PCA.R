#PCA by Eva, May 23rd 2019
#capplying fold change to continous variables 
Fold_Change = NCI_TPW_gep_treated - NCI_TPW_gep_untreated
Treated = data.frame(NCI_TPW_gep_treated)
Untreated = data.frame(NCI_TPW_gep_untreated)


pca <- prcomp(t(Fold_Change), scale = TRUE)
#transform df by function t() because prcomp expects rows to be samples, which is not the case for Fold_Change
#goal: show relation between related samples

plot(pca$x[,1], pca$x[,2]) 

pca.var <- pca$sdev^2  ##sdev calculates variation each PC accounts for
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
#since percentages make more sense then normal variation values
#calculate % or variation, which is much more interesing

barplot(pca.var.per, main = "Scree plot", xlab = "Principal Components", ylab = "% variation")
plot(pca.var.per, type = "l", xlab = "Principal Components", ylab = "% variation")
plot(pca.var.per[1:15], type = "l", xlab = "Principal Components", ylab = "% variation")
plot(cumsum(pca.var.per[1:150]), type = "l", xlab = "Principal Components", ylab = "% variation")

library(ggplot2) ##for fancy pca plots with a lot of information
#format data for ggplot

pca.data <- data.frame(Sample = rownames(pca$x),
                       x = pca$x[,1],
                       y = pca$x[,2])
pca.data
View(pca.data)

rownames(pca.data)

###This is a test for color coding by drug, not yet finsihed 
###pls ignore until code runs perfectly

#for all drugs -> 15 differet colors for plotting pca later
## Azacytidine samples
azacytidine <- pca.data[1:39,1]
cb1 <- rownames(pca.data) %in% azacytidine

#ifelse(cb1$1, "TRUE", "forestgreen", "red")
#View(cb1)

#Bortezomib samples
bortezomib <- pca.data[40:96,1]
cb2 <- rownames(pca.data) %in% bortezomib

#Cisplatin samples
cisplatin <- pca.data[97:149,1]
cb3 <- rownames(pca.data) %in% cisplatin

#Dasatinib samples
dasatinib <- pca.data[150:198,1]
cb4 <- rownames(pca.data) %in% dasatinib

#Doxorubicin samples
doxorubicin <- pca.data[199:248,1]
cb5 <- rownames(pca.data) %in% doxorubicin

#Erlotinib samples 
erlotinib <- pca.data[249:307,1]
cb6 <- rownames(pca.data) %in% erlotinib

#Geldanamycin samples
geldanamycin <- pca.data[308:364,1]
cb7 <- rownames(pca.data) %in% geldanamycin

#Gemcitibine samples
gemcitibine <- pca.data[365:420,1]
cb8 <- rownames(pca.data) %in% gemcitibine

#Lapatinib samples
lapatinib <- pca.data[421:474,1]
cb9 <- rownames(pca.data) %in% lapatinib

#Paclitaxel samples
paclitaxel <- pca.data[475:533,1]
cb10 <- rownames(pca.data) %in% paclitaxel

#Sirolimus samples
sirolimus <- pca.data[534:589,1]
cb11 <- rownames(pca.data) %in% sirolimus

#Sorafenib samples
sorafenib <- pca.data[590:646,1]
cb12 <- rownames(pca.data) %in% sorafenib

#Sunitinib samples
sunitinib <- pca.data[647:702,1]
cb13 <- rownames(pca.data) %in% sunitinib

#Topotecan sample
topotecan <- pca.data[703:760,1]
cb14 <- rownames(pca.data) %in% topotecan

#Vorinostat samples
vorinostat <- pca.data[761:819,1]
cb15 <- rownames(pca.data) %in% vorinostat

cb <- cbind('azacytidine' = cb1, 'bortezomib' = cb2, 'cisplatin' = cb3, 'dasatinib' = cb4,
            'doxorubicin' = cb5, 'erlotinib' = cb6, 'geldanamycin' = cb7, 'gemcitibine' = cb8,
            'lapatinib' = cb9, 'paclitaxel' = cb10, 'sirolimus' = cb11, 'sorafenib' = cb12,
            'sunitinib' = cb13, 'topotecan' = cb14, 'vorinostat' = cb15)

View(cb)
new.cb <- as.data.frame(cb)


  
a <- ifelse(cb[,1] == "TRUE", "green", "red")

View(a)




plot(pca$rotation[,2], pca$rotation[,3], col = a, pch = 19, xlab= "PC1", ylab = "PC2")
###End of test


###PCA with ggplot not working, idk why
###either correction or discard

ggplot(pca.data, aes(x=pca$x[,1], y=pca$x[,2])) #+
 # geom_text() +
#  xlab(paste("PC1", pca.var.per[,1], "%", sep = "")) +
#  ylab(paste("PC2", pca.var.per[,2], "%", sep = "")) +
#  theme_bw() +
#  ggtitle("My PCA")

loading_scores <- pca$rotation[,1:6]
ranked_pca1 <- sort(loading_scores[,1])
View(ranked_pca1)

gene_score <- abs(loading_scores) ## sort magnitude
mean(gene_score)
max(gene_score)
min(gene_score)

gene_score_ranked <- sort(gene_score, decreasing = TRUE)

head(gene_score_ranked)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes



pca$rotation[top_10_genes,1] ##show scores and +/- sign
barplot(pca$rotation[, 1], horiz = TRUE, main = "PC1")


heatmap(var(Fold_Change))




