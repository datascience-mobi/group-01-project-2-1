#PCA by Eva, May 23rd 2019
#capplying fold change to continous variables 
Fold_Change = NCI_TPW_gep_treated - NCI_TPW_gep_untreated
Treated = data.frame(NCI_TPW_gep_treated)
Untreated = data.frame(NCI_TPW_gep_untreated)


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
plot(pca.var.per[1:15], type = "l", xlab = "Principal Components", ylab = "% variation")
plot(cumsum(pca.var.per[1:150]), type = "l", xlab = "Principal Components", ylab = "% variation") 



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



  
azac <- ifelse(cb[,1] == "TRUE", "green", "black")
bort <- ifelse(cb[,2] == "TRUE", "violet", "black")
cisp <- ifelse(cb[,3] == "TRUE", "blue", "black")
dasa <- ifelse(cb[,4] == "TRUE", "magenta", "black")
doxo <- ifelse(cb[,5] == "TRUE", "yellow", "black")
erlo <- ifelse(cb[,6] == "TRUE", "red", "black")
geld <- ifelse(cb[,7] == "TRUE", "orange", "black")
gemc <- ifelse(cb[,8] == "TRUE", "lightblue", "black") 
lapa <- ifelse(cb[,9] == "TRUE", "aquamarine", "black")
pacl <- ifelse(cb[,10] == "TRUE", "chocolate", "black")
siro <- ifelse(cb[,11] == "TRUE", "coral", "black")
sora <- ifelse(cb[,12] == "TRUE", "azure", "black")
suni <- ifelse(cb[,13] == "TRUE", "coral", "black")
topo <- ifelse(cb[,14] == "TRUE", "lightseagreen", "black")
vori <- ifelse(cb[,15] == "TRUE", "plum", "black")

cb_col <- cbind(azac,bort,cisp,dasa,doxo,erlo,geld,gemc,lapa,pacl,siro,sora,suni,topo,vori)



plot(pca$rotation[,1], pca$rotation[,2], col = cb_col, pch = 19, xlab= "PC1", ylab = "PC2")

# Problem: black not ideal as else color -> is there a way to give no color at all (evtl. for loop?)



## get names of top 10 genes that contribute most to pc1
loading_scores_1 <- pca$rotation[,1]

loading_scores <- pca$rotation[,1:6]
ranked_pca1 <- sort(loading_scores[,1], decreasing = TRUE)
View(ranked_pca1)

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




