
#biomarker matrix for breast cancer cell lines 
breast_celllines <- c("MCF7", "MDA-MB-231", "MDA-MB-468", "HS-578T", "BT-549","T-47D" )
biomarker_genes <-c("ESR1", "ERBB2", "ERBB3", "ERBB4", "KRT75", "PGR")
Matrix_biomarker_breast=CCLE_basalexpression[biomarker_genes, breast_celllines]
summary(Matrix_biomarker_breast)

min(Matrix_biomarker_breast)
max(Matrix_biomarker_breast)

data_biomarker=as.matrix(Matrix_biomarker_breast)
head(data_biomarker_breast)
heatmap(data_biomarker_breast)


#biomarker matrix for all cell lines 
biomarker_genes <-c("ESR1", "ERBB2" ,"ERBB3", "ERBB4", "KRT75", "PGR")
Matrix_biomarker=CCLE_basalexpression[biomarker_genes, ]
summary(Matrix_biomarker)

min(Matrix_biomarker)
max(Matrix_biomarker)

data_biomarker=as.matrix(Matrix_biomarker)
head(data_biomarker)
heatmap(data_biomarker)


#setting a threshold, matrix of the highest gene expression
rmv.rows = apply(CCLE_basalexpression, 1, function(x) {
 sum(x<8)})
which(rmv.rows <8)
highest_expression = CCLE_basalexpression[-which(rmv.rows <8), ]  
rm(rmv.rows)
summary(highest_expression)

#genome wide paired t-test
before=as.matrix(LapatimibUnTreat)
after=as.matrix(LapatimibTreat)
t.test(before, after, paired=TRUE)



#funktioniert noch nicht, ist aber vielleicht eine bessere Variante

library(genefilter)
install.packages("matrixTests")

col_t_paired(LapatimibUnTreat, LapatimibTreat, alternative = "two.sided", mu = 0,conf.level = 0.95)


res2 <- vector(nrow(NCI_TPW_gep_treated), mode="list")
 
for(i in 1:nrow(NCI_TPW_gep_treated)) {
    res2[[i]] <- t.test(LapatimibTreat[i,], LapatimibUnTreat[i,], paired = TRUE)
 }
View(res2)
res2[1:2]
