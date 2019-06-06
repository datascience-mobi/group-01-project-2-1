
#biomarker matrix for breast cancer cell lines 
breast_celllines <- c("MCF7", "MDA-MB-231", "MDA-MB-468", "HS-578T", "BT-549","T-47D" )
biomarker_genes <-c("ESR1", "ERBB2", "ERBB3", "ERBB4", "KRT75", "PGR")
Matrix_biomarker_breast<-CCLE_basalexpression[biomarker_genes, breast_celllines]


min(Matrix_biomarker_breast)
max(Matrix_biomarker_breast)

#Variance analysis
Variance_Breast <- apply(Matrix_biomarker_breast, 1, var)
summary(Variance_Breast)

boxplot(Matrix_biomarker_breast)

data_biomarker_breast<-as.matrix(Matrix_biomarker_breast)
head(data_biomarker_breast)
heatmap(data_biomarker_breast)

#pretty heatsmap
install.packages("pheatmap")
library(pheatmap)
pheat_breast<-pheatmap(Matrix_biomarker_breast, show_rownames = TRUE, show_colnames = FALSE, cutree_cols = 4, cutree_rows = 3, drop_levels = TRUE, clustering_method = "ward.D2", scale="row")

#biomarker matrix for all cell lines 
biomarker_genes <-c("ESR1", "ERBB2" ,"ERBB3", "ERBB4", "KRT75", "PGR")
Matrix_biomarker<-CCLE_basalexpression[biomarker_genes, ]
summary(Matrix_biomarker)

min(Matrix_biomarker)
max(Matrix_biomarker)

boxplot(Matrix_biomarker)


data_biomarker<-as.matrix(Matrix_biomarker)
head(data_biomarker)
heatmap(data_biomarker)

#pretty heatsmap
pheat_all_genes<-pheatmap(Matrix_biomarker, show_rownames = TRUE, show_colnames = FALSE, cutree_cols = 4, cutree_rows = 3, drop_levels = TRUE, clustering_method = "ward.D2", scale="row")


#Ergänzung
install.packages("ggpubr")
library("ggpubr")

T_Matrix_biomarker_breast<-t(Matrix_biomarker_breast)


#setting a threshold, matrix of the highest gene expression
rmv.rows = apply(CCLE_basalexpression, 1, function(x) {
 sum(x<8)})
which(rmv.rows <8)
highest_expression = CCLE_basalexpression[-which(rmv.rows <8), ]  
rm(rmv.rows)
summary(highest_expression)
pheat_highest_expression<-pheatmap(highest_expression, show_rownames = TRUE, show_colnames = TRUE, cutree_cols = 4, cutree_rows = 3, drop_levels = TRUE, clustering_method = "ward.D2", scale="row")

#genome wide paired t-test
before=as.matrix(LapatimibUnTreat)
after=as.matrix(LapatimibTreat)
t.test(before, after, paired=TRUE)


install.packages("matrixTests")
library(matrixTests)

#t-test over each column -> cell lines
col_t_paired(LapatimibUnTreat, LapatimibTreat, alternative = "two.sided", mu = 0,conf.level = 0.95)

#t-test over each row -> genes
row_t_test<-row_t_paired(LapatimibUnTreat, LapatimibTreat, alternative = "two.sided", mu = 0,conf.level = 0.95)
summary(row_t_test$pvalue)


