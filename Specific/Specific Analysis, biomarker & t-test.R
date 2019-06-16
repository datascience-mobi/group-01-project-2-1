
Treated_t<-data.frame(t(Treated))
Untreated_t<-data.frame(t(Untreated))

library(viridis)
boxplot(Untreated_t$ESR1, Treated_t$ESR1, 
        Untreated_t$ERBB2, Treated_t$ERBB2,
        Untreated_t$ERBB3, Treated_t$ERBB3, 
        Untreated_t$ERBB4, Treated_t$ERBB4,
        Untreated_t$KRT75, Treated_t$KRT75, 
        Untreated_t$PGR, Treated_t$PGR,
        col = viridis(12),
        names = c("ESR1 before", "ESR1 after", "ERBB2 before", "ERBB2 after", "ERBB3 before", "ERBB3 after", "ERBB4 before","ERBB4 after", "KRT75 before", "KRT75 after", "PGR before",  "PGR after"),
        main="Expression of Biomakers before and after treatment with Lapatinib", 
        xlab="Biomarker",
        ylab="expression")


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

#ANOVA
cells<-c(rep('MCF7',6),rep('MDA-MB-231',6),rep('MDA-MB-468',6),rep('HS-578T',6),rep('BT-549',6),rep('T-47D',6))
expression<-c(Matrix_biomarker_breast$MCF7, Matrix_biomarker_breast$`MDA-MB-231`, Matrix_biomarker_breast$`MDA-MB-468`, Matrix_biomarker_breast$`HS-578T`, Matrix_biomarker_breast$`BT-549`, Matrix_biomarker_breast$`T-47D`)
df_breast<-data.frame(cells,expression)
plot(expression ~ cells, data=df_breast)
breast.aov<-aov(expression ~ cells, data=df_breast)
summary(breast.aov)

#ggboxplot of cell lines
library(ggpubr)
ggboxplot(df_breast, x = "cells", y = "expression", color="cells")


#ggboxplot of biomarkers for breast
t_breast<-t(Matrix_biomarker_breast)
t_breast_df<-data.frame(t_breast)
Marker_Breast<-c(rep('ESR1',6), rep('ERBB2',6), rep('ERBB3',6), rep('ERBB4',6), rep('KRT75',6), rep('PGR',6))
expression_breast<-c(t_breast_df$`ESR1`, t_breast_df$`ERBB2`, t_breast_df$`ERBB3`, t_breast_df$`ERBB4`, t_breast_df$`KRT75`, t_breast_df$`PGR`)
df_breast_marker<-data.frame(expression_breast, Marker_Breast)
plot(expression_breast~Marker_Breast, data=df_breast_marker)
ggboxplot(df_breast_marker, x="Marker_Breast", y="expression_breast", color = "Marker_Breast",
          add = "jitter", legend = "none")+
          rotate_x_text(angle = 45)+
          geom_hline(yintercept = mean(Biomarker_df$ERBB2), linetype = 2)+ # Add horizontal line at base mean
          stat_compare_means(method = "anova", label.y = 12)+        # Add global annova p-value
          stat_compare_means(label = "p.signif", method = "t.test",
          ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all


data_biomarker_breast<-as.matrix(Matrix_biomarker_breast)
head(data_biomarker_breast)
heatmap(data_biomarker_breast)

#pretty heatsmap
install.packages("pheatmap")
install.packages("viridis")

library(pheatmap)
library(viridis)
pheat_breast<-pheatmap(data_biomarker_breast, 
                       show_rownames = TRUE, 
                       show_colnames = FALSE, 
                       cutree_cols = 4, 
                       cutree_rows = 3, 
                       color = viridis(12),
                       drop_levels = TRUE, 
                       clustering_method = "ward.D2", 
                       scale="row")

#biomarker matrix for all cell lines 
biomarker_genes <-c("ESR1", "ERBB2" ,"ERBB3", "ERBB4", "KRT75", "PGR")
Matrix_biomarker<-CCLE_basalexpression[biomarker_genes, ]
summary(Matrix_biomarker)

min(Matrix_biomarker)
max(Matrix_biomarker)

boxplot(Matrix_biomarker)

#ggboxplot
library(ggpubr)
Biomarker_df<-data.frame(t(Matrix_biomarker))
Marker<-c(rep('ESR1',45), rep('ERBB2',45), rep('ERBB3',45), rep('ERBB4',45), rep('KRT75',45), rep('PGR',45))
expression<-c(Biomarker_df$`ESR1`, Biomarker$`ERBB2`, Biomarker$`ERBB3`, Biomarker$`ERBB4`, Biomarker$`KRT75`, Biomarker$`PGR`)
df_marker<-data.frame(expression, Marker)
plot(expression~Marker, data=df_marker)
ggboxplot(df_marker, x="Marker", y="expression", color = "Marker",
        add = "jitter", legend = "none")+
        rotate_x_text(angle = 45)+
        geom_hline(yintercept = mean(Biomarker_df$ERBB2), linetype = 2)+ # Add horizontal line at base mean
        stat_compare_means(method = "anova", label.y = 12)+        # Add global annova p-value
        stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all


data_biomarker<-as.matrix(Matrix_biomarker)
head(data_biomarker)
heatmap(data_biomarker)

#pretty heatsmap
pheat_all_genes<-pheatmap(Matrix_biomarker, 
                          show_rownames = TRUE, 
                          show_colnames = FALSE, 
                          cutree_cols = 4, 
                          cutree_rows = 3, 
                          color = viridis(4),
                          drop_levels = TRUE, 
                          clustering_method = "ward.D2", 
                          scale="row")


#Erg?nzung
install.packages("ggpubr")
library("ggpubr")
Matrix_biomarker_breast$MCF7<-as.factor(Matrix_biomarker_breast$MCF7)

T_Matrix_biomarker_breast<-t(Matrix_biomarker_breast)
ggplot(Matrix_biomarker_breast, aes(x=MCF7, y="MDA-MB-231"))  
     geom_boxplot(outlier.colour="red")

#setting a threshold, matrix of the highest gene expression
rmv.rows = apply(CCLE_basalexpression, 1, function(x) {
 sum(x<14)})
which(rmv.rows <14)
highest_expression = CCLE_basalexpression[-which(rmv.rows <14), ]  
rm(rmv.rows)
summary(highest_expression)
pheat_highest_expression<-pheatmap(highest_expression, 
                                   show_rownames = TRUE, 
                                   show_colnames = TRUE, 
                                   cutree_cols = 4, 
                                   cutree_rows = 3, 
                                   drop_levels = TRUE, 
                                   clustering_method = "ward.D2", 
                                   scale="row")

#searching for our own biomarkers
basalexpression <- data.frame(CCLE_basalexpression)
b<-apply(basalexpression,1, max)
highest_basalexpression<-sort(b, decreasing = TRUE)
head(highest_basalexpression)

#boxplot of our biomarkers before and after treatment
library(viridis)
boxplot(Untreated_t$COX6C, Treated_t$COX6C, 
        Untreated_t$'2', Treated_t$'2',
        Untreated_t$RPS11, Treated_t$RPS11, 
        Untreated_t$GAPDH, Treated_t$GAPDH,
        Untreated_t$MT2A, Treated_t$MT2A, 
        Untreated_t$TUBA1B, Treated_t$TUBA1B,
        col = viridis(12),
        names = c("COX6C before", "COX6C after", 
                  "2 before", "2 after", 
                  "RPS11 before", "RPS11 after", 
                  "GAPDH before","GAPDH after", 
                  "$MT2A before", "$MT2A after", 
                  "TUBA1B before",  "TUBA1B after"),
        main="Expression of Biomakers before and after treatment with Lapatinib", 
        xlab="Biomarker",
        ylab="expression")

#genome wide paired t-test
before=as.matrix(LapatimibUnTreat)
after=as.matrix(LapatimibTreat)
t.test(before, after, paired=TRUE)


#variance analysis before and after treatment
my_data <- data.frame( 
               group = rep(c("before", "after")),
                expression = c(before,  after)
               )
print(my_data)

ggboxplot(my_data, x = "group", y = "expression", 
                    color = "group", palette = c("#00AFBB", "#E7B800"),
                    order = c("before", "after"),
                    ylab = "expression", xlab = "group")



install.packages("matrixTests")
library(matrixTests)

#t-test over each column -> cell lines
col_t_paired(LapatimibUnTreat, LapatimibTreat, alternative = "two.sided", mu = 0,conf.level = 0.95)

#t-test over each row -> genes
row_t_test<-row_t_paired(LapatimibUnTreat, LapatimibTreat, alternative = "two.sided", mu = 0,conf.level = 0.99)
summary(row_t_test$pvalue)
row_t_test

Treated_t<-data.frame(t(Treated))
Untreated_t<-data.frame(t(Untreated))
