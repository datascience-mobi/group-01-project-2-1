
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
#breast_celllines <- c("MCF7", "MDA-MB-231", "MDA-MB-468", "HS-578T", "BT-549","T-47D" )

breastcells<-subset(Metadata, tissue == "Breast")
cellnames<-data.frame(cells=rownames(CCLE_basalexpression))
breast<-subset(CCLE_basalexpression, breastcells$cell %in%  cells)

tbreastcells<-data.frame(t(breast))
biomarker_genes <-c("ESR1", "ERBB2", "ERBB3", "ERBB4", "KRT75", "PGR")

breast_marker <- c(tbreastcells$ESR1,tbreastcells$ERBB2, tbreastcells$ERBB3, tbreastcells$ERBB4, tbreastcells$KRT75, tbreastcells$PGR )
breast_marker_matrix <- data.frame(tbreastcells$ESR1,tbreastcells$ERBB2, tbreastcells$ERBB3, tbreastcells$ERBB4, tbreastcells$KRT75, tbreastcells$PGR )

#Variance analysis
Variance_Breast <- apply(breast, 1, var)
summary(Variance_Breast)

boxplot(breast)

#ANOVA
cells<-c(rep('ESR1',45),rep('ERBB2',45),rep('ERBB3',45),rep('ERBB4',45),rep('KRT75',45),rep('PGR',45))
expression<-c(tbreastcells$ESR1,tbreastcells$ERBB2, tbreastcells$ERBB3, tbreastcells$ERBB4, tbreastcells$KRT75, tbreastcells$PGR )
df_breast<-data.frame(cells,expression)
plot(expression ~ cells, data=df_breast)
breast.aov<-aov(expression ~ cells, data=df_breast)
summary(breast.aov)



#ggboxplot of cell lines
library(ggpubr)
ggboxplot(df_breast, x = "cells", y = "expression", color="cells")



#ggboxplot of biomarkers for breast
library(ggpubr)
ggboxplot(df_breast, x = "cells", y = "expression", color="cells",
          add = "jitter", legend = "none")+
          rotate_x_text(angle = 45)+
          geom_hline(yintercept = mean(Biomarker_df$ERBB2), linetype = 2)+ # Add horizontal line at base mean
          stat_compare_means(method = "kruskal.test", label.y = 12)+        # Add global p-value
          stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all





#pretty heatsmap
#install.packages("pheatmap")
#install.packages("viridis")

library(pheatmap)
library(viridis)

pheat_breast<-pheatmap(breast_marker_matrix ,
                       show_rownames = FALSE, 
                       show_colnames = TRUE, 
                       cutree_cols = 4, 
                       cutree_rows = 3, 
                       color = viridis(4),
                       drop_levels = TRUE, 
                       clustering_method = "ward.D2", 
                       scale="row")
pheat_breast

#pheatmap of all genes and breast cells
all_genes_cells<-pheatmap(breast, 
         show_rownames = FALSE, 
         show_colnames = TRUE, 
         color = viridis(4),
         drop_levels = TRUE, 
         clustering_method = "ward.D2", 
         scale="row")
all_genes_cells

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
t_basalexpression<-data.frame(t(CCLE_basalexpression))
biomarker_genes <-c("ESR1", "ERBB2" ,"ERBB3", "ERBB4", "KRT75", "PGR")
Matrix_biomarker<-c(t_basalexpression$ESR1, t_basalexpression$ERBB2, t_basalexpression$ERBB3, t_basalexpression$ERBB4, t_basalexpression$KRT75, t_basalexpression$PGR)
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
  stat_compare_means(method = "kruskal.test", label.y = 12)+        # Add global annova p-value
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



library(pheatmap)
library(viridis)
pheatmap(highest_expression, 
                            show_rownames = FALSE, 
                            show_colnames = FALSE, 
                            color=viridis(4),
                            drop_levels = TRUE, 
                            clustering_method = "ward.D2", 
                            scale="row")

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

library(ggpubr)
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

row_t_test$pvalue


pvalues<-sort(row_t_test$pvalue,decreasing=FALSE)
head(pvalues)

basalexpression <- data.frame(CCLE_basalexpression)
b<-apply(basalexpression,1, max)
highest_basalexpression<-sort(b, decreasing = TRUE)
head(highest_basalexpression)

row_t_test

Treated_t<-data.frame(t(Treated))
Untreated_t<-data.frame(t(Untreated))

