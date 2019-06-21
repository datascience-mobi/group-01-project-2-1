library(dplyr)
library(data.table)
library(matrixTests)


install.packages("plotly")
suppressPackageStartupMessages(library("plotly"))

update.packages("ggplot2")
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

L_fc <- select(Fold_Change, contains("Lapa"))
L_fc <- as.data.frame(t(L_fc))
rownames(Metadata) <- Metadata$sample


L_treated <- select(Treated, contains("Lapa"))
L_treated <- t(L_treated)
L_untreated <- select(Untreated, contains("Lapa"))
L_untreated <- t(L_untreated)



# selecting breast Lapatinib samples
breast <- Metadata[Metadata[,'tissue']=="Breast",]
rownames(breast) <- breast$sample
rownames(breast) <- gsub(x = rownames(breast), pattern = "-", replacement = ".")  

breastFC <- subset(L_fc, rownames(L_fc) %in% rownames(breast))
breastTreated <- subset(L_treated, rownames(L_treated) %in% rownames(breast))
breastUntreated <- subset(L_untreated, rownames(L_untreated) %in% rownames(breast))#


# selecting CNS Lapatinib samples
cns <- Metadata[Metadata[,'tissue']=="CNS",]
rownames(cns) <- cns$sample
rownames(cns) <- gsub(x = rownames(cns), pattern = "-", replacement = ".")

cnsFC <- subset(L_fc, rownames(L_fc) %in% rownames(cns))
cnsTreated <- subset(L_treated, rownames(L_treated) %in% rownames(cns))
cnsUntreated <- subset(L_untreated, rownames(L_untreated) %in% rownames(cns))


t_test_cns <- col_t_paired(cnsTreated, cnsUntreated, alternative = "two.sided", mu = 0,conf.level = 0.95)
t_test_breast <- col_t_paired(breastTreated, breastUntreated, alternative = "two.sided", mu = 0,conf.level = 0.95)
pval_cns <- t_test_cns$pvalue
pval_breast <- t_test_breast$pvalue

fdr_cns <- p.adjust(pval_cns, "BH")
fdr_breast <- p.adjust(pval_breast, "BH")


breastFCm <- as.numeric(colMeans(breastFC))
cnsFCm <- as.numeric(colMeans(cnsFC))
genes <- colnames(breastFC)

## breast volcano plot

diff_df_breast <- data.frame(gene = genes, Fold = breastFCm, FDR = fdr_breast)
head(diff_df_breast)


# add a grouping column; default value is "not significant"
diff_df_breast["group"] <- "NotSignificant"

diff_df_breast["group"] <- if_else((-log10(diff_df_breast['FDR']) < 0.5 & abs(diff_df_breast['Fold'] > 0.4)),"group"<- "Significant&FoldChange", "group"<- "NotSignificant")



# change the grouping for the entries a large enough Fold change but not a low enough p value
#diff_df_breast[(abs(diff_df_breast['Fold']) > 1.5 ),"group"] <- "FoldChange"

# change the grouping for the entries with both significance and large enough fold change
#diff_df_breast[which(diff_df_breast['FDR'] < 0.05 & abs(diff_df_breast['Fold']) > 4 ),"group"] <- "Significant&FoldChange"



# Find and label the top peaks..
top_peaks_breast <- diff_df_breast[with(diff_df_breast, order(Fold, FDR)),][1:5,]
top_peaks_breast <- rbind(top_peaks_breast, diff_df_breast[with(diff_df_breast, order(-Fold, FDR)),][1:5,])


# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
a <- list()
for (i in seq_len(nrow(top_peaks_breast))) {
  m <- top_peaks_breast[i, ]
  a[[i]] <- list(
    x = m[["Fold"]],
    y = -log10(m[["FDR"]]),
    text = m[["gene"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}

plot_breat <- plot_ly(data = diff_df_breast, x = diff_df_breast$Fold, y = -log10(diff_df_breast$FDR), text = diff_df_breast$gene, mode = "markers", color = diff_df_breast$group) %>% 
  layout(title ="Volcano Plot of Lapatinib breast cancer samples") %>%
  layout(annotations = a)
plot_breast


## CNS volcano plot

diff_df_cns <- data.frame(gene = genes, Fold = cnsFCm, FDR = fdr_cns)
head(diff_df_cns)


# add a grouping column; default value is "not significant"
diff_df_cns["group"] <- "NotSignificant"

diff_df_cns["group"] <- if_else((-log10(diff_df_cns['FDR']) < 0.5 & abs(diff_df_cns['Fold'] > 0.4)),"group"<- "Significant&FoldChange", "group"<- "NotSignificant")


# Find and label the top peaks..
top_peaks_cns <- diff_df_cns[with(diff_df_cns, order(Fold, FDR)),][1:5,]
top_peaks_cns <- rbind(top_peaks_cns, diff_df_cns[with(diff_df_cns, order(-Fold, FDR)),][1:5,])


# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
a <- list()
for (i in seq_len(nrow(top_peaks_cns))) {
  m <- top_peaks_cns[i, ]
  a[[i]] <- list(
    x = m[["Fold"]],
    y = -log10(m[["FDR"]]),
    text = m[["gene"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}

plot_cns <- plot_ly(data = diff_df_cns, x = diff_df_cns$Fold, y = -log10(diff_df_cns$FDR), text = diff_df_cns$gene, mode = "markers", color = diff_df_cns$group) %>% 
  layout(title ="Volcano Plot of Lapatinib CNS cancer samples") %>%
  layout(annotations = a)
plot_cns



#EnhancedVolcano(diff_df_breast,
#                lab = diff_df_breast$gene,
#                x = 'Fold',
#                y = 'FDR',
#                xlim = c(-0.5, 1))


#install.packages("manhattanly")
#library(manhattanly)
#volcanoly(diff_df_breast, FC = c("Fold"), pval = c("FDR"))