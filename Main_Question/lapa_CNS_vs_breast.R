library(dplyr)
library(data.table)
library(matrixTests)
library(plotly)
library(ggrepel)
library(compare)
library(viridis)
library(pheatmap)


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
diff_df_breast$absFold <- abs(diff_df_breast$Fold)
head(diff_df_breast)


# add a grouping column; default value is "not significant"
diff_df_breast$group <- "NotSignificant"



# change the grouping for the entries with significance but not a large enough Fold change
diff_df_breast[which(diff_df_breast['FDR'] < 0.5 & (diff_df_breast['absFold']) < 0.2 ),"group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
diff_df_breast[which(diff_df_breast['FDR'] > 0.5 & (diff_df_breast['absFold']) > 0.2 ),"group"] <- "FoldChange"

# change the grouping for the entries with both significance and large enough fold change
diff_df_breast[which(diff_df_breast['FDR'] < 0.5 & (diff_df_breast['absFold']) > 0.2 ),"group"] <- "Significant&FoldChange"


View(diff_df_breast)



# Find and label the top peaks..
top_peaks_breast <- diff_df_breast[with(diff_df_breast, order(Fold, FDR)),][1:10,]
top_peaks_breast <- rbind(top_peaks_breast, diff_df_breast[with(diff_df_breast, order(-Fold, FDR)),][1:10,])


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

plot_breast <- plot_ly(data = diff_df_breast, x = diff_df_breast$Fold, y = -log10(diff_df_breast$FDR), type = "scatter", text = diff_df_breast$gene, mode = "markers", color = diff_df_breast$group) %>% 
  layout(title ="Volcano Plot of Lapatinib breast cancer samples", 
         xaxis = list(title="log2 Fold Change"),
         yaxis = list(title="FDR")) %>%
  layout(annotations = a)
plot_breast





## CNS volcano plot

diff_df_cns <- data.frame(gene = genes, Fold = cnsFCm, FDR = fdr_cns)
diff_df_cns$absFold <- abs(diff_df_cns$Fold)
head(diff_df_cns)


# add a grouping column; default value is "not significant"
diff_df_cns$group <- "NotSignificant"


# change the grouping for the entries with significance but not a large enough Fold change
diff_df_cns[which(diff_df_cns['FDR'] < 0.5 & (diff_df_cns['absFold']) < 0.2 ),"group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
diff_df_cns[which(diff_df_cns['FDR'] > 0.5 & (diff_df_cns['absFold']) > 0.2 ),"group"] <- "FoldChange"

# change the grouping for the entries with both significance and large enough fold change
diff_df_cns[which(diff_df_cns['FDR'] < 0.5 & (diff_df_cns['absFold']) > 0.2 ),"group"] <- "Significant&FoldChange"


View(diff_df_cns)



# Find and label the top peaks..
top_peaks_cns <- diff_df_cns[with(diff_df_cns, order(Fold, FDR)),][1:10,]
top_peaks_cns <- rbind(top_peaks_cns, diff_df_cns[with(diff_df_cns, order(-Fold, FDR)),][1:10,])


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

plot_cns <- plot_ly(data = diff_df_cns, x = diff_df_cns$Fold, y = -log10(diff_df_cns$FDR),type = "scatter", text = diff_df_cns$gene, mode = "markers", color = diff_df_cns$group) %>% 
  layout(title ="Volcano Plot of Lapatinib CNS cancer samples",
         xaxis = list(title="log2 Fold Change"),
         yaxis = list(title="FDR"))%>%
  layout(annotations = a)
plot_cns



# selecet top peak genes common in cns and breast tissue
tpb_comparison <- subset(top_peaks_breast, gene %in% top_peaks_cns$gene)
View(tpb_comparison)
tpc_comparison <- subset(top_peaks_cns, gene %in% top_peaks_breast$gene)
View(tpc_comparison)


# order common genes alphabetically
tpb_comparison <- tpb_comparison[order(tpb_comparison$gene),]
tpc_comparison <- tpc_comparison[order(tpc_comparison$gene),]




## creating heat map of FCs to compare values 
cor_mat <- cbind("breast" = tpb_comparison$Fold, "cns" = tpc_comparison$Fold)
rownames(cor_mat) <- tpb_comparison$gene



pheatmap(
  mat               = cor_mat,
  color             = magma(10),
  border_color      = "black",
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Comparison:
  FC levels of CNS and breast top peak genes"
)


#correlation between top peak gene expression patterns in breast and cns tissues treated by lapatinib
cor_mat <- as.data.frame(cor_mat)
dif.FC.BC = data.frame(breast = cor_mat$breast -  mean(t(breastFC)), cns = cor_mat$cns -  mean(t(cnsFC))) #breast data - mean value, cns data - mean value
dif.FC.BC$gene = rownames(cor_mat)



# PLot

ggplot(dif.FC.BC,aes(x = gene, y = breast)) +
  geom_bar(stat = "identity", fill = "skyblue") + 
  geom_text(aes(label = round(breast, 2)), vjust = -0.5, color = "black", size = 3) + 
  coord_flip() + labs(title = "Mean graph plot of Fold Change for breast top peak genes")


ggplot(dif.FC.BC,aes(x = gene, y = cns)) +
  geom_bar(stat = "identity", fill = "skyblue") + 
  geom_text(aes(label = round(cns, 2)), vjust = -0.5, color = "black", size = 3) + 
  coord_flip() + labs(title = "Mean graph plot of Fold Change for CNS top peak genes")

#correlation
cor(cor_mat$breast, cor_mat$cns, method = "pearson")
