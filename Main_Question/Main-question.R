  install.packages("dplyr")
  library(dplyr)
  

lapa<-data.frame(Metadata[which(Metadata[,'drug'] == "lapatinib"), ])
erlo<-data.frame(Metadata[which(Metadata[,'drug'] == "erlotinib"), ])

el<-right_join(lapa,erlo, by="cell")
el

rmv.rows = apply(el, 1, function(x) {
  sum(is.na(x))
})  # Go through each row and sum up all missing values
row.names(rmv.rows)



fc<-(Treated-Untreated)
fc<-data.frame(scale(fc))
all<-data.frame(fc[grep("lapatinib|erlotinib", colnames(fc))])
all

#since rlotinip contains more columns than lapatinib, we have to remove these columns
library(dplyr)
all<-data.frame(all %>% select(-"CCRF.CEM_erlotinib_0nM_24h", 
                    -"HL.60_erlotinib_0nM_24h", 
                    -"HT29_erlotinib_0nM_24h", 
                    -"K.562_erlotinib_0nM_24h", 
                    -"LOX_erlotinib_0nM_24h",
                    -"SR_erlotinib_0nM_24h",
                    -"COLO.205_lapatinib_0nM_24h"))
                   

la<-data.frame(all[grep("lapatinib", colnames(all))])
ncol(la)

er<-data.frame(all[grep("erlotinib", colnames(all))])
ncol(er)

erla<-data.frame(er,la)

ncol(all) #to prove if the columns are removed


#ggplot(lapatinib, erlotinib)

drug<-c(rep('Erlotinib',53), rep('Lapatinib',53))
expression_drug<-apply(erla, MARGIN = 2, sum)
df_drug<-data.frame(expression_drug, drug)

boxplot(erlotinib, lapatinib)
library(ggpubr)
ggboxplot (data = df_drug, x="drug", y="expression_drug", color = "drug",    
          add = "jitter", legend = "none")+
          rotate_x_text(angle = 45)+
          geom_hline(yintercept = mean(lapatinib$MCF7_lapatinib_0nM_24h), linetype = 2)+ # Add horizontal line at base mean
          stat_compare_means(method = "anova")+        # Add global annova p-value
          stat_compare_means(label = "p.signif", method = "t.test",
          ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all



