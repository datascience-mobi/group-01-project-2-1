#importing data
NCI_TPW_gep_treated <- readRDS("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/NCI_TPW_gep_treated.rds")
dat = NCI_TPW_gep_treated[, 421:474]
LapatimibTreat <- data.frame(dat)
rm(dat)
NCI_TPW_gep_untreated <- readRDS("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/NCI_TPW_gep_untreated.rds")
dat = NCI_TPW_gep_untreated[, 421:474]
LapatimibUnTreat <- data.frame(dat)
rm(dat)
Metadata <- read.delim("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/NCI_TPW_metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
Cellline_Annotation <- read.delim("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/cellline_annotation.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
Drug_Annotation <- read.delim("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/drug_annotation.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
CCLE_mutations <- readRDS("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/CCLE_mutations.rds")
CCLE_copynumber <- readRDS("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/CCLE_copynumber.rds")
NegLogGI50 <- readRDS("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/NegLogGI50.rds")
CCLE_basalexpression <- readRDS("C:/Users/schentarra/Documents/Uni/4. Fachsemester/Bioinfo/Rohdaten/CCLE_basalexpression.rds")


  
library(ggplot2)
ggplot(CCLE_mutations, aes(x=Chromosome, fill=Variant_Classification)) + geom_bar()
#ggplot(Cellline_Annotation, aes(x=Cancer_type, fill=Microsatellite_instability_status)) + geom_bar()
#ggplot(Cellline_Annotation, aes(x=Cell_Line_Name, y=Doubling_Time, z=Inoculation_Density)) + geom_contour()
#Computation failed in `stat_contour()`:
#no proper 'z' matrix specified 




#c_annotation <- log2(Biobase::exprs(Cellline_Annotation))
#PCA_raw <- prcomp(t(Cellline_Annotation), scale. = FALSE)

#percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
#sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
 #                    Cancer = pData(raw_data)$Cancer_type,
  #                   MIS = pData(raw_data)$Microsatellite_instability_status,
   #                  Cellline = pData(raw_data)$Cell_Line_Name.)

#ggplot(dataGG, aes(PC1, PC2)) +
#  geom_point(aes(shape = Disease, colour = Phenotype)) +
#  ggtitle("PCA plot of the log-transformed raw expression data") +
#  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
#  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
#  theme(plot.title = element_text(hjust = 0.5))+
#  coord_fixed(ratio = sd_ratio) +
#  scale_shape_manual(values = c(4,15)) + 
#  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

colnames(CCLE_basalexpression)
dim(CCLE_basalexpression)
colnames(CCLE_copynumber)
dim(CCLE_copynumber)
colnames(CCLE_mutations)
dim(CCLE_mutations)
#hugo symbol is a gene nomenclatur
colnames(Cellline_Annotation)
Cellline_Annotation[,1]
dim(Cellline_Annotation)
#microsatellite instability status is a predisposition to mutation
colnames(Drug_Annotation)
Drug_Annotation[,1]
dim(Drug_Annotation)
colnames(Metadata)
dim(Metadata)
colnames(LapatimibTreat)
colnames(LapatimibUnTreat)

