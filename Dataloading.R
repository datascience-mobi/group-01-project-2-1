
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

NCI_TPW_gep_treated = readRDS(paste0(wd, "/Data/NCI_TPW_gep_treated.rds"))
NCI_TPW_gep_untreated = readRDS(paste0(wd, "/Data/NCI_TPW_gep_untreated.rds"))
Metadata = read.delim(paste0(wd, "/Data/NCI_TPW_metadata.tsv"), header = TRUE, sep = "\t", stringsAsFactors = TRUE)
Cellline_Annotation = read.delim(paste0(wd, "/Data/cellline_annotation.tsv"), header = TRUE, sep = "\t", stringsAsFactors = TRUE)
Drug_Annotation = read.delim(paste0(wd, "/Data/drug_annotation.tsv"), header = TRUE, sep = "\t", stringsAsFactors = TRUE)
CCLE_mutations = readRDS(paste0(wd, "/Data/CCLE_mutations.rds"))
CCLE_copynumber = readRDS(paste0(wd, "/Data/CCLE_copynumber.rds"))
CCLE_basalexpression = readRDS(paste0(wd, "/Data/CCLE_basalexpression.rds"))
NegLogGI50 = as.data.frame(readRDS(paste0(wd, "/Data/NegLogGI50.rds")))
Treated = data.frame(NCI_TPW_gep_treated)
Untreated = data.frame(NCI_TPW_gep_untreated)


#For scaling the data (taken from and analysed  in "Broad")

list = list(Treated,Untreated)
nlist = lapply(list,scale)
Treated = as.data.frame(nlist[[1]])
Untreated = as.data.frame(nlist[[2]])
Fold_Change = Treated - Untreated
Fold_Change <- data.frame(Fold_Change)
rm(NCI_TPW_gep_treated,NCI_TPW_gep_untreated,list,nlist)
