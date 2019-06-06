
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
Fold_Change = NCI_TPW_gep_treated - NCI_TPW_gep_untreated
Treated = data.frame(NCI_TPW_gep_treated)
Untreated = data.frame(NCI_TPW_gep_untreated)
Fold_Change <- data.frame(Fold_Change)