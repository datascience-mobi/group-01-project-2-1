library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

NCI_TPW_gep_treated = readRDS(paste0(wd, "/Rohdaten/NCI_TPW_gep_treated.rds"))
dat = NCI_TPW_gep_treated[, 421:474]
LapatimibTreat <- data.frame(dat)
rm(dat)
NCI_TPW_gep_untreated = readRDS(paste0(wd, "/Rohdaten/NCI_TPW_gep_untreated.rds"))
dat = NCI_TPW_gep_untreated[, 421:474]
LapatimibUnTreat <- data.frame(dat)
rm(dat)

Metadata = read.delim(paste0(wd, "/Rohdaten/NCI_TPW_metadata.tsv"))
Cellline_Annotation = read.delim(paste0(wd, "/Rohdaten/cellline_annotation.tsv"))
Drug_Annotation = read.delim(paste0(wd, "/Rohdaten/drug_annotation.tsv"))
CCLE_mutations = readRDS(paste0(wd, "/Rohdaten/CCLE_mutations.rds"))
CCLE_copynumber = readRDS(paste0(wd, "/Rohdaten/CCLE_copynumber.rds"))
CCLE_basalexpression = readRDS(paste0(wd, "/Rohdaten/CCLE_basalexpression.rds"))
NegLogGI50 = readRDS(paste0(wd, "/Rohdaten/NegLogGI50.rds"))

