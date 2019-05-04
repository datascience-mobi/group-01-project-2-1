NCI_TPW_gep_treated <- readRDS("C:/Users/tobia/Desktop/R-Codes/Daten/NCI_TPW_gep_treated.rds")
dat = NCI_TPW_gep_treated[, 421:474]
LapatimibTreat <- data.frame(dat)
rm(dat)
NCI_TPW_gep_treated <- readRDS("C:/Users/tobia/Desktop/R-Codes/Daten/NCI_TPW_gep_untreated.rds")
dat = NCI_TPW_gep_untreated[, 421:474]
LapatimibUnTreat <- data.frame(dat)
rm(dat)
Metadata <- read.delim("C:/Users/tobia/Desktop/R-Codes/Daten/NCI_TPW_metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
Cellline_Annotation <- read.delim("C:/Users/tobia/Desktop/R-Codes/Daten/cellline_annotation.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
Drug_Annotation <- read.delim("C:/Users/tobia/Desktop/R-Codes/Daten/drug_annotation.tsv", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
CCLE_mutations <- readRDS("C:/Users/tobia/Desktop/R-Codes/Daten/CCLE_mutations.rds")
CCLE_copynumber <- readRDS("C:/Users/tobia/Desktop/R-Codes/Daten/CCLE_copynumber.rds")
NegLogGI50 <- readRDS("C:/Users/tobia/Desktop/R-Codes/Daten/NegLogGI50.rds")