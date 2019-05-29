# Load data
wd = dirname(rstudioapi::getSourceEditorContext()$path)
Cellline_Annotation = read.delim(paste0(wd, "/Data/cellline_annotation.tsv"), header = TRUE, sep = "\t", stringsAsFactors = TRUE)

# Pearson correlation+plot

cor(Cellline_Annotation$Doubling_Time, Cellline_Annotation$Inoculation_Density, method = "pearson")
# 0,3209821

plot(Cellline_Annotation$Doubling_Time, Cellline_Annotation$Inoculation_Density, pch= 16, col= "blue", main = "Pearson correlation between Doubling Time and Inoculation Density", xlab = "Doubling Time", ylab = "Inoculation Density")

lm(Cellline_Annotation$Inoculation_Density~ Cellline_Annotation$Doubling_Time)

#   Call:
#   lm(formula = Cellline_Annotation$Inoculation_Density ~ Cellline_Annotation$Doubling_Time)
#
#   Coefficients:
#   (Intercept)  Cellline_Annotation$Doubling_Time  
#       9110.9                              169.4   

abline(9110.9, 169.4, col= "red", lwd =2)

# better:

abline(lm(Cellline_Annotation$Inoculation_Density ~ Cellline_Annotation$Doubling_Time), col = "red", lwd = 2)