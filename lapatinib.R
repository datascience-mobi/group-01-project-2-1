Metadata_Lapatinib_treated= Metadata[421:474, ]
Metadata_Lapatinib_untreated= Metadata[1240:1293, ]
summary(NegLogGI50[9, ])
mycol <- c("blue", "red", "green", "yellow", "violet", "lightblue", "dark violet","dark green","orange")
pie(table(Metadata_Lapatinib_treated [, "tissue"]), main = "Tissue types Lapatinib", radius = 1, cex = 1.5, col = mycol)