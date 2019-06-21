library(dplyr)

Fold_ChangeLap = select(Fold_Change, contains("Lapa"))
NegLogGI50Lap = NegLogGI50[9,]
means = colMeans(Fold_ChangeLap)
Fold_Changemeans = as.data.frame(t(means))

asdf = gsub(x = colnames (Fold_Changemeans), pattern = "_lapatinib_10000nM_24h", replacement = "")
colnames(Fold_Changemeans) = asdf

asdfa = gsub(x = asdf, pattern = "X7", replacement = "7")
colnames(Fold_Changemeans) = asdfa


asd = gsub(x = colnames (NegLogGI50Lap), pattern = "-", replacement = ".")
colnames(NegLogGI50Lap) = asd

c1 = rbind(asd,NegLogGI50Lap)
c2 = rbind(asdfa,Fold_Changemeans)

c1 = t(c1)
c2 = t(c2)

c1 =as.data.frame(c1)
c2 =as.data.frame(c2)



c3 = subset(c1, `1` %in% intersect(c1$`1`, c2$V1))
c4 = as.numeric(as.character(c3$lapatinib))
adjustedNeglogI50Lap = as.data.frame(c4)
#adjustedNeglogI50Lap = as.data.frame(t(adjustedNeglogI50Lap))


Fold_Changemeans = as.data.frame(t(Fold_Changemeans))

combined = cbind(adjustedNeglogI50Lap, Fold_Changemeans)

names = c( "NegLogI50Lap","Fold_Changemeans")
colnames(combined) = names
                      
lmFold = lm(NegLogI50Lap ~ Fold_Changemeans, data = combined)

qqnorm(lmFold$residuals)
qqline(lmFold$residuals)

plot(combined$NegLogI50Lap, lmFold$fitted.values, pch = 20, col = "blue", xlab = "Real values", 
     ylab = "Predicted values")
abline(0, 1, col = "red")
