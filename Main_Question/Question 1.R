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

combined1 = cbind(adjustedNeglogI50Lap, Fold_Changemeans)

names1 = c( "NegLogI50Lap","Fold_Changemeans")
colnames(combined1) = names1
                      
lmFold = lm(NegLogI50Lap ~ Fold_Changemeans, data = combined1)

summary(lmFold)

qqnorm(lmFold$residuals, main = "Test for normaldistribution of residuals")
qqline(lmFold$residuals)

plot(combined1$NegLogI50Lap, lmFold$fitted.values, pch = 20, col = "blue", xlab = "Real values", 
     ylab = "Predicted values", main = "Comparison: real and predicted values ~ linear regression (Fold_Changemeans)")
abline(0, 1, col = "red")

cor(combined1$NegLogI50Lap,combined1$Fold_Changemeans)




#Same wth doubling-time

NegLogGI50Lap = NegLogGI50[9,]

#Sort by Cellline-Name
df = arrange(Cellline_Annotation, Cell_Line_Name)
Doublingtime = cbind.data.frame (df$Cell_Line_Name, df$Doubling_Time)

c21 = as.data.frame(t(NegLogGI50Lap))

combined2 = cbind(c21, Doublingtime$`df$Doubling_Time`)
names2 = c( "NegLogI50Lap","Doubling_Time")
colnames(combined2) = names2

combined2 =na.omit(combined2)

lmDouble = lm(NegLogI50Lap ~ Doubling_Time, data = combined2)

summary(lmDouble)

qqnorm(lmDouble$residuals, main = "Test for normaldistribution of residuals")
qqline(lmDouble$residuals)

plot(combined2$NegLogI50Lap, lmDouble$fitted.values, pch = 20, col = "blue", xlab = "Real values", 
     ylab = "Predicted values", main = "Comparison: real and predicted values ~ linear regression (Doubling-Time)")
abline(0, 1, col = "red")

cor(combined2$NegLogI50Lap,combined2$Doubling_Time)






#Multiple regression

qwer = gsub(x =Doublingtime$`df$Cell_Line_Name`, pattern = "-", replacement = ".")
Doublingtime1 =  rbind(qwer,Doublingtime$`df$Doubling_Time`)
Doublingtime1 = as.data.frame(t(Doublingtime1)) 

c31 = subset(Doublingtime1, qwer %in% intersect(Doublingtime1$qwer, c2$V1))
c41 = as.numeric(as.character(c31$V2))
adjustedDoubling_Time = as.data.frame(c41)

combined3 = cbind(adjustedNeglogI50Lap, Fold_Changemeans, adjustedDoubling_Time)
names3 = c( "NegLogI50Lap","Fold_Changemeans","Doubling_Time")
colnames(combined3) = names3

mlr = lm(NegLogI50Lap ~ ., data = combined3)

summary(mlr)

qqnorm(mlr$residuals, main = "Test for normaldistribution of residuals")
qqline(mlr$residuals)

plot(combined3$NegLogI50Lap, mlr$fitted.values, pch = 20, col = "blue", xlab = "Real values", 
     ylab = "Predicted values" , main = "Comparison: real and predicted values ~ multiple regression")
abline(0, 1, col = "red")
