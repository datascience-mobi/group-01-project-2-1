library(dplyr)

Fold_ChangeLap = select(Fold_Change, contains("Lapa"))
NegLogGI50Lap = NegLogGI50[9,]
means = colMeans(Fold_ChangeLap)
Fold_Changemeans = as.data.frame(t(means))

a2 = gsub(x = colnames (Fold_Changemeans), pattern = "_lapatinib_10000nM_24h", replacement = "")
colnames(Fold_Changemeans) = a2

a3 = gsub(x = a2, pattern = "X7", replacement = "7")
colnames(Fold_Changemeans) = a3


a1 = gsub(x = colnames (NegLogGI50Lap), pattern = "-", replacement = ".")
colnames(NegLogGI50Lap) = a1

c1 = rbind(a1,NegLogGI50Lap)
c2 = rbind(a3,Fold_Changemeans)

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

#Split the data (Training - Testing)

n = nrow(combined1)
rmse1 = sqrt(1/n * sum(lmFold$residuals^2))
rmse1

i1.train = sample(1:nrow(combined1), 44)

dat1.train = combined1[i1.train, ]
dat1.test = combined1[-i1.train, ]

l1.train = lm(NegLogI50Lap ~ Fold_Changemeans, data = dat1.train)
summary(l1.train)

n = nrow(dat1.train)
rmse1.train = sqrt(1/n * sum(l1.train$residuals^2))
rmse1.train

pred1 = predict(l1.train, newdata = dat1.test)

n = nrow(dat1.test)
residuals = dat1.test$NegLogI50Lap - pred1
rmse1.test1 = sqrt(1/n * sum(residuals^2))
rmse1.test1






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

#Split the data (Training - Testing)

n = nrow(combined2)
rmse2 = sqrt(1/n * sum(lmDouble$residuals^2))
rmse2

i2.train = sample(1:nrow(combined2), 48)

dat2.train = combined2[i2.train, ]
dat2.test = combined2[-i2.train, ]

l2.train = lm(NegLogI50Lap ~ Doubling_Time, data = dat2.train)
summary(l2.train)

n = nrow(dat2.train)
rmse2.train = sqrt(1/n * sum(l2.train$residuals^2))
rmse2.train

pred2 = predict(l2.train, newdata = dat2.test)

n = nrow(dat1.test)
residuals = dat2.test$NegLogI50Lap - pred2
rmse2.test = sqrt(1/n * sum(residuals^2))
rmse2.test




#Multiple regression

b1 = gsub(x =Doublingtime$`df$Cell_Line_Name`, pattern = "-", replacement = ".")
Doublingtime1 =  rbind(b1,Doublingtime$`df$Doubling_Time`)
Doublingtime1 = as.data.frame(t(Doublingtime1)) 

c31 = subset(Doublingtime1, b1 %in% intersect(Doublingtime1$b1, c2$V1))
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



#Split the data (Training - Testing)

n = nrow(combined3)
rmse3 = sqrt(1/n * sum(mlr$residuals^2))
rmse3

i3.train = sample(1:nrow(combined2), 44)

dat3.train = combined3[i3.train, ]
dat3.test = combined3[-i3.train, ]

l3.train = lm(NegLogI50Lap ~ ., data = dat3.train)
summary(l3.train)

n = nrow(dat3.train)
rmse3.train = sqrt(1/n * sum(l3.train$residuals^2))
rmse3.train

pred3 = predict(l3.train, newdata = dat3.test)

n = nrow(dat3.test)
residuals = dat3.test$NegLogI50Lap - pred3
rmse3.test = sqrt(1/n * sum(residuals^2))
rmse3.test