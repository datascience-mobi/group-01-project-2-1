# T.test
#Daten
#list = list(Treated,Untreated)
#nlist = lapply(list,scale)
#Treated = as.data.frame(nlist[[1]])

#Lapatinib col = 421:474, erlotinib col= 249:307
lap.treated = grep("lapatinib", colnames(Treated), value = FALSE)
lapatinib.treated = Treated[, lap.treated]

erl.treated = grep("erlotinib", colnames(Treated), value = FALSE)
erlotinib.treated = Treated[, erl.treated]

erl = data.frame(position = 1:819, celllines = colnames(Treated))[249:307,]
lap = data.frame(position = 1:819, celllines = colnames(Treated))[421:474,]
