library(CCA)
library(candisc) 

rm(list = ls())
bioData <- read.csv("http://msekce.karlin.mff.cuni.cz/~maciak/NMST539/bioData.csv", header = T)
chemData <- read.csv("http://msekce.karlin.mff.cuni.cz/~maciak/NMST539/chemData.csv", header = T)

head(bioData)
head(chemData)

ind <- match(chemData[,1], bioData[,1])

data <- data.frame(bioData[ind, ], chemData[, 2:8])

write.csv(data, file = "my_data_frame.csv", row.names = FALSE)

View(data)
X <- data[,2:9]
Y <- data[,19:25]

#Correlation between two sets
matcor(X,Y)

cc_model <- cc(X,Y)
cc_model$cor

par(mfrow = c(1,1))

barplot(cc_model$cor, main = "Canonical correlations for 'cancor()'", col = "gray")

library(CCP)

rho <- cc_model$cor
n <- dim(X)[1]
p <- dim(X)[2]
q <- dim(Y)[2]

#Wilk's test
p.asym(rho,n,p,q,tstat = "Wilks")

#Lawley-Hotelling test
p.asym(rho,n,p,q,tstat = "Hotelling")

#Pillai's test
p.asym(rho,n,p,q,tstat = "Pillai")

#Roy's largest root test
p.asym(rho,n,p,q,tstat = "Roy")

#Test the independence between two sets of variables
Wilks(cancor(X,Y))


#Significant canonical correlations
cc_model$cor[1:3]

#Squared canoncal correlations
cc_model$cor[1:3]^2

cc_model$xcoef

cc_model$ycoef

loadings <- comput(X,Y,cc_model)

loadings$corr.X.xscores

loadings$corr.Y.yscores

loadings$corr.X.yscores

loadings$corr.Y.xscores
