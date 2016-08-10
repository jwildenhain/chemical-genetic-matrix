# ML_AlternativeModels.R
#
# Copyright 2015 Jan Wildenhain
#
# This file is part of the Chemical Genetic Matrix 
# Scientific Data submission.
#
# This software is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published 
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.

# This software is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License 
# along with Foobar. If not, see http://www.gnu.org/licenses/.
#
# This script reads data for machine learning, loaded files are extracted 
# from ChemGRID database (http://chemgrid.org)
#
# Last modified: 31.07.2016
#
# Script performs:
# * read in data files (NB strain model, cryptic matrix )
# * run RF model with 10-fold cross validation
# * run RF regression model using Bliss Independence values with 10-fold cross validation
# * compare model with J48 tree to support MeanDecreaseGini and sentinel behaviour 
# * compare prediction with the naive Bayes classifier in R

kYourUserFolder <- "~/Documents/GitHub/chemical-genetic-matrix/Prediction"
setwd(kYourUserFolder)

# install required packages
source("../PackageInstallation.R") 

# load parameters for RF training
source("../LoadConfigProperties.R")
# working directory

# to shortcut the Random Forest learner you can load the cgm_RFModel_*** file here
# d2.rf <- read.csv(file="../Data/cgm_RFModel_***.csv.gz",head=TRUE,sep=",")
# 
d2.rf <- read.csv(file="../Data/RFModel_train_on_ext_data.csv.gz", head=TRUE, sep=",")
d2.rf <- d2.rf[ ,-1]
names(d2.rf)

balance.none = which(d2.rf$Synergy == "none")
row.none     = sample(balance.none, round(length(balance.none) * use.per.data))
length(row.none)
balance.syn  = which(d2.rf$Synergy == "synergist")
row.syn      = sample(balance.syn, round(length(balance.syn) * use.per.data))
length(row.syn)
balance.ant  = which(d2.rf$Synergy == "antagonist")
row.ant      = sample(balance.ant, round(length(balance.ant) * use.per.data))
length(row.ant)

data.train <- d2.rf[c(row.ant, row.none, row.syn), ]
dim(data.train)
data.test  <- d2.rf[c(row.ant, row.none, row.syn) * -1, ]
dim(data.test)
names(data.train)

dim(d2.rf)
names(d2.rf[,c(-2, -103, -204)])
fit.rf <- randomForest(factor(synergy) ~.,
                       data       = d2.rf[ ,c(-2, -103, -204)], 
                       sampsize   = c(110, 1400, 700),
                       importance = T,  
                       ntree      = rf.t,
                       nodesize   = rf.ns, 
                       mtry       = rf.mtry,
                       CV         = T)  #
print(fit.rf)  # view results 
importance(fit.rf)  # estimate importance of each sentinetal strain
varImpPlot(fit.rf, 
           type  = 2, 
           sort  = TRUE,
           n.var = 85,
           main  = "Sentinel Strains Importance")
dev.print(pdf, 
          width  = 10,
          height = 20,
          file   = paste0("Plots/NB-RF_v", rf.version,"_MeanDecreaseGini",
                          "_t", rf.t,
                          "_ns", rf.ns,
                          "_mtry", rf.mtry, ".pdf"))

tnum     <- unclass(d2.rf$synergy)
head(tnum)
proba.rf <- predict(fit.rf, type = "prob")[ ,2]
sc.rf    <- Logit(predict(fit.rf, type = "prob")[ ,2])
min(sc.rf)
max(sc.rf)
# Filter set Inf values to max min
rm.min = which(sc.rf == -Inf)
rp.min = min(sc.rf[rm.min * -1])
sc.rf[rm.min] <- rp.min

rm.max = which(sc.rf == Inf)
rp.max = min(sc.rf[rm.max * -1])
sc.rf[rm.max] <- rp.max
min(sc.rf)
max(sc.rf)

overlapDensity(sc.rf[tnum == 2] * -1,
               sc.rf[tnum == 3] * -1,
               gpnames = c("Mutual", "Synergistic"))
dev.print(pdf, 
          width  = 10, 
          height =  8,
          file   = paste0("Plots/NB-RF_v", rf.version,
                          "_OverlapDensity",
                          "_t", rf.t,
                          "_ns", rf.ns,
                          "_mtry", rf.mtry, ".pdf"))

ppp  <- sc.rf * -1
lll  <- d2.rf$bliss > bliss.syn.cutoff  #
pred <- prediction(ppp, lll)
perf <- performance(pred, "tpr", "fpr")
auc  <- performance(pred, "auc")@y.values[[1]]
auc
plot(perf, avg= "threshold", colorize = T, lwd = 3, main= "Random Forest Model ROC")
legend("bottomright", 
       legend = c(paste("Random Forest (AUC = ", formatC(auc, digits = 4, format="f"), ")", sep="")),
       col    = c("red"),
       lty    = 1)
dev.print(pdf, 
          width  = 10,
          height =  7,
          file   = paste0("Plots/NB-RF_v", rf.version, "_TPRvsFPR",
                          "_t", rf.t,
                          "_ns", rf.ns,
                          "_mtry", rf.mtry, ".pdf"))

perf <- performance(pred, measure = "acc", x.measure = "cutoff")
perf
# Now let's get the cutoff for the best accuracy
bestAccInd <- which.max(perf@"y.values"[[1]])
bestMsg    <- print(paste0("Best accuracy=", perf@"y.values"[[1]][bestAccInd],
                       " at cutoff=", round(perf@"x.values"[[1]][bestAccInd], 4)))
plot(perf, sub = bestMsg)
dev.print(pdf, 
          width  = 7,
          height = 7,
          file   = paste0("Plots/NB-RF_v", rf.version, "_ROC",
                          "_t", rf.t,
                          "_ns", rf.ns,
                          "_mtry", rf.mtry, ".pdf"))

pred <- prediction(ppp, lll)
perf <- performance(pred, "prec", "rec")
plot(perf, colorize=T)
dev.print(pdf, 
          width  = 10, 
          height = 7, 
          file   = paste0("Plots/NB-RF_v",rf.version,"_PrecisionRecall",
                          "_t", rf.t,
                          "_ns", rf.ns,
                          "_mtry", rf.mtry, ".pdf"))

plot(sc.rf[tnum >= 2] * -1, d2.rf$bliss[tnum >= 2], 
     pch  = 16,
     col  = alpha("darkgrey",0.5),
     xlab = "Synergy Score",
     ylab = "Bliss Independence")
dev.print(pdf, 
          width  = 10,
          height = 8, 
          file   = paste0("Plots/NB-RF_v",rf.version,"_Scatter", 
                          "_t", rf.t,
                          "_ns", rf.ns,
                          "_mtry", rf.mtry, ".pdf"))

cor(sc.rf[tnum >= 2]* -1, d2.rf$bliss[tnum >= 2])
summary(lm((sc.rf[tnum >= 2]* -1) ~ d2.rf$bliss[tnum >= 2]))

#
# RF regression classification
dim(d2.rf)
names(d2.rf[,c(-2, -103)])
length(which(d2.rf$synergy == "none" | d2.rf$synergy == "synergist"))
dim(d2.rf[which(d2.rf$synergy == "none" | d2.rf$synergy == "synergist"),c(-1, -2, -103)])

fit.rf.regression <- randomForest(bliss ~.,
                                  data       = d2.rf[which(d2.rf$synergy == "none" | d2.rf$synergy == "synergist"), c(-1,-2,-103)],
                                  sampsize   = c(10000),
                                  importance = T,
                                  ntree      = rf.t,
                                  nodesize   = rf.ns,
                                  mtry       = rf.mtry,
                                  CV         = T)  #
print(fit.rf.regression)
rf.regcor <- cor(d2.rf$bliss[which(d2.rf$synergy == "none" | d2.rf$synergy == "synergist")],
                 fit.rf.regression$predicted)
reg = lm(fit.rf.regression$predicted ~ d2.rf$bliss[which(d2.rf$synergy == "none" | d2.rf$synergy == "synergist")]) 
rf.rsq    <- summary(reg)$r.squared
plot(d2.rf$bliss[which(d2.rf$synergy == "none" | d2.rf$synergy == "synergist")], fit.rf.regression$predicted,
     pch  = 19, 
     cex  = 0.75,
     col  = alpha("grey", 0.5),
     xlab = "Bliss Independence",
     ylab = "RF Regression")
legend("bottomright",
       legend = c(paste0("P=", formatC(rf.regcor, digits = 4, format = "f"),
                        " R2=", formatC(rf.rsq, digits = 4, format = "f"), ")")),
       col    = c("red"),
       lty    = 0)
dev.print(pdf, 
          width  = 8,
          height = 8,
          file   = paste0("Plots/NB-RFregr_v", rf.version, "_Scatter",
                         "_t", rf.t,
                         "_ns", rf.ns, 
                         "_mtry", rf.mtry, ".pdf"))

varImpPlot(fit.rf.regression,
           type  = 2,
           sort  = TRUE,
           n.var = 85,
           main  = "RF Regression Model Node Purity")
dev.print(pdf, 
          width  = 10,
          height = 30,
          file   = paste0("Plots/NB-RFregr_v", rf.version,"_NodePurity",
                         "_t", rf.t,
                         "_ns", rf.ns,
                         "_mtry", rf.mtry, ".pdf"))

ppp  <- fit.rf.regression$predicted
lll  <- d2.rf$bliss[which(d2.rf$synergy == "none" | d2.rf$synergy == "synergist")] > bliss.syn.cutoff
pred <- prediction(ppp, lll)
perf <- performance(pred, "tpr", "fpr")
auc  <- performance(pred, "auc") 
auc
plot(perf, 
     avg      = "threshold",
     colorize = T,
     lwd      = 3,
     main     = "RF Regression Model ROC")
dev.print(pdf, 
          width  = 7,
          height = 7,
          file   = paste0("Plots/NB-RFregr_v",rf.version,"_ROC",
                         "_t", rf.t, 
                         "_ns", rf.ns,
                         "_mtry", rf.mtry, ".pdf"))

perf <- performance(pred, "prec", "rec")
plot(perf, colorize=T)
dev.print(pdf, 
          width  = 7,
          height = 7,
          file   = paste0("Plots/NB-RFregr_v", rf.version, "_PrecisionRecall",
                         "_t", rf.t,
                         "_ns", rf.ns,
                         "_mtry", rf.mtry, ".pdf"))

# Build J48 tree and compare with MeanDecreaseGini how sentinels trend towards importance
# in predicting synergy, additivity or antagonism

# Balance dataset to gain split information about sentinels
balance.none = which(d2.rf$synergy == "none")
row.none     = sample(balance.none, 1400)
row.none
balance.syn  = which(d2.rf$synergy == "synergist")
row.syn      = sample(balance.syn, 700)
row.syn

names(d2.rf)
fit <- ctree(factor(synergy) ~., 
             data = d2.rf[c(row.none, row.syn), c(-2, -103, -204)]) 

plot(fit, main="")
dev.print(pdf, 
          width  = 100,
          height = 34, 
          file   = paste0("Plots/J48_singletree_v", rf.version, ".pdf"))
summary(fit)

# Test model using the native Naive Bayes algorithm
dim(d2.rf)

names(d2.rf)

dim(d2.rf[c(row.none,row.syn),c(-1, -2, -103, -204)])

synergy <- d2.rf$synergy[c(row.none, row.syn)]

fit.nBmodel <- naiveBayes(synergy ~. ,
                          data = d2.rf[c(row.none, row.syn), c(-1, -2, -103, -204)])
fit.nBmodel
fit.nB <- predict(fit.nBmodel, d2.rf[c(-1, -2, -103, -204)])
fit.nB

tab <- table(fit.nB, d2.rf$synergy)
tab

# show ROC curve:
ppp  <- d2.rf$bliss # [c(row.none, row.syn)]
lll  <- as.vector(fit.nB) 
pred <- prediction(ppp, lll)
perf <- performance(pred, "tpr","fpr")
auc  <- performance(pred,"auc") 
auc

plot(perf, 
     avg      = "threshold",
     colorize = T,
     lwd      = 3,
     main     = "NB-NB Model ROC performance")
dev.print(pdf,
          width  = 7,
          height = 7,
          file   = paste0("Plots/NB-NB_performance_v", rf.version, ".pdf"))

