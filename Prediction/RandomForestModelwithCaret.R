# RandomForestModelwithCaret.R
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
# having that matrix written to a file you can run the mapping on the syngery
# matrix und then next build the RF model 
# 
# * read in data files (NB strain model, cryptic matrix )
# * run RF model for prediction using caret
# * compare model results with the randomForest package
# * compare prediction with the naive Bayes classifier in R

kYourUserFolder <- "~/Documents/GitHub/Chemical_Genetic_Matrix/Prediction"
setwd(kYourUserFolder)

# source required libraries
source("../PackageInstallation.R")

# load parameters for RF training
source("../LoadConfigProperties.R")
# working directory


# Filter for NB models that actually have an AUC
nb.auc <- read.csv("../Data/ECFP4_2013_nbtable.tab.csv.gz",
                   header      = TRUE, 
                   check.names = FALSE,
                   sep         = "\t") # choose an adjacency
head(nb.auc)
nb.auc.string = do.call(paste, 
                       c(as.list(nb.auc$Strain[nb.auc$AUC > 0.6]),
                       sep = "','"))
nb.auc.string = paste("'", nb.auc.string, "'", sep = "")
nb.auc.string

# load full NB model dataset 6.4 MB
d <- read.csv(file = "../Data/ChemGRID_sentinels_NBextract_2013.csv.gz",
              head = TRUE, 
              sep  = ",")
colnames(d)[1] <- "cid"
d$cid <- as.character(d$cid)
head(d)
dim(d)

# Data contains matrix of compound identifiers (rows) and strains (columns)
given <- read.csv(file = "../Data/cgm_cryptic_combinations.csv.gz",
                  head = TRUE,
                  sep  = ",")
given$X <- NULL
head(given)
dim(given)

# merge cids, bliss and sentinel data
df.nBtm1 <- merge(d, given, by.x = "cid", by.y = "cid1")
dim(df.nBtm1)
df.nBtm  <- merge(d, df.nBtm1, by.x = "cid", by.y = "cid2")
dim(df.nBtm)

#
# add the follow-up data to extend synergistic combos for learners
d              <- read.csv(file = "../Data/ChemGRID_sentinels_NBextract_ext.csv.gz",
                           head = TRUE,
                           sep  = ",")
colnames(d)[1] <- "cid"
d$cid <- as.character(d$cid)
head(d)
dim(d)

#
# cryptic matrix with containing 128 Spectrum compounds 
given   <- read.csv(file = "../Data/cgm_cryptic_combinations_ext.csv.gz",
                    head = TRUE,
                    sep  = ",")
given$X <- NULL
head(given)
dim(given)


#
# merge cids, bliss and sentinel data
df.nBetm1 <- merge(d, given, by.x = "cid", by.y = "cida")
dim(df.nBetm1)
df.nBetm <- merge(d, df.nBetm1, by.x = "cid",by.y = "cidb")
dim(df.nBetm)
head(df.nBetm)

df.nBetm[ ,c(1, 102, 203)]
dim(df.nBetm)

bp.data = boxplot(df.nBtm$bliss)
bp.data$stats

# check sentinel dimensions and order order
dim(df.nBtm)
dim(df.nBetm)

head(df.nBtm)
head(df.nBetm)
names(df.nBtm)
names(df.nBetm)

df.nBtm.ext = rbind(df.nBtm, df.nBetm)
dim(df.nBtm.ext)
df.nBtm     = df.nBtm.ext


# Apply statistics / machine learning algorithms (do PCA, RF and J48)
synergy <- rep("none", nrow(df.nBtm))
length(which(synergy == "none"))
synergy[df.nBtm$bliss > bliss.syn.cutoff] <- "synergist"
length(which(synergy == "synergist"))
synergy[df.nBtm$bliss < bliss.ant.cutoff] <- "antagonist"
length(which(synergy == "antagonist"))

d2 <- cbind(synergy, df.nBtm)
colnames(d2)

#
# Do random forest analysis using caret
d2.rev <- d2[ ,c(1, 103:203, 2:102, 204)]  # create reverse obj
colnames(d2.rev) <- colnames(d2)
head(d2)
d2.rf  <- rbind(d2, d2.rev)  # make symmetric matrix
syn    <- d2.rf$synergy  # make Syn vector
syn
names(d2.rf)
factor(syn)
head(d2.rf[, c(-1, -2, -103, -204)])  # remove obvious non informative columns for learner

dim(d2.rf)
write.csv(d2.rf, 
          file      = "../Data/RFModel_train_on_ext_data_vtmp.csv",
          quote     = TRUE,
          row.names = T)
system("GZIP = -9 gzip ../Data/RFModel_train_on_ext_data_vtmp.csv");
system("mv ../Data/RFModel_train_on_ext_data_vtmp.csv.gz ../Data/RFModel_train_on_ext_data.csv.gz");

# to shortcut the Random Forest learner you can load the cgm_RFModel_*** file here
# d2.rf <- read.csv(file="../Data/cgm_RFModel_***.csv.gz",head=TRUE,sep=",")
# d2.rf <- read.csv(file="../Data/RFModel_train_on_ext_data.csv.gz",head=TRUE,sep=",")

balance.none = which(d2.rf$synergy == "none")
row.none     = sample(balance.none, round(length(balance.none) * use.per.data))
length(row.none)
balance.syn  = which(d2.rf$synergy == "synergist")
row.syn      = sample(balance.syn, round(length(balance.syn) * use.per.data))
length(row.syn)
balance.ant  = which(d2.rf$synergy == "antagonist")
row.ant      = sample(balance.ant, round(length(balance.ant) * use.per.data))
length(row.ant)

data.train <- d2.rf[c(row.ant, row.none, row.syn), ]
data.test  <- d2.rf[c(row.ant, row.none, row.syn)*-1, ]
names(data.train)

# RF classification optimisation
ctrl <- trainControl(method          = "CV", 
                     classProbs      = TRUE, 
                     summaryFunction = twoClassSummary,
                     number          = 5,
                     repeats         = 10)

rf.estimate <- train(factor(synergy) ~., 
                     data=data.train[which(data.train$synergy == "none" | data.train$synergy == "synergist"), c(-2,-103,-204)],
                     method     = "rf",
                     metric     = "ROC", 
                     ntree      = rf.t,
                     nodesize   = rf.ns,
                     tuneGrid   = rf.grid, 
                     trControl  = ctrl, 
                     importance = T)

# produce figures for analysis on training and testset

# rfE500 <- rfEstimate
trellis.par.set(caretTheme())
plot(rf.estimate, ylim = c(0.5, 0.93))
dev.print(pdf, width  = 8,
               height = 8, 
               file   = paste0("Plots/NB-RFcaret_v", kParamVersion, "_Estimate",
                              "_p", use.per.data,
                              "_ns", rf.ns, 
                              "_t", rf.t, ".pdf" ))
# Figures for Train dataset
proba.rf <- predict(rf.estimate, newdata = data.train, type = "prob")
head(proba.rf)
dim(proba.rf)

tnum  <- unclass(data.train$synergy)
sc.rf <- Logit(proba.rf[ ,2])

plot(sc.rf, data.train$bliss, 
     pch  = 16,
     main = "Bliss Independence and Synergy Score", 
     ylab = "Bliss Independence", 
     xlab = "Synery Score", 
     col  = alpha("navy", 0.5))
dev.print(pdf, 
          width  = 8, 
          height = 8,
          file   = paste0("Plots/NB-RFcaret_train_v", kParamVersion, "_Scatter",
                         "_p", use.per.data,
                         "_ns", rf.ns, 
                         "_t", rf.t, ".pdf" ))
overlapDensity(sc.rf[tnum == 2],
               sc.rf[tnum == 3],
               gpnames = c("Mutual", "Synergistic"))
dev.print(pdf, 
          width  = 8, 
          height = 8, 
          file   = paste0("Plots/NB-RFcaret_train_v", 
                         kParamVersion, "_Density",
                         "_p", use.per.data,
                         "_ns", rf.ns, 
                         "_t", rf.t, ".pdf"))

# prepare data for area under the curve visualisation (estimation)
ppp  <- sc.rf
lll  <- data.train$bliss > bliss.syn.cutoff  # bliss cutoff for synergy
pred <- prediction(ppp, lll)
perf <- performance(pred, "tpr","fpr")
auc  <- performance(pred, "auc")@y.values[[1]] 
auc
plot(perf, avg = "threshold", colorize = T, lwd = 3, main = "NB-RF caret train ROC")
legend("bottomright", 
       legend = c(paste("Random Forests (AUC=", formatC(auc, digits = 4, format = "f"), ")", sep = "")),
          col = c("red"),
          lty = 1)
dev.print(pdf, 
          width  = 10,
          height = 7,
          file   = paste0("Plots/NB-RFcaret_train_v", kParamVersion, "_ROC",
                          "_p", use.per.data,
                          "_ns", rf.ns,
                          "_t", rf.t, ".pdf"))

perf <- performance(pred, measure = "acc", x.measure = "cutoff")
perf
# get the cutoff for the best accuracy
bestAccInd <- which.max(perf@"y.values"[[1]])
bestMsg    <- print(paste("Best accuracy =",
                          perf@"y.values"[[1]][bestAccInd],
                          " at cutoff = ", round(perf@"x.values"[[1]][bestAccInd], 4),
                          sep=""))
plot(perf,sub=bestMsg)
dev.print(pdf, 
          width  = 7, 
          height = 7, 
          file   = paste0("Plots/NB-RFcaret_train_v",
                          kParamVersion, "_BestAccuracy",
                          "_p", use.per.data,
                          "_ns", rf.ns,
                          "_t", rf.t, ".pdf"))

perf <- performance(pred, "prec", "rec")
plot(perf, colorize = T)
dev.print(pdf, 
          width  = 10,
          height = 7, 
          file   = paste0("Plots/NB-RFcaret_train_v", kParamVersion, "_PrecisionRecall",
                          "_p", use.per.data,
                          "_ns", rf.ns,
                          "_t", rf.t, ".pdf"))

# run analysis on the 1/3 training dataset
proba.rf <- predict(rf.estimate, newdata = data.test, type = "prob")
head(proba.rf)
dim(proba.rf)

tnum <- unclass(data.test$synergy)
sc.rf <- Logit(proba.rf[ ,2])

plot(sc.rf, data.test$bliss, 
     pch  = 16,
     main = "Bliss Independence and Synergy Score", 
     ylab = "Bliss Independence", 
     xlab = "Synery Score", 
     col  = alpha("navy", 0.5))
dev.print(pdf, 
          width  = 8, 
          height = 8,
          file   = paste("Plots/NB-RFcaret_test_v", kParamVersion, "_Scatter",
                         "_p", use.per.data,
                         "_ns", rf.ns, 
                         "_t", rf.t, ".pdf",
                         sep = ""))

overlapDensity(sc.rf[tnum == 2], 
               sc.rf[tnum == 3],
               gpnames = c("Mutual", "Synergistic"))
dev.print(pdf, 
          width  = 8,
          height = 8, 
          file = paste0("Plots/NB-RFcaret_test_v", kParamVersion, "_Density", 
                       "_p", use.per.data, 
                       "_ns", rf.ns, 
                       "_t", rf.t, ".pdf"))


ppp  <- sc.rf
lll  <- data.test$bliss > bliss.syn.cutoff
pred <- prediction(ppp, lll)
perf <- performance(pred, "tpr", "fpr")
auc  <- performance(pred, "auc")@y.values[[1]] 
auc
plot(perf, avg = "threshold", colorize = T, lwd = 3, main = "NB-RF caret test")
legend("bottomright", 
       legend = c(paste("Random Forest (AUC = ", 
                        formatC(auc, digits = 4, format = "f"), ")",
                        sep = "")),
          col = c("red"),
          lty = 1)
dev.print(pdf, 
          width  = 10,
          height = 7,
          file   = paste0("Plots/NB-RFcaret_test_v", kParamVersion, "_ROC",
                          "_p", use.per.data, 
                          "_ns", rf.ns, 
                          "_t", rf.t, ".pdf"))

perf <- performance(pred, measure = "acc", x.measure = "cutoff")
perf
# Now let's get the cutoff for the best accuracy
bestAccInd <- which.max(perf@"y.values"[[1]])
bestMsg    <- print(paste0("Best accuracy = ", perf@"y.values"[[1]][bestAccInd],
                          " at cutoff = ", round(perf@"x.values"[[1]][bestAccInd], 4)))
plot(perf, sub = bestMsg)
dev.print(pdf,
          width  = 7,
          height = 7,
          file   = paste0("Plots/NB-RFcaret_test_v", kParamVersion, "_BestAccuracy",
                          "_p", use.per.data, 
                          "_ns", rf.ns, 
                          "_t", rf.t, ".pdf"))
          
perf <- performance(pred, "prec", "rec")
plot(perf, colorize = T)
dev.print(pdf, 
          width  = 10,
          height =  7, 
          file   = paste0("Plots/NB-RFcaret_test_v", kParamVersion, "_PrecisionRecall",
                    "_p", use.per.data, 
                    "_ns", rf.ns, 
                    "_t", rf.t, ".pdf"))
