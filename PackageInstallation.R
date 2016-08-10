# PackageInstallation.R
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
# This file installs the packages required to run all of the 
# Chemical Genetic Matrix R code

if (!require("pacman")) install.packages("pacman")
pacman::p_load(caret, 
               pROC, 
               randomForest, 
               FactoMineR, 
               DAAG, 
               ggplot2, 
               ROCR, 
               plyr, 
               party, 
               e1071)

# load packages required for scripts
library(gplots)
library(RColorBrewer)
library(rpart)
library(randomForest)
library(MASS)
library(splines)
library(DAAG)
library(FactoMineR)
library(caret)
library(ggplot2)
library(ROCR)
library(party)