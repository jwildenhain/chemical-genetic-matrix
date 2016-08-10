# LoadConfigProperties.R
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
# This property file contains the run information for each verion number of 
# per data processing. To showcase the analysis procedure only one run is 
# represented here to obtain numbers for the different simulations performed 
# to generate the used in the manuscript:
# 
# Prediction of Synergism from Chemical-Genetic Interactions by Machine Learning
# Wildenhain et al. Cell Systems 2015
# 

# working directory
kYourUserFolder  <- "~/Documents/GitHub/Chemical_Genetic_Matrix/Prediction"
kParamVersion <- "01"

# Bliss Independence cutoff
bliss.syn.cutoff <- 0.25  # cutoff for synergistic combinations
bliss.ant.cutoff <- -0.18  # cutoff for antagonistic combinations

# parameters for RF learner
set.seed(16543)
use.per.data <- .33  # set size of the training data set
rf.version   <- kParamVersion  # track data simulation version 
rf.ns        <- 14  # node size (tree depth) 
rf.t         <- 512  # number of trees
rf.mtry      <- 14  # number of variables to try at each split point 
rf.grid      <- expand.grid(mtry = c(3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43))

#
# Functions
#

Logit <- function(p, offset = 0.001) log((p + offset)/(1 + offset - p))
