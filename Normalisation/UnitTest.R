# UnitTest.R
#
# 2009 Jan Wildenhain
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
# This file contains process function to normalise and calculate
# statistics for the HTS screens used in ChemGRID.org
#
# These functions need access to the ChemGRID MySQL database

setwd("~/Documents/GitHub/chemical-genetic-matrix/Normalisation")

source("../PackageInstallation.R")

Sys.setenv(R_DB_HOST = "")
Sys.setenv(R_DB_NAME = "")
Sys.setenv(R_DB_USER = "")
Sys.setenv(R_DB_PWD = "")

source("ChemGridDataProcessingFunctions.R")

UnitTestScriptVariables <- function() {
  
  d <- data.frame(
    id = "WT972-1st",  # Unique Screen Identifier
    sp = 4896,          # Unique Species Identifier
    lib = "Lopac",   # Unique Libary Identifier
    kZscoreBioactiveThreshold = -4,  # minimum Z Score to be classified as active
    kZscoreToxicThreshold = 1,  # Z Score to be considered highly active
    kPValueThreshold = 0.001,  # P value cutoff to be classified as active 
    kRunLevel = 1,
    kUniqueFilename = "UnitTest"
  )
  return(d)
}

# load test variables
d <- UnitTestScriptVariables()
d

##################################################
# do analysis for a specific library             #
##################################################

distributeControlsToTables(d$lib,d$id,d$sp)
dbset <- exp_normalize(d$lib,d$id,d$sp)
dbset <- get_zscores(dbset)
dbset <- undermine_find_outliers(dbset)
insert_stat2db(dbset)
platestat2db(dbset)
expstat2db(dbset)

#####################
# run clusters      #
#####################

if (RUNLEVEL == 2) {
  lib <- "Arya"
  sp <- 4932
  create_cluster_table(lib,sp)  # create heatmap db table
}

#####################
# generate graphs   #
#####################

# raw plot
# normalized plot
# distribution
# correlation

# Lopac
dbset.A <- exp_normalize("Lopac","Tep1","4932")
dbset.A <- get_zscores(dbset.A)
dbset.A
# insert_stat2db(dbset.A)

dbset.B <- exp_normalize("Lopac","Gim3",4932)
dbset.B <- get_zscores(dbset.B)
dbset.A.sam <- get_sam_score(dbset.A,dbset.B)


##############################
# build objects for analysis #
##############################

# positve cycloheximide control
w.r1.pc <- GetCorrectionValues("Lopac","WT",4932,1,1)
w.r1.pc
# negative dmso control
w.r1.nc <- GetCorrectionValues("Lopac","WT",4932,1,1)
w.r1.nc

# OD values 
wildtype <- get_norm_values("Lopac","WT",4932,1,w.r1.nc,w.r1.pc)


