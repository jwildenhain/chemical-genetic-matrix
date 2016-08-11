# HeatmapLibrary.R
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
# Script creates cryptagen activity heatmap for each library
#
library(gplots)

kYourUserFolder <- "~/Documents/GitHub/chemical-genetic-matrix/Heatmaps"

# Load Yeast Bioactive 1 library
library        <- "Bioactive1"
heatmap.width  <- 26
heatmap.height <- 48
filename       <- paste0("cgm_cryptagen_",library,".csv")
source("cgm_Heatmap.R")

# Load Yeast Bioactive 2 library
library        <- "Bioactive2"
heatmap.width  <- 26
heatmap.height <- 48
filename       <- paste0("cgm_cryptagen_",library,".csv")
source("cgm_Heatmap.R")

# Load Maybridge Hitskit 1000 library
library        <- "Hitskit"
heatmap.width  <- 26
heatmap.height <- 48
filename       <- paste0("cgm_cryptagen_",library,".csv")
source("cgm_Heatmap.R")

# Load Sigma Lopac library
library        <- "Lopac"
heatmap.width  <- 14
heatmap.height <- 22
filename       <- paste0("cgm_cryptagen_",library,".csv")
source("cgm_Heatmap.R")

# Load Microsource Spectrum library
library        <- "Spectrum"
heatmap.width  <- 36
heatmap.height <- 58
filename       <- paste0("cgm_cryptagen_",library,".csv")
source("cgm_Heatmap.R")
