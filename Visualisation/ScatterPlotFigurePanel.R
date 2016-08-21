# ScatterPlotFigurePanel.R
#
# 2008 Jan Wildenhain
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

setwd("~/Documents/GitHub/chemical-genetic-matrix/Visualisation")
source("../PackageInstallation.R")

# database access
Sys.setenv(R_DB_HOST = "")
Sys.setenv(R_DB_NAME = "")
Sys.setenv(R_DB_USER = "")
Sys.setenv(R_DB_PWD = "")

source("../Normalisation/ChemGridDataProcessingFunctions.R")

# set plot type:
# 1:  Create scatter plot (plates and screens)
# 2:  Creates histogram and boxplot of plates and screens 
# 3:  Do image plot of the RAW data plate
# 4:  Do replicate plot of technical/biological replicates (screen)
# 5:  Plot similar compounds with strain bioactivity profiles in different colours
# 6:  Plot raw data for each replicate in a screen
# 7:  Do scatter plot between two screens, reading from the database (+ tooltips)   
# 8:  Plot single compound with strain bioactivity profiles
# 20: Do scatter plot between two screens calculates from raw data 

Sys.setenv(R_DATA_SET_NUM = "")  
# set compound ID if you like to plot a single 
# compoud across different sentinel strains
# R_DATA_SET_NUM has to be  
Sys.setenv(R_COMPOUND  = "")  # compound identifier
Sys.setenv(R_STRAIN    = "")  # unique screen identifier 
Sys.setenv(R_LIBRARY   = "")  # library identifier
Sys.setenv(R_SPECIES   = "")  # species identifier
Sys.setenv(R_REPLICATE = "")  # replicate identifier
Sys.setenv(R_P_NUM     = "0")  # plate identifier 0 = all plates in a library
# comparisions
Sys.setenv(R_STRAIN_2  = "") # used in screen A vs B
Sys.setenv(R_SPECIES_2 = "") # used in screen A vs B 
# file storage
Sys.setenv(R_PLOT_FILE_NAME = "") # file name
Sys.setenv(R_PLOT_FILE_SPEC = "") # file name addition (suffix)
# statistic cutoffs
Sys.setenv(R_Z_BIO_THRHOLD  = "-5") # Z score cutoff for bioactivity
Sys.setenv(R_Z_TOX_THRHOLD  = "-44") # Z-Score cutoff to define highly active compounds
Sys.setenv(R_ALPHA_THRHOLD  = "0.001") # Pvalue cutoff based on N(1,sd) fit
# visualisation specifics
Sys.setenv(R_PLOT_TYPE            = "pdf")  # pdf, svg, ps, csv (data sheet), screen (svg inlet on website)
Sys.setenv(R_DATA_PLOT_COLUMN     = "z_score")  # show z_score column value, z_score, z (z-factor) etc.
Sys.setenv(R_DB_ORDERBY_STATEMENT = " a.plate_number,a.plate_row,a.plate_column ")
Sys.setenv(R_SVG_HYPERLINK = "http://chemgrid.org/cgm/tmp_compound.php?")

########################################################
# show representative ChemGRID.org data visualisations #
########################################################

# make comparison scatterplot with recalulation
# passing identifiers to generate plot
Sys.setenv(R_STRAIN    = "RPL9B")    # unique screen identifier 
Sys.setenv(R_LIBRARY   = "Spectrum") # library identifier
Sys.setenv(R_SPECIES   = "4932")     # species identifier
Sys.setenv(R_STRAIN_2  = "MT1448-2") # used in screen A vs B
Sys.setenv(R_SPECIES_2 = "4932")     # used in screen A vs B
Sys.setenv(R_DATA_SET_NUM = "20")    # plot type DSN 
Sys.setenv(R_PLOT_FILE_NAME = "Plots/TestPlotDSN20") # output file
source("PlotsOnChemGRID.R")

# show comparision as svg
Sys.setenv(R_DATA_SET_NUM = "7") 
Sys.setenv(R_PLOT_FILE_NAME = "Plots/TestPlotDSN7")
source("PlotsOnChemGRID.R")

# show scatter plot between technical/biological replicates 
Sys.setenv(R_DATA_SET_NUM = "4")
Sys.setenv(R_PLOT_FILE_NAME = "Plots/TestPlotDSN4")
source("PlotsOnChemGRID.R")

# show raw data example
Sys.setenv(R_DATA_SET_NUM = "6")
Sys.setenv(R_PLOT_FILE_NAME = "Plots/TestPlotDSN6")
source("PlotsOnChemGRID.R")

# show normalised data example
Sys.setenv(R_DATA_SET_NUM = "1")
Sys.setenv(R_PLOT_FILE_NAME = "Plots/TestPlotDSN1")
source("PlotsOnChemGRID.R")

# show histogram/boxplot example
Sys.setenv(R_DATA_SET_NUM = "2")
Sys.setenv(R_PLOT_FILE_NAME = "Plots/TestPlotDSN2")
source("PlotsOnChemGRID.R")

# show single compounds example
Sys.setenv(R_COMPOUND  = "LOPAC 00608")
Sys.setenv(R_DATA_SET_NUM = "8")
Sys.setenv(R_PLOT_FILE_NAME = "Plots/TestPlotDSN8")
source("PlotsOnChemGRID.R")

# show similar compounds example
Sys.setenv(R_COMPOUND  = "LOPAC 00608")
Sys.setenv(R_DATA_SET_NUM = "5")
Sys.setenv(R_PLOT_FILE_NAME = "Plots/TestPlotDSN5")
source("PlotsOnChemGRID.R")

