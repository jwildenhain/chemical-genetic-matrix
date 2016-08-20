# ShowSentinelZFactors.R
#
# 2016 Jan Wildenhain
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
#
library(DBI)
library(RMySQL)
library(ggplot2)

kDB_HOST <- Sys.getenv("R_DB_HOST")
kDB_NAME <- Sys.getenv("R_DB_NAME")
kDB_USER <- Sys.getenv("R_DB_USER")
kDB_PWD  <- Sys.getenv("R_DB_PWD")

kMyDBHandShake <- dbConnect(MySQL(), user = kDB_USER, dbname = kDB_NAME, host = kDB_HOST, password = kDB_PWD)

sql  <- "SELECT * FROM `res_exp`"
rs   <- dbSendQuery(kMyDBHandShake, sql)
data <- fetch(rs, n = -1) # get all datapoints

data$lib[data$lib == 'Spectrum'] <- "Spectrum1"
data$lib[data$lib == 'SPECMTS3'] <- "Spectrum2"
data$lib[data$lib == 'Spectrum_ED'] <- "Spectrum3"
data$lib[data$lib == 'Maybridge1000'] <- "Maybridge"
data$lib[data$lib == 'Bioactive'] <- "Bioactive1"
data$lib[data$lib == 'Cytotoxic'] <- "Bioactive2"
data$lib

df <- data.frame()

df1 <-
  data.frame(lib = rep("Lopac", length(data$z[data$lib == 'Lopac' &
                                                data$z > 0])), data = data$z[data$lib == 'Lopac' & data$z > 0])
df2 <-
  data.frame(lib = rep("Spectrum1", length(data$z[data$lib == 'Spectrum1' &
                                                   data$z > 0])), data = data$z[data$lib == 'Spectrum1' & data$z > 0])
df3 <-
  data.frame(lib = rep("Spectrum2", length(data$z[data$lib == 'Spectrum2' &
                                                    data$z > 0])), data = data$z[data$lib == 'Spectrum2' & data$z > 0])
df4 <-
  data.frame(lib = rep("Spectrum3", length(data$z[data$lib == 'Spectrum3' &
                                                      data$z > 0])), data = data$z[data$lib == 'Spectrum3' &
                                                                                     data$z > 0])
df5 <-
  data.frame(lib = rep("Maybridge", length(data$z[data$lib == 'Maybridge' &
                                                    data$z > 0])), data = data$z[data$lib == 'Maybridge' &
                                                                                   data$z > 0])
df6 <-
  data.frame(lib = rep("Bioactive1", length(data$z[data$lib == 'Bioactive1' &
                                                    data$z > 0])), data = data$z[data$lib == 'Bioactive1' & data$z > 0])
df7 <-
  data.frame(lib = rep("Bioactive2", length(data$z[data$lib == 'Bioactive2' &
                                                    data$z > 0])), data = data$z[data$lib == 'Bioactive2' & data$z > 0])
df <- rbind(df1, df2, df3, df4, df5, df6, df7)
df

p <- ggplot(df, aes(factor(lib), data))

p + geom_violin(fill = "grey70", colour = "grey70") + 
  geom_jitter(height = 0) + ggtitle("Sentinel HTS Quality Assessment") +
  labs(x = "Library", y = "Z-Factor") 
