# cgm_Heatmap.R
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
# histogram and heatmap for figures
#
xx <- read.table(filename, row.names = 1, header = TRUE, sep= ",")
head(xx[1:5], )
colnames(xx)
dim(xx)

# keep cryptagens 
avec <- rep(0,nrow(xx))
for (i in 1:nrow(xx))  {
  avec[i] <- sum(xx[i,] < -4)
}

cryptagen <- which(avec > 4 & avec <= (0.66 * ncol(xx)))
length(cryptagen)
xy <- as.matrix(xx[cryptagen, ])

# draw histograms and heatmaps with cryptagen
hist(as.numeric(xy), breaks=50, xlab="Z-Score", main=library, col="grey")
dev.print(pdf, 
          width  = 5,
          height = 5,
          file   = paste0("Plots/Histogram_", library, ".pdf"))

# cleanup to improve visual appearance
xy[xy < -35] <- -35
xy[xy > 35] <- 35
colors = c(seq(-35, -5, length=10), seq(-5, 5, length = 10), seq(5, 20, length = 10))
colors
length(colors)
my_palette <- colorRampPalette(c("blue", "gray95", "darkgreen"))(n = 29)

heatmap.2(xy, 
          trace        = "none",
          col          = my_palette,
          symm         = F,
          symkey       = F,
          symbreaks    = T,
          density.info = 'none', 
          scale        = "none",
          main         = library) 
dev.print(pdf, 
          width  = heatmap.width,
          height = heatmap.height,
          file   = paste0("Plots/Heatmap_", library, ".pdf"))
