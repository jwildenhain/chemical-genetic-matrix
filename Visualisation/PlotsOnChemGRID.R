# PlotsOnChemGRID.R
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
# This file contains visualisation functions used on ChemGRID.org
#
# This file will require access to the ChemGRID MySQL database

setwd("~/Documents/GitHub/chemical-genetic-matrix/Visualisation")

# Settings transfered as system variables from PHP to R cmd line
dsn      <- as.numeric(Sys.getenv("R_DATA_SET_NUM"))
cid      <- Sys.getenv("R_COMPOUND")
strain   <- (Sys.getenv("R_STRAIN"))
lib      <- (Sys.getenv("R_LIBRARY"))
sp       <- (Sys.getenv("R_SPECIES"))
replicate <- as.numeric(Sys.getenv("R_REPLICATE"))
pn        <- as.numeric(Sys.getenv("R_P_NUM"))
# File settings
adm_fpath <- (Sys.getenv("R_ADM_PATH"))
uid_name  <- (Sys.getenv("R_PLOT_FILE_NAME"))
uid_spec  <- (Sys.getenv("R_PLOT_FILE_SPEC"))
# Statistic cutoff settings
z_bio_thr <- as.numeric(Sys.getenv("R_Z_BIO_THRHOLD"))
z_tox_thr <- as.numeric(Sys.getenv("R_Z_TOX_THRHOLD"))
z_pval_thr <- as.numeric(Sys.getenv("R_ALPHA_THRHOLD"))
kPlotOutputType       <- (Sys.getenv("R_PLOT_TYPE")) # PDF, SVG, PS, SCREEN
# show Z_score column value, z_score, z (z-factor) etc.
kShowDataColumn   <- (Sys.getenv("R_DATA_PLOT_COLUMN")) 
kOrderByStatement <- (Sys.getenv("R_DB_ORDERBY_STATEMENT"))
kRSVGHyperlink    <- Sys.getenv("R_SVG_HYPERLINK")

library("RSvgDevice")  # used initially in 2005 
library("RSVGTipsDevice")  # was introduced later 2007 and used from then
library(DBI)
library(RMySQL)

DB_HOST <- Sys.getenv("R_DB_HOST")
DB_NAME <- Sys.getenv("R_DB_NAME")
DB_USER <- Sys.getenv("R_DB_USER")
DB_PWD  <- Sys.getenv("R_DB_PWD")

mycon <-
  dbConnect(
    MySQL(),
    user     = DB_USER,
    dbname   = DB_NAME,
    host     = DB_HOST,
    password = DB_PWD
  )

do_replicate_plot <- function(dbset, plcount, modlabel, viewcol) {
  # show plot between biological/technical replciates
  # Args: 
  #   dbset: data list object
  #   plcount: replicate ID
  #   modlabel: label for axis
  #   viewcol: selected data column
  # Returns: 
  #   Scatter plot
  if (viewcol == "z_score") {
    dbset$av <- dbset$z_score
  }
  
  replicateA <- paste("Replicate ", (plcount - 1), "", sep = "")
  replicateB <- paste("Replicate ", plcount, "", sep = "")
  title <-
    paste(
      "OD Scores for ",
      noquote(dbset$strain),
      " ",
      noquote(dbset$library),
      " (",
      noquote(dbset$species),
      ")",
      sep = ""
    )
  plot(
    dbset$av[, (plcount - 1)],
    dbset$av[, plcount],
    xlab = replicateA,
    ylab = replicateB,
    main = title
  )
  
  fontsize <- 0.5
  
  for (i in 1:dbset$datapoints) {
    if (modlabel) {
      plotlabel <-
        paste(dbset$db[i, 1], sep = "")
    } else {
      plotlabel <-
        paste(dbset$db[i, 1], dbset$strain, dbset$species, sep = ";")
    }
    mypos <- 3
    if ((dbset$av[i, (plcount - 1)] / dbset$av[i, plcount]) < 1) {
      mypos <- 2
    } else {
      mypos <- 4
    }
    if (dbset$av[i, (plcount - 1)] > (max(dbset$av[, (plcount - 1)]) - 0.15)) {
      mypos <- 2
    }
    if (dbset$av[i, plcount] < (min(dbset$av[, plcount]) + 0.3)) {
      mypos <- 4
    }
    if (dbset$av[i, (plcount - 1)] < (min(dbset$av[, (plcount - 1)]) + 0.3)) {
      mypos <- 4
    }
    if ((dbset$oav[i, 1] < 1) &
        (
          (dbset$z_score[i, (plcount - 1)] < z_bio_thr &
           dbset$z_score[i, (plcount - 1)] >= z_tox_thr) &
          (dbset$z_score[i, plcount] < z_bio_thr &
           dbset$z_score[i, plcount] >= z_tox_thr)
        ) &
        (dbset$z_pvalue[i, (plcount - 1)] < z_pval_thr &
         dbset$z_pvalue[i, plcount] < z_pval_thr)) {
      points(dbset$av[i, (plcount - 1)], dbset$av[i, plcount], pch = 20, col =
               "blue")
      text(dbset$av[i, (plcount - 1)],
           dbset$av[i, plcount],
           cex = fontsize,
           pos = mypos ,
           plotlabel)
    }
    else if (dbset$oav[i, 1]) {
      points(dbset$av[i, (plcount - 1)], dbset$av[i, plcount], pch = 20, col =
               "red")
      text(dbset$av[i, (plcount - 1)],
           dbset$av[i, plcount],
           cex = fontsize,
           pos = mypos ,
           plotlabel)
    }
    else if (dbset$z_score[i, (plcount - 1)] < z_tox_thr &
             dbset$z_score[i, plcount] < z_tox_thr) {
      points(dbset$av[i, (plcount - 1)], dbset$av[i, plcount], pch = 20, col =
               "yellow")
      text(dbset$av[i, (plcount - 1)],
           dbset$av[i, plcount],
           cex = fontsize,
           pos = mypos ,
           plotlabel)
    }
    else if (dbset$z_score[i, (plcount - 1)] > abs(z_bio_thr) &
             dbset$z_score[i, plcount] > abs(z_bio_thr)) {
      points(dbset$av[i, (plcount - 1)], dbset$av[i, plcount], pch = 20, col =
               "green")
      text(dbset$av[i, (plcount - 1)],
           dbset$av[i, plcount],
           cex = fontsize,
           pos = mypos ,
           plotlabel)
    }
    else if ((dbset$z_score[i, (plcount - 1)] < z_tox_thr |
              dbset$z_score[i, plcount] < z_tox_thr) &
             dbset$z_score[i, (plcount - 1)] < z_bio_thr &
             dbset$z_score[i, plcount] < z_bio_thr) {
      points(dbset$av[i, (plcount - 1)], dbset$av[i, plcount], pch = 20, col =
               "blue")
      text(dbset$av[i, (plcount - 1)],
           dbset$av[i, plcount],
           cex = fontsize,
           pos = mypos ,
           plotlabel)
    }
    
  }
  
}


do_betw_replicate_plot2 <- function(kPlotOutputType) {
  # produce A vs. B scatter plot
  # Args:
  #   kPlotOutputType: svg or csv
  strainA <- (Sys.getenv("R_STRAIN"))
  strainB <- (Sys.getenv("R_STRAIN_2"))
  spA <- (Sys.getenv("R_SPECIES"))
  spB <- (Sys.getenv("R_SPECIES_2"))
  lib <- (Sys.getenv("R_LIBRARY"))
  cid <- Sys.getenv("R_COMPOUND")

  kOrderByStatement <- " a.plate_number,a.plate_row,a.plate_column "
  dbtablecolumn     <- paste(strainA, "versus", strainB)
  filename <-
    paste0(noquote(uid_name), noquote(uid_spec), ".svg")
  
  # add queries go get dataset a and dataset b, add link column to direct single dataset links
  sqlstr <-
    paste0(
      "SELECT a.supplier_obj_id, a.value, b.value, a.z_score, b.z_score, a.z_factor, b.z_factor, a.per_inh, b.per_inh, a.outlier, b.outlier FROM `Z_norm` a join `Z_norm` b on a.species='",
      noquote(spA),
      "' and a.strain='",
      noquote(strainA),
      "' and b.species='",
      noquote(spB),
      "' and b.strain='",
      noquote(strainB),
      "' and a.library=b.library and a.plate_row =b.plate_row and a.plate_column = b.plate_column and a.plate_number = b.plate_number and a.library ='",
      noquote(lib),
      "' order by ",
      noquote(kOrderByStatement)
    )
  
  rs        <- dbSendQuery(mycon, sqlstr)
  res       <- fetch(rs, n = -1)
  df.length <- nrow(res)
  df.color  <- rep("grey", df.length)
  df.color[which(res[, 10] == 1 | res[, 11] == 1)] <- "red"
  df.sign   <- rep(19, df.length)
  mydat     <- data.frame(res, df.color, df.sign)
  
  devSVGTips(
    filename,
    height = 7,
    width = 8,
    toolTipMode = 1,
    toolTipFontSize = 8,
    title = "Plate test, tooltips are title + 1 line",
    xmlHeader = TRUE
  )
  if (kShowDataColumn  == "value") {
    colA   <- 2
    colB   <- 3
    xlabel <- paste("Normalised Value", strainA)
    ylabel <- paste("Normalised Value", strainB)
    xaxis  <- c(0, 1.1)
    yaxis  <- c(0, 1.1)
  } else if (kShowDataColumn  == "per_inh") {
    colA   <- 8
    colB   <- 9
    xlabel <- paste("Percent Inhibition", strainA)
    ylabel <- paste("Percent Inhibition", strainB)
    xaxis  <- c(0, max(mydat[, colA]))
    yaxis  <- c(0, max(mydat[, colB]))
  } else if (kShowDataColumn  == "z_score") {
    colA   <- 4
    colB   <- 5
    xlabel <- paste("Z Score", strainA)
    ylabel <- paste("Z Score", strainB)
    xaxis  <- c(min(mydat[, colA]), max(mydat[, colA]))
    yaxis  <- c(min(mydat[, colA]), max(mydat[, colB]))
  } else if (kShowDataColumn  == "z_factor") {
    colA   <- 6
    colB   <- 7
    xlabel <- paste("Z Factor", strainA)
    ylabel <- paste("Z Factor", strainB)
    xaxis  <- c(-5, 1)
    yaxis  <- c(-5, 1)
  }
  
  par(mar = c(8, 8, 10, 10), cex = .75)
  plot(
    mydat[, colA],
    mydat[, colB],
    pch  = (as.numeric(mydat$df.sign)),
    col  = as.character(mydat$df.color),
    type = "n",
    bty  = "n",
    lwd  = 2,
    xlab = xlabel,
    ylab = ylabel,
    xlim = xaxis,
    ylim = yaxis,
    main = dbtablecolumn
  )
  
  for (i in 1:df.length) {
    setSVGShapeToolTip(title = mydat[i, "supplier_obj_id"], desc = paste(mydat[i, colA]))
    
    #symbol=ARNT2&row=F&col=17&p=100004208&library=dharmacon
    setSVGShapeURL(
      paste(
        kRSVGHyperlink,
        "cid=",
        mydat[i, "supplier_obj_id"],
        "&s=",
        strainA,
        "&sp=",
        spA,
        sep = ""
      ),
      target = "_top"
    )
    if (cid == res[i, "supplier_obj_id"]) {
      points(
        mydat[i, colA],
        y   = mydat[i, colB],
        pch = (as.numeric(mydat[i, "df.sign"])),
        cex = 2,
        col = "black"
      )
    } else {
      points(
        mydat[i, colA],
        y   = mydat[i, colB],
        pch = (as.numeric(mydat[i, "df.sign"])),
        cex = 1,
        col = (as.character(mydat[i, "df.color"]))
      )
    }
  }

  dev.off()
  
  if (kPlotOutputType == "csv") {
    filename <- paste0(noquote(uid_name), ".csv")
    x        <- data.frame(mydat)
    write.table(
      x,
      file = filename,
      sep  = ",",
      row.names = FALSE,
      qmethod   = "double"
    )
  }
} 


do_betw_replicate_plot <- function(dbsetA, dbsetB, modlabel) {
  # do replicate plot between two screens
  # Args:
  #  dbsetA: data object A
  #  dbsetB: data object B
  #  modlabel: label for the plot
  # Returns:
  #  Plot pdf, svg, ...
  replicateA <- paste("Strain ", dbsetA$strain, "", sep = "")
  replicateB <- paste("Strain ", dbsetB$strain, "", sep = "")
  title <-
    paste(
      kShowDataColumn ,
      " for ",
      noquote(dbsetA$strain),
      " vs. ",
      noquote(dbsetB$strain),
      " ",
      noquote(dbset$library),
      sep = ""
    )
  
  Aav <- apply(dbsetA$av, 1, median, na.rm = T)
  Bav <- apply(dbsetB$av, 1, median, na.rm = T)
  Az  <- apply(dbsetA$z_score, 1, median, na.rm = T)
  Bz  <- apply(dbsetB$z_score, 1, median, na.rm = T)
  Ap  <- apply(dbsetA$z_pvalue, 1, median, na.rm = T)
  Bp  <- apply(dbsetB$z_pvalue, 1, median, na.rm = T)
  plot(Aav,
       Bav,
       xlab = replicateA,
       ylab = replicateB,
       main = title)
  
  fontsize <- 0.5
  
  for (i in 1:dbset$datapoints) {
    if (modlabel) {
      plotlabel <- "" # paste(dbset$db[i,1],sep="")
    } else {
      plotlabel <-
        paste(dbset$db[i, 1], dbset$strain, dbset$species, sep = ";")
    }
    mypos   <- 3
    if ((Aav[i] / Bav[i]) < 1) {
      mypos <- 2
    } else {
      mypos <- 4
    }
    if (Aav[i] > (max(Aav) - 0.15)) {
      mypos <- 2
    }
    if (Bav[i] < (min(Bav) + 0.3)) {
      mypos <- 4
    }
    if (Aav[i] < (min(Aav) + 0.3)) {
      mypos <- 4
    }
    if ((dbsetA$oav[i, 1] == 0) & (dbsetB$oav[i, 1] == 0) &
        ((Az[i] <= z_bio_thr &
          Az[i] >= z_tox_thr) |
         (Bz[i] < z_bio_thr &
          Bz[i] >= z_tox_thr)) &
        (Az[i] >= z_tox_thr & Bz[i] >= z_tox_thr) &
        (Ap[i] < z_pval_thr |
         Bp[i] < z_pval_thr)) {
      points(Aav[i], Bav[i], pch = 20, col = "blue")
      text(Aav[i], Bav[i], cex = fontsize, pos = mypos , plotlabel)
    }
    else if (dbset$oav[i, 1] == 1 |
             dbsetB$oav[i, 1] == 1 |
             ((Az[i] <= z_bio_thr &
               Ap[i] > z_pval_thr) |
              (Bz[i] <= z_bio_thr & Bp[i] > z_pval_thr))) {
      points(Aav[i], Bav[i], pch = 20, col = "red")
      text(Aav[i], Bav[i], cex = fontsize, pos = mypos , plotlabel)
    }
    else if (Az[i] < z_tox_thr &
             Bz[i] < z_tox_thr & Ap[i] <= z_pval_thr & Bp[i] <= z_pval_thr) {
      points(Aav[i], Bav[i], pch = 20, col = "yellow")
      text(Aav[i], Bav[i], cex = fontsize, pos = mypos , plotlabel)
    }
    else if (((Az[i] < z_tox_thr &
               Bz[i] < z_bio_thr) |
              (Az[i] < z_bio_thr &
               Bz[i] < z_tox_thr)) &
             Ap[i] <= z_pval_thr & Bp[i] <= z_pval_thr) {
      points(Aav[i], Bav[i], pch = 20, col = "orange")
      text(Aav[i], Bav[i], cex = fontsize, pos = mypos , plotlabel)
    }
    else if ((Az[i] < z_tox_thr |
              Bz[i] < z_tox_thr)  &
             (Ap[i] <= z_pval_thr | Bp[i] <= z_pval_thr)) {
      points(Aav[i], Bav[i], pch = 20, col = "magenta")
      text(Aav[i], Bav[i], cex = fontsize, pos = mypos , plotlabel)
    }
    else if (Az[i] > abs(z_bio_thr) &
             Bz[i] > abs(z_bio_thr) &
             Ap[i] <= z_pval_thr & Bp[i] <= z_pval_thr) {
      points(Aav[i], Bav[i], pch = 20, col = "green")
      text(Aav[i], Bav[i], cex = fontsize, pos = mypos , plotlabel)
    }
    else if ((Az[i] < z_tox_thr |
              Bz[i] < z_tox_thr) &
             Az[i] < z_bio_thr &
             Bz[i] < z_bio_thr & Ap[i] > z_pval_thr & Bp[i] > z_pval_thr) {
      points(Aav[i], Bav[i], pch = 20, col = "blue")
      text(Aav[i], Bav[i], cex = fontsize, pos = mypos , plotlabel)
    }
    
  }
  
}



plot_raw_data <- function(dbset, plnr, modlabel) {
  title <-
    paste(
      kShowDataColumn ,
      " for ",
      noquote(dbset$strain),
      " ",
      noquote(dbset$library),
      " (",
      noquote(dbset$species),
      ")",
      sep = ""
    )
  plot(
    1:dbset$datapoints,
    dbset$av[, plnr],
    xlab = "Datapoints",
    ylab = "OD Measurement",
    main = title
  )
  
}


set_plot_options <- function(kShowDataColumn , data) {
  # adjustments to optimize axis for visualisation
  # Args: 
  #   kShowDataColumn: MySQL data column
  #   data: the data object
  # Returns: 
  #   Dataframe with x and y axis ranges
  if (kShowDataColumn  == "per_inh") {
    ylabel <- paste("Percent Inhibition")
    maxpoint <- 100
    if (max(data) < 100) {
      maxpoint <- 100
    }
    else {
      maxpoint <- max(data)
    }
    xaxis <- c(0, max(data))
    yaxis <- c(0, max(data))
    xadj <- c(2, 16)
  } else if (kShowDataColumn  == "z_factor") {
    ylabel <- paste("Z Factor")
    xaxis <- c(-2, 1)
    yaxis <- c(-2, 1)
    xadj <- c(0.1, 1)
  } else if (kShowDataColumn  == "value") {
    ylabel <- paste("Normalised Measurement (OD)")
    maxpoint <- 1
    if (max(data) < 1) {
      maxpoint <- 1
    }
    else {
      maxpoint <- max(data)
    }
    xaxis <- c(min(data), max(data))
    yaxis <- c(min(data), maxpoint)
    xadj <- c(0.04, 0.3)
  } else if (kShowDataColumn  == "z_score") {
    ylabel <- paste("Z Score")
    xaxis <- c(min(data), max(data))
    yaxis <- c(min(data), max(data))
    xadj <- c(2, 16)
  }
  return(data.frame(
    yLabel = ylabel,
    xRange = xaxis,
    yRange = yaxis,
    xLabAdj = xadj
  ))
}


plot_similar_compounds2 <- function(dbset, cid, pdfps_plot) {
  # Plots compound identifier across several strains that have 90%
  # structural similarity
  # Args:
  #   dbset: data object
  #   cid: compound identifier
  #   pdfps_plot: plot type 
  #num_rowcols <- dim(dbset)
  kRSVGHyperlink <- "http://chemgrid.org/cgm/tmp_compound.php?"
  str <- factor(dbset$strain)
  sup <- factor(dbset$supplier_obj_id)
  ptitle <- paste("", noquote(cid), " and analogs", sep = "")
  dbset$s <- as.integer(str)
  dbset$c <- as.integer(sup)
  
  #    sqlstring <- paste("SELECT species, strain, value, z_score, p_value, outlier, library FROM Z_norm, may_keys WHERE supplier_obj_id = '",noquote(cid),"' and b=supplier_obj_id order by species, strain, supplier_obj_id ",sep="")
  #    rs <- dbSendQuery(mycon, sqlstring)
  #    dbset <- fetch(rs, n = -1) # get datapoints per library
  # Set margins to make room for x axis labels
  #par(mar = c(12, 4, 4, 2) + 0.1)
  # plot
  #par(cex.axis=0.1)
  par(mar = c(15, 3.8, 4.2, 2), cex = .75)
  SPO <- set_plot_options(kShowDataColumn , dbset$value)
  plot(
    value ~ value,
    xaxt = "n",
    type = "n",
    bty = "n",
    title = ptitle,
    xlab = "",
    ylim = SPO$yRange,
    ylab = SPO$yLabel[1],
    data = dbset,
    pch = 20,
    col = c,
    cex = 2,
    ann = T
  )
  
  test = axis(1, at = seq(1, nrow(dbset)), lab = F)
  if (nrow(dbset) <= 30) {
    text(
      seq(1, nrow(dbset)),
      rep(par("usr")[3] - SPO$xLabAdj[1], nrow(dbset)),
      srt = 45,
      adj = 1,
      labels = dbset$strain,
      xpd = T,
      cex = 0.7
    )
    text(
      nrow(dbset) / 2,
      par("usr")[3] - SPO$xLabAdj[2],
      adj = 1,
      labels = "Experiments",
      xpd = T,
      cex = 1
    )
  } else {
    text(
      seq(1, nrow(dbset)),
      rep(par("usr")[3] - SPO$xLabAdj[1], nrow(dbset)),
      srt = 45,
      adj = 1,
      labels = dbset$strain,
      xpd = T,
      cex = 0.9
    )
    text(
      nrow(dbset) / 2,
      par("usr")[3] - SPO$xLabAdj[2],
      adj = 1,
      labels = "If there are too many experiments, please use tooltip for identification.",
      xpd = T,
      cex = .7,
      col = "red"
    )
  }
  #axis(1, at = seq(1,nrow(dbset)), srt=45, adj=1, xpd=T, cex=0.5, labels = dbset$strain)
  title(
    ptitle,
    cex.main = 1.5,
    font.main = 1,
    col.main = "darkgrey"
  )
  for (i in 1:length(dbset$value)) {
    setSVGShapeToolTip(title = sup[i],
                       desc = paste(dbset$value[i], "in", str[i]))
    
    #symbol=ARNT2&row=F&col=17&p=100004208&library=dharmacon
    setSVGShapeURL(
      paste(
        kRSVGHyperlink,
        "cid=",
        sup[i],
        "&s=",
        str[i],
        "&sp=",
        dbset$species[i],
        sep = ""
      ),
      target = "_top"
    )
    #if(cid == sup[i]) {
    #         points(i,y=dbset$value[i], col=((colvec[i])), pch=20,cex=2)
    #} else {
    #         points(i,y=dbset$value[i], col=((colvec[i])), pch=20,cex=1 )
    #}
    if (dbset$outlier[i] == 0 &&
        dbset$z_score[i] <= (z_tox_thr + dbset$tox[i]) &&
        dbset$p_value[i] <= z_pval_thr) {
      points(
        i,
        y = dbset$value[i],
        pch = 20,
        cex = 2,
        col = "gold"
      )
    } else if (dbset$outlier[i] == 0 &&
               dbset$z_score[i] <= z_bio_thr && dbset$p_value[i] <= z_pval_thr) {
      points(
        i,
        y = dbset$value[i],
        pch = 20,
        cex = 2,
        col = "blue"
      )
    } else if (dbset$outlier[i] == 0 &&
               dbset$z_score[i] >= abs(z_bio_thr) &&
               dbset$p_value[i] <= z_pval_thr) {
      points(
        i,
        y = dbset$value[i],
        pch = 20,
        cex = 2,
        col = "green"
      )
    } else if (dbset$outlier[i] ||
               (abs(dbset$z_score[i]) >= abs(z_bio_thr) &&
                dbset$p_value[i] > z_pval_thr)) {
      points(
        i,
        y = dbset$value[i],
        pch = 20,
        cex = 2,
        col = "red"
      )
    } else {
      points(
        i,
        y = dbset$value[i],
        pch = 20,
        cex = 2,
        col = "grey"
      )
    }
  }
  # x labels
  #axis(1:num_rowcols[1], labels = FALSE)
  #text(1:num_rowcols[1], par("usr")[3] - 0.055, srt = 90, adj = 0.5, xpd = TRUE)
  #text(1:num_rowcols[1], par("usr")[3] - 0.055, srt = 90, adj = 0.5, labels = dbset$strain, xpd = TRUE)
  
  return(dbset)
}



plot_similar_compounds <- function(dbset, cid, pdfps_plot) {
  #num_rowcols <- dim(dbset)
  str <- factor(dbset$strain)
  sup <- factor(dbset$supplier_obj_id)
  ptitle <-
    paste(kShowDataColumn , " for ", noquote(cid), " and analogs", sep = "")
  dbset$s <- as.integer(str)
  
  dbset$c <- as.integer(sup)
  
  #    sqlstring <- paste("SELECT species, strain, value, z_score, p_value, outlier, library FROM Z_norm, may_keys WHERE supplier_obj_id = '",noquote(cid),"' and b=supplier_obj_id order by species, strain, supplier_obj_id ",sep="")
  #    rs <- dbSendQuery(mycon, sqlstring)
  #    dbset <- fetch(rs, n = -1) # get datapoints per library
  # Set margins to make room for x axis labels
  par(mar = c(4, 4, 3, 2) + 0.1)
  # plot
  #par(cex.axis=0.1)
  #par(mar=c(15, 3.8, 4.2, 2),cex = .75)
  plot(
    value ~ value,
    xaxt = "n",
    title = ptitle,
    xlab = "Experiments",
    ylab = "OD Values",
    data = dbset,
    col = c
  )
  
  # plot(dbset$value, xaxt = "n", xlab = "Experiments",ylab="OD Values", groups="dbset$supplier_obj_id");
  # legend(0.5,2,group=dbset$supplier_obj_id, pch=c(1,2), col=c("green","brown")) # add legend
  newcol = 0
  ymax <- max(dbset$value)
  fontsize <- 0.5
  for (i in 1:length(dbset$value)) {
    if (pdfps_plot) {
      mystring <- paste(dbset$strain[i], sep = "")
    } else {
      mystring <-
        paste(dbset$supplier_obj_id[i],
              dbset$strain[i],
              dbset$species[i],
              sep = ";")
    }
    mypos <- 3
    if (dbset$value[i] > (ymax * .9)) {
      mypos <- 1
    }
    if (i < length(dbset$value) * .1) {
      mypos <- 4
    }
    if (i > (length(dbset$value) * .9)) {
      mypos <- 2
    }
    if (dbset$outlier[i] == 0 &&
        dbset$z_score[i] <= (z_tox_thr + dbset$tox[i]) &&
        dbset$p_value[i] <= z_pval_thr) {
      points(
        i,
        dbset$value[i],
        pch = 20,
        cex = 1,
        col = "yellow",
        bg = newcol
      )
      text(i,
           dbset$value[i],
           cex = fontsize,
           pos = mypos ,
           mystring)
      #text(i,dbset$value[i], cex=1, pos=mypos, srt = 90, adj = 0.5 ,dbset$strain[i])
    } else if (dbset$outlier[i] == 0 &&
               dbset$z_score[i] <= z_bio_thr &&
               dbset$p_value[i] <= z_pval_thr) {
      points(
        i,
        dbset$value[i],
        pch = 20,
        cex = 1,
        col = "blue",
        bg = newcol
      )
      text(i,
           dbset$value[i],
           cex = fontsize,
           pos = mypos ,
           mystring)
      #text(i,dbset$value[i], cex=1, pos=mypos, srt = 90, adj = 0.5 ,dbset$strain[i])
    } else if (dbset$outlier[i] == 0 &&
               dbset$z_score[i] >= abs(z_bio_thr) &&
               dbset$p_value[i] <= z_pval_thr) {
      points(
        i,
        dbset$value[i],
        pch = 20,
        cex = 1,
        col = "green",
        bg = newcol
      )
      text(i,
           dbset$value[i],
           cex = fontsize,
           pos = mypos ,
           mystring)
      #text(i,dbset$value[i], cex=1, pos=mypos, srt = 90, adj = 0.5 ,dbset$strain[i])
    } else if (dbset$outlier[i] ||
               (abs(dbset$z_score[i]) >= abs(z_bio_thr) &&
                dbset$p_value[i] > z_pval_thr)) {
      points(
        i,
        dbset$value[i],
        pch = 20,
        cex = 1,
        col = "red",
        bg = newcol
      )
      text(i,
           dbset$value[i],
           cex = fontsize,
           pos = mypos ,
           mystring)
      #text(i,dbset$value[i], cex=1, pos=mypos, srt = 90, adj = 0.5 ,dbset$strain[i])
    } else {
      if (!pdfps_plot) {
        points(i, dbset$value[i], pch = 1, col = dbset$c[i])
        text(i,
             dbset$value[i],
             cex = fontsize,
             pos = mypos ,
             mystring)
      }
    }
  }
  #mtext("text to place", 2)
  
  # x labels
  #axis(1:num_rowcols[1], labels = FALSE)
  #text(1:num_rowcols[1], par("usr")[3] - 0.055, srt = 90, adj = 0.5, xpd = TRUE)
  #text(1:num_rowcols[1], par("usr")[3] - 0.055, srt = 90, adj = 0.5, labels = dbset$strain, xpd = TRUE)
  
  return(dbset)
}

plot_norm_scores <- function(dbset, z_bio_thr, z_tox_thr, z_pval_thr, cid, modlabel) {
  # Plot normalised data and highlight current molecule (chemgrid.org)
  # Args:
  #   dbset: data object
  #   z_bio_thr: Z Score cutoffs
  #   z_tox_thr: Highly active cutoff (draws yellow datapoints)
  #   z_pval_thr: Fitted N(0,IQR) p-value cutoff
  #   cid: compound identifier
  #   modlabel: label for plot
    num_rowcols <- dim(dbset)
    #par(cex.axis=0.1)
    SPO <- set_plot_options(kShowDataColumn , dbset$value)
    
    par(mar = c(8, 8, 4.2, 0.2))
    plot(
      1:num_rowcols[1],
      dbset$value,
      ylim = SPO$yRange,
      xlab = "Datapoints",
      ylab = "OD Values",
      main = title
    )
    #lines(lowess(1:num_rowcols[1],dbset$value, f=1/10))   #lowess smooth with 30% window
    text(
      nrow(dbset) / 2,
      par("usr")[3] - SPO$xLabAdj[2],
      adj = 1,
      labels = "Too many experiments, use tooltip for identification.",
      xpd = T,
      cex = 1,
      col = "red"
    )
    
    ymax <- max(dbset$value)
    fontsize <- 0.5
    
    for (i in 1:num_rowcols[1]) {
      if (modlabel) {
        label <-
          paste(dbset$supplier_obj_id[i], sep = "")
      } else {
        label <- paste(dbset$supplier_obj_id[i], strain, sp, sep = ";")
      }
      mypos <- 3
      if (dbset$value[i] > (ymax * .90)) {
        mypos <- 1
      }
      if (i < num_rowcols[1] * 0.10) {
        mypos <- 4
      }
      if (i > (num_rowcols[1] * .90)) {
        mypos <- 2
      }
      
      if (dbset$outlier[i] < 1 &
          dbset$z_score[i] < z_bio_thr &
          dbset$z_score[i] >= z_tox_thr  & dbset$p_value[i] < z_pval_thr) {
        if (dbset$supplier_obj_id[i] == noquote(cid)) {
          points(i,
                 dbset$value[i],
                 pch = 20,
                 cex = 2,
                 col = "black")
          text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
        } else {
          points(i,
                 dbset$value[i],
                 pch = 20,
                 cex = 2,
                 col = "blue")
          text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
        }
      }
      if (dbset$outlier[i] ||
          (abs(dbset$z_score[i]) >= abs(z_bio_thr) &&
           dbset$p_value[i] > z_pval_thr)) {
        if (dbset$supplier_obj_id[i] == noquote(cid)) {
          points(i,
                 dbset$value[i],
                 pch = 20,
                 cex = 2,
                 col = "black")
          text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
        } else {
          points(i,
                 dbset$value[i],
                 pch = 20,
                 cex = 2,
                 col = "red")
          text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
        }
      }
      if (dbset$outlier[i] == 0 &
          dbset$z_score[i] <= z_tox_thr  & dbset$p_value[i] <= z_pval_thr) {
        if (dbset$supplier_obj_id[i] == noquote(cid)) {
          points(i,
                 dbset$value[i],
                 pch = 20,
                 cex = 2,
                 col = "black")
          text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
        } else {
          points(i,
                 dbset$value[i],
                 pch = 20,
                 cex = 2,
                 col = "yellow")
          text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
        }
      }
      if (dbset$outlier[i] == 0 &
          dbset$z_score[i] >= abs(z_bio_thr)  &
          dbset$p_value[i] <= z_pval_thr) {
        if (dbset$supplier_obj_id[i] == noquote(cid)) {
          points(i,
                 dbset$value[i],
                 pch = 20,
                 cex = 2,
                 col = "black")
          text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
        } else {
          points(i,
                 dbset$value[i],
                 pch = 20,
                 cex = 2,
                 col = "green")
          text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
        }
      }
      if (dbset$supplier_obj_id[i] == noquote(cid)) {
        points(i,
               dbset$value[i],
               pch = 20,
               cex = 2,
               col = "black")
        text(i, dbset$value[i], cex = fontsize, pos = mypos , label)
      }
    }
    
  }


myImagePlot <- function(x, ...) {
  # Plots Heatmap of a screening plate or plate stack
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <- c()
  # check for additional function arguments
  if (length(list(...))) {
    Lst <- list(...)
    if (!is.null(Lst$zlim)) {
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if (!is.null(Lst$yLabels)) {
      yLabels <- c(Lst$yLabels)
    }
    if (!is.null(Lst$xLabels)) {
      xLabels <- c(Lst$xLabels)
    }
    if (!is.null(Lst$title)) {
      title <- Lst$title
    }
  }
  # check for null values
  if (is.null(xLabels)) {
    xLabels <- c(1:ncol(x))
  }
  if (is.null(yLabels)) {
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(
    data = c(1, 2),
    nrow = 2,
    ncol = 1
  ),
  widths = c(1, 1),
  heights = c(4, 1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb(seq(0, 1, length = 256),
                   # Red
                   seq(0, 1, length = 256),
                   # Green
                   seq(1, 0, length = 256))  # Blue
  ColorLevels <- seq(0, 1, length = length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x):1
  yLabels <- yLabels[reverse]
  x <- x[reverse, ]
  
  # Data Map
  # par(mar = c(1,2.5,2.5,2))
  par(mar = c(2, 5, 5, 4))
  image(
    1:length(xLabels),
    1:length(yLabels),
    t(x),
    col = ColorRamp,
    xlab = "",
    ylab = "",
    axes = FALSE,
    zlim = c(min, max)
  )
  if (!is.null(title)) {
    title(main = title)
  }
  axis(
    BELOW <- 1,
    at = 1:length(xLabels),
    labels = xLabels,
    cex.axis = 1.4
  )
  axis(
    LEFT <-
      2,
    at = 1:length(yLabels),
    labels = yLabels,
    las = HORIZONTAL <- 1,
    cex.axis = 1.4
  )
  
  # Color Scale
  #par(mar = c(3,2.5,2.5,2))
  par(mar = c(6, 5, 3, 4))
  
  image(
    ColorLevels,
    1,
    matrix(
      data = ColorLevels,
      nrow = length(ColorLevels),
      ncol = 1
    ),
    col = ColorRamp,
    xlab = "",
    ylab = "",
    yaxt = "n",
    cex.axis = 1.4
  )
  
  layout(1)
}


# runs the hierachical tree to generate different graphics
if (dsn != 5 && dsn != 8) {
  sql_exp <-
    paste(
      "SELECT tox FROM `res_exp` WHERE sp='",
      noquote(sp),
      "' and s='",
      noquote(strain),
      "' and lib ='",
      noquote(lib),
      "'",
      sep = ""
    )
  exprs <- dbSendQuery(mycon, sql_exp)
  exp_tox <- fetch(exprs, n = -1)
  z_tox_thr <- z_tox_thr + exp_tox
}

if (pn) {
  substr <- paste("and plate_number = ", pn, sep = "")
  title <-
    paste0(
      kShowDataColumn ,
      " for ",
      noquote(strain),
      " for ",
      noquote(lib),
      " (",
      noquote(sp),
      ") Plate Nr: ",
      pn
    )
  num_of_br <- 25
  
} else {
  substr <- paste0("")
  title <-
    paste0(kShowDataColumn ,
          " for ",
          noquote(strain),
          " for ",
          noquote(lib),
          " (",
          noquote(sp),
          ")")
  num_of_br <- 55
}

if (dsn == 1) {
  # scatter plot of the plates
  
  sqlstring <-
    paste0(
      "SELECT Z_norm.supplier_obj_id, ",
      noquote(kShowDataColumn ),
      " as value, z_score, outlier, p_value, `value` as norm_value, per_inh FROM `Z_norm` WHERE species='",
      noquote(sp),
      "' and strain='",
      noquote(strain),
      "' ",
      noquote(substr),
      " and library ='",
      noquote(lib),
      "' order by ",
      noquote(kOrderByStatement)
    )
  
  rs <- dbSendQuery(mycon, sqlstring)
  dbset <- fetch(rs, n = -1)
  
  filename <- paste0(noquote(uid_name), ".", kPlotOutputType)
  if (kPlotOutputType == "ps") {
    
    postscript(file = filename, horizontal = TRUE)
    plot_norm_scores(dbset, z_bio_thr, z_tox_thr, z_pval_thr, cid, 1)
    
  } else if (kPlotOutputType == "pdf") {
    
    pdf(file = filename,
        width = 12,
        height = 8)
    plot_norm_scores(dbset, z_bio_thr, z_tox_thr, z_pval_thr, cid, 1)
    
  } else if (kPlotOutputType == "svg") {
    
    devSVG( 
      file = filename,
      width = 12,
      height = 8,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    plot_norm_scores(dbset, z_bio_thr, z_tox_thr, z_pval_thr, cid, 0)
    
  } else if (kPlotOutputType == "screen") {
    
    filename <- paste0(noquote(uid_name), noquote(uid_spec), ".svg")
    par(mar = c(15, 6.8, 4.2, 0.2))
    devSVG(
      file = filename,
      width = 8,
      height = 7,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    plot_norm_scores(dbset, z_bio_thr, z_tox_thr, z_pval_thr, cid, 0)
    
  } else if (kPlotOutputType == "csv") {

    filename <-
      paste0(noquote(uid_name), noquote(uid_spec), ".csv")
    write.table(
      dbset,
      file = filename,
      sep = ",",
      row.names = FALSE,
      qmethod = "double"
    )
  }
  dev.off()
  
  
} else if (dsn == 2) {
  
  sqlstring <-
    paste(
      "SELECT ",
      noquote(kShowDataColumn ),
      " as value, plate_number `plate`, plate_row `row`, plate_column `col`, supplier_obj_id as cid FROM `Z_norm` WHERE species='",
      noquote(sp),
      "' and library ='",
      noquote(lib),
      "' and strain='",
      noquote(strain),
      "' ",
      noquote(substr),
      " order by ",
      noquote(kOrderByStatement),
      sep = ""
    )
  rs <- dbSendQuery(mycon, sqlstring)
  data <- fetch(rs, n = -1) # get datapoints per library
  filename <- paste0(noquote(uid_name), ".", kPlotOutputType)
  
  if (kPlotOutputType == "ps") {
    
    postscript(file = filename, horizontal = TRUE) # 
    par(mfrow = c(2, 1))
    hist(
      data[, 1],
      main = title,
      xlab = "",
      breaks = num_of_br,
      prob = TRUE,
      col = "lightgrey"
    )
    lines(density(data$value, bw = "SJ"), col = 'red') # Use SJ bandwidth, in red
    boxplot(data[, 1], xaxis = FALSE, horizontal = TRUE)
    
  } else if (kPlotOutputType == "pdf") {
    
    pdf(file = filename,
        width = 12,
        height = 8)
    par(mfrow = c(2, 1))
    hist(
      data[, 1],
      main = title,
      xlab = "",
      breaks = num_of_br,
      prob = TRUE,
      col = "lightgrey"
    )
    lines(density(data$value, bw = "SJ"), col = 'red') # Use SJ bandwidth, in red
    boxplot(data[, 1], xaxis = FALSE, horizontal = TRUE)
    
  } else if (kPlotOutputType == "svg") {
    
    devSVG(
      file = filename,
      width = 12,
      height = 8,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    par(mfrow = c(2, 1))
    hist(
      data[, 1],
      main = title,
      xlab = "",
      breaks = num_of_br,
      prob = TRUE,
      col = "lightgrey"
    )
    lines(density(data$value, bw = "SJ"), col = 'red') # Use SJ bandwidth, in red
    par(bty = 'n')
    boxplot(
      data[, 1],
      xaxis = FALSE,
      horizontal = TRUE,
      las = 2,
      col = "sienna"
    )
    
  } else if (kPlotOutputType == "screen") {
    
    filename <- paste(noquote(uid_name), noquote(uid_spec), ".svg", sep = "")
    devSVG(
      file = filename,
      width = 8,
      height = 7,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    par(
      mfrow = c(2, 1),
      oma = c(2, 2, 2, 0) + 0.1,
      mar = c(4, 5, 4, 0) + 0.1,
      xpd = NA
    )
    hist(
      data[, 1],
      main = title,
      xlab = "",
      breaks = num_of_br,
      prob = TRUE,
      col = "lightgrey"
    )
    lines(density(data$value, bw = "SJ"), col = 'red') # Use SJ bandwidth, in red
    par(bty = 'n')
    boxplot(
      data[, 1],
      xaxis = FALSE,
      horizontal = TRUE,
      las = 2,
      col = "sienna"
    )
    
  } else if (kPlotOutputType == "csv") {
    ## Not run:
    ## To write a CSV file for input to Excel one might use
    ##dbset = I("dbset \" quote")
    write.table(
      data,
      file = filename,
      sep = ",",
      row.names = FALSE,
      qmethod = "double"
    )
  }
  dev.off()
  
} else if (dsn == 3) {
  # Do image plot of the RAW data plate
  
  if (pn > 0) {
    sqlstring <-
      paste(
        "select plate_value from `exp_info` where organism='",
        noquote(sp),
        "' and expid='",
        noquote(strain),
        "' and library = '",
        noquote(lib),
        "' and replicate_nr = ",
        replicate ,
        " and plate_number = ",
        pn,
        " order by plate_number, plate_column, plate_row",
        sep = ""
      )
    #sqlstring <- paste("select AVG(plate_value) from `exp_info` where expid='",noquote(s),"' and library = '",noquote(db),"' and replicate_nr = ", rn ," group by plate_column, plate_row order by plate_number,  plate_column, plate_row",sep="")
    #sqlstring <- paste("SELECT value FROM `Z_norm`,`may_info` WHERE strain='",noquote(strain),"' and Z_norm.plate_number = ",pn," and species='",noquote(sp),"' and Z_norm.plate_number=may_info.plate_number and Z_norm.plate_row=may_info.plate_row and Z_norm.plate_column=may_info.plate_column and db='",noquote(lib),"' and Z_norm.supplier_obj_id= may_info.supplier_obj_id order by Z_norm.plate_number, Z_norm.plate_column, Z_norm.plate_row ",sep="")
    rs <- dbSendQuery(mycon, sqlstring)
    data <- fetch(rs, n = -1) # get datapoints per library
    data[, 1] <- as.numeric(data[, 1])
    a <- matrix(data[, 1], 8, 12)
    rows <- c("A", "B", "C", "D", "E", "F", "G", "H")
    
    k <- 1
    l <- 12
    rownames(a) <- rows
    colnames(a) <- seq(k, l)
    
  } else if (pn == 0) {
    
    sqlstring <-
      paste(
        "select AVG(value) from `Z_norm` where species='",
        noquote(sp),
        "' and library = '",
        noquote(lib),
        "' and strain='",
        noquote(strain),
        "' group by plate_column, plate_row order by ",
        noquote(kOrderByStatement),
        sep = ""
      )
    rs <- dbSendQuery(mycon, sqlstring)
    data <- fetch(rs, n = -1) # get datapoints per library
    a <- matrix(data[, 1], 8, 10)
    rows <- c("A", "B", "C", "D", "E", "F", "G", "H")
    
    k <- 2
    l <- 11
    rownames(a) <- rows
    colnames(a) <- seq(k, l)
    
  }
  
  filename <- paste0(noquote(uid_name), ".", kPlotOutputType)
  if (kPlotOutputType == "ps") {
    postscript(file = filename, horizontal = TRUE) #
    myImagePlot(a,
                yLabels = rows,
                xLabels = c(k:l),
                title = title)
  } else if (kPlotOutputType == "pdf") {
    
    pdf(file = filename,
        width = 12,
        height = 8)
    myImagePlot(a,
                yLabels = rows,
                xLabels = c(k:l),
                title = title)
    
  } else if (kPlotOutputType == "svg") {
    
    devSVG(
      file = filename,
      width = 12,
      height = 8,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    myImagePlot(a,
                yLabels = rows,
                xLabels = c(k:l),
                title = title)
    
  } else if (kPlotOutputType == "screen") {
    
    filename <- paste0(noquote(uid_name), noquote(uid_spec), ".svg")
    devSVG(
      file = filename,
      width = 8,
      height = 7,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    myImagePlot(a,
                yLabels = rows,
                xLabels = c(k:l),
                title = title)
    
  } else if (kPlotOutputType == "csv") {
    
    ## Not run:
    ## To write a CSV file for input to Excel one might use
    ##dbset = I("dbset \" quote")
    write.table(
      a,
      file = filename,
      sep = ",",
      row.names = TRUE,
      qmethod = "double"
    )
    
  }
  dev.off()
  
} else if (dsn == 4 && pn == 0) {
  
  dbset <- exp_normalize(lib, strain, sp)
  dbset <- get_zscores(dbset)
  dbset <- undermine_find_outliers(dbset)
  plcount <- dbset$replicates
  
  while ((mod(dbset$replicates, 2) == 0) && (plcount >= 2)) {
    filename <- paste0(noquote(uid_name), "_rpc", plcount, ".", kPlotOutputType)
    
    if (kPlotOutputType == "ps") {
      postscript(file = filename, horizontal = TRUE)
      do_replicate_plot(dbset, plcount, 1, kShowDataColumn )
      
    } else if (kPlotOutputType == "pdf") {
      
      pdf(file = filename,
          width = 12,
          height = 8)
      do_replicate_plot(dbset, plcount, 1, kShowDataColumn )
      
    } else if (kPlotOutputType == "svg") {
      
      devSVG(
        file = filename,
        width = 12,
        height = 8,
        onefile = TRUE,
        xmlHeader = TRUE
      )
      do_replicate_plot(dbset, plcount, 0, showpcol)
      
    } else if (kPlotOutputType == "screen") {
      
      filename <-
        paste(noquote(uid_name),
              "_rpc",
              plcount,
              noquote(uid_spec),
              ".svg",
              sep = "")
      devSVG(
        file = filename,
        width = 8,
        height = 7,
        onefile = TRUE,
        xmlHeader = TRUE
      )
      do_replicate_plot(dbset, plcount, 0, kShowDataColumn )
      
    } else if (kPlotOutputType == "csv") {
      
      ## Not run:
      ## To write a CSV file for input to Excel one might use
      ##dbset = I("dbset \" quote")
      x <-
        data.frame(
          dbset$db_identifier,
          dbset$av,
          dbset$z_score,
          dbset$z_pvalue,
          dbset$rav,
          dbset$oav
        )
      colnames(x)[5:ncol(x)] <-
        c(
          "norm value x",
          "norm value y",
          "Z-score x",
          "Z-score y",
          "p-value x",
          "p-value y",
          "raw x",
          "raw y",
          "non-replicate",
          "var",
          "Z-factor"
        )
      write.table(
        x,
        file = filename,
        sep = ",",
        row.names = FALSE,
        qmethod = "double"
      )
    }
    
    dev.off()
    plcount <- plcount - 2
  }
  
} else if (dsn == 8) {
  
  sqlstring <-
    paste0(
      "SELECT supplier_obj_id, Z_norm.strain as strain, species, Z_norm.",
      noquote(kShowDataColumn ),
      " as value, z_score, p_value, outlier, Z_norm.library as library FROM Z_norm join res_exp  on strain =s and species = sp and library =lib join exp_info c on organism = species and Z_norm.library = c.library and expid = strain and c.replicate_nr=0 and Z_norm.plate_number = c.plate_number and Z_norm.plate_row = Z_norm.plate_row and Z_norm.plate_column = c.plate_column WHERE Z_norm.library = '",
      noquote(lib),
      "' and species = '",
      noquote(sp),
      "' and Z_norm.supplier_obj_id = '",
      noquote(cid),
      "' order by ",
      noquote(kOrderByStatement)
    )
  
  sqlstring <-
    paste0(
      "SELECT supplier_obj_id, Z_norm.strain as strain, species, ",
      noquote(kShowDataColumn ),
      " as value, z_score, p_value, outlier, Z_norm.library as library, res_exp.* FROM Z_norm join res_exp  on strain =s and species = sp and library =lib join exp_info c on c.organism = species and Z_norm.library = c.library and c.expid = strain and c.replicate_nr = 1 and Z_norm.plate_number = c.plate_number and Z_norm.plate_row = c.plate_row and Z_norm.plate_column = c.plate_column WHERE Z_norm.library = '",
      noquote(lib),
      "' and species = '",
      noquote(sp),
      "' and Z_norm.supplier_obj_id = '",
      noquote(cid),
      "' order by ",
      noquote(kOrderByStatement)
    )
  #echo "sqlstring";
  rs <- dbSendQuery(mycon, sqlstring)
  dbset <- fetch(rs, n = -1) # get datapoints per library
  
  filename <- paste(noquote(uid_name), ".", kPlotOutputType, sep = "")
  
  if (kPlotOutputType == "ps") {
    
    postscript(file = filename, horizontal = FALSE)
    plot_similar_compounds2(dbset, cid, 1)
    
  } else if (kPlotOutputType == "pdf") {
    
    pdf(file = filename,
        width = 12,
        height = 8)
    plot_similar_compounds2(dbset, cid, 1)
    
  } else if (kPlotOutputType == "svg") {
    
    devSVG(
      file = filename,
      width = 12,
      height = 8,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    plot_similar_compounds2(dbset, cid, 0)
    
  } else if (kPlotOutputType == "screen") {
    
    filename <-
      paste0(noquote(uid_name), noquote(uid_spec), ".svg")

    devSVGTips(
      file = filename,
      width = 8,
      height = 7,
      onefile = TRUE,
      toolTipMode = 1,
      toolTipFontSize = 8,
      title = "",
      xmlHeader = TRUE
    )
    plot_similar_compounds2(dbset, cid, 0)
    
  } else if (kPlotOutputType == "csv") {
    
    write.table(
      dbset,
      file = filename,
      sep = ",",
      row.names = FALSE,
      qmethod = "double"
    )
  }
  dev.off()
  
  
} else if (dsn == 5) {
  
  sqlstring <-
    paste(
      "SELECT supplier_obj_id, species, strain, ",
      noquote(kShowDataColumn ),
      " as value, z_score, p_value, outlier, exp_info.library, exp_info.plate_number, exp_info.plate_row, exp_info.plate_column, date, temperature, absorbance, tox FROM Z_norm JOIN res_exp ON strain = s and species = sp and Z_norm.library = lib JOIN exp_info on expid = strain and species = organism and Z_norm.library = exp_info.library and  exp_info.replicate_nr= 1 and exp_info.plate_number = Z_norm.plate_number and exp_info.plate_row =Z_norm.plate_row and exp_info.plate_column = Z_norm.plate_column WHERE supplier_obj_id in ( Select p from may_keys where b = '",
      noquote(cid),
      "') or supplier_obj_id = '",
      noquote(cid),
      "' order by supplier_obj_id, species, date",
      sep = ""
    )
  
  rs <- dbSendQuery(mycon, sqlstring)
  dbset <- fetch(rs, n = -1) # get datapoints per library
  
  filename <- paste(noquote(uid_name), ".", kPlotOutputType, sep = "")
  if (kPlotOutputType == "ps") {
    postscript(file = filename, horizontal = FALSE)
    plot_similar_compounds(dbset, cid, 1)
  } else if (kPlotOutputType == "pdf") {
    pdf(file = filename,
        width = 12,
        height = 8)
    plot_similar_compounds(dbset, cid, 1)
  } else if (kPlotOutputType == "svg") {
    devSVG(
      file = filename,
      width = 12,
      height = 8,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    plot_similar_compounds(dbset, cid, 0)
  } else if (kPlotOutputType == "screen") {
    filename <- paste(noquote(uid_name), noquote(uid_spec), ".svg", sep = "")
    #par(mfrow=c(1,1),oma = c(1,1,0,1) + 0.1,mar=c(2,2,2,1)+0.1,xpd=NA,cex.lab = 1.4)
    par(mar = c(15, 6.8, 4.2, 0.2))
    devSVG(
      file = filename,
      width = 8,
      height = 7,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    plot_similar_compounds(dbset, cid, 0)
  } else if (kPlotOutputType == "csv") {
    write.table(
      dbset,
      file = filename,
      sep = ",",
      row.names = FALSE,
      qmethod = "double"
    )
  }
  dev.off()
  
} else if (dsn == 6) {
  # select max replicates
  # var replicate
  # query
  
  sqlstring <-
    paste(
      "SELECT replicate_nr, exp_info.plate_number, exp_info.plate_column, exp_info.plate_row, plate_value, library, expid FROM `exp_info` WHERE organism = '",
      noquote(sp),
      "' and library='",
      noquote(lib),
      "' and expid='",
      noquote(strain),
      "' ",
      noquote(substr),
      " ORDER BY replicate_nr, exp_info.plate_number, exp_info.plate_column, exp_info.plate_row",
      sep = ""
    )
  rs <- dbSendQuery(mycon, sqlstring)
  exp_plot_data <- fetch(rs, n = -1)
  m_size <- dim(exp_plot_data)
  num_rows <- m_size[1] / replicate
  # split into matrix
  
  colvec <- rep("grey", nrow(exp_plot_data))
  colvec[exp_plot_data$plate_column == 1 |
           exp_plot_data$plate_column == 12] <- "brown"
  yx <-
    matrix(as.numeric(exp_plot_data$plate_value),
           num_rows,
           replicate,
           byrow = FALSE) #,dimnames = list(dr$supplier_obj_id, dc$strain))
  
  filename <- paste(noquote(uid_name), ".", kPlotOutputType, sep = "")
  
  if (kPlotOutputType == "ps") {
    
    postscript(file = filename, horizontal = FALSE)
    par(mfrow = c(replicate, 1), mar = c(0, 4, 0, 0))
    
  } else if (kPlotOutputType == "pdf") {
    
    pdf(file = filename,
        width = 12,
        height = 8)
    par(mfrow = c(replicate, 1), mar = c(0, 4, 0, 0))
    
  } else if (kPlotOutputType == "svg") {
    
    devSVG(
      file = filename,
      width = 12,
      height = 8,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    par(mfrow = c(replicate, 1), mar = c(0, 4, 0, 0))
    
  } else if (kPlotOutputType == "screen") {
    
    # plot each column
    filename <-
      paste(noquote(uid_name), noquote(uid_spec), ".svg", sep = "")
    devSVGTips(
      file = filename,
      width = 8,
      height = 7,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    par(mfrow = c(replicate, 1), mar = c(0, 8, 0, 0))
    
  } else if (kPlotOutputType == "csv") {
    write.table(
      exp_plot_data,
      file = filename,
      sep = ",",
      row.names = FALSE,
      qmethod = "double"
    )
  }
  for (i in 1:replicate) {
    
    #repstr <- paste("Experiment",i)
    plot(
      1:num_rows,
      yx[, i],
      xaxt = "n",
      ylab = paste("Raw data", i) ,
      col = colvec,
      pch = 20
    )
  }
  
  dev.off()
  
  # produce output data strain replicate
  #filename <- paste("/tmp/exp_plot_data_",noquote(strain),"_r",replicate,".png",sep="")
  #bitmap(file=filename, type = "png256", height = 6, width=6,res=96)
  #plot(1:mmax[1],exp_plot_data[,2],main=strain, xlab= replicate,ylab= "OD value")
  #dev.off()
  
} else if (dsn == 7 && pn == 0) {
  
  do_betw_replicate_plot2(kPlotOutputType)
  
} else if (dsn == 20 && pn == 0) {
  
  spB <- (Sys.getenv("R_SPECIES_2"))
  strainB <- (Sys.getenv("R_STRAIN_2"))
  
  #strain <- "arp6"
  #strainB <- "cla4"
  #lib <- "Arya"
  #sp <- "S. cerevisiae"

  dbset <- exp_normalize(lib, strain, sp)
  dbset <- get_zscores(dbset)
  dbset <- undermine_find_outliers(dbset)
  
  dbsetB <- exp_normalize(lib, strainB, spB)
  dbsetB <- get_zscores(dbsetB)
  dbsetB <- undermine_find_outliers(dbsetB)
  
  filename <- paste(noquote(uid_name), ".", kPlotOutputType, sep = "")
  if (kPlotOutputType == "ps") {
    postscript(file = filename, horizontal = TRUE)
    do_betw_replicate_plot(dbset, dbsetB, 1)
    
  } else if (kPlotOutputType == "pdf") {
    
    pdf(file = filename,
        width = 12,
        height = 8)
    do_betw_replicate_plot(dbset, dbsetB, 1)
    
  } else if (kPlotOutputType == "svg") {
    
    devSVG(
      file = filename,
      width = 12,
      height = 8,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    do_betw_replicate_plot(dbset, dbsetB, 0)
    
  } else if (kPlotOutputType == "screen") {
    
    filename <- paste(noquote(uid_name), noquote(uid_spec), ".svg", sep = "")
    devSVG(
      file = filename,
      width = 8,
      height = 7,
      onefile = TRUE,
      xmlHeader = TRUE
    )
    do_betw_replicate_plot(dbset, dbsetB, 0)
    
  } else if (kPlotOutputType == "csv") {
    
    x <-
      data.frame(
        dbset$db_identifier,
        dbset$av,
        dbset$z_score,
        dbset$z_pvalue,
        dbset$rav,
        dbset$oav,
        dbsetB$db_identifier,
        dbsetB$av,
        dbsetB$z_score,
        dbsetB$z_pvalue,
        dbsetB$rav,
        dbsetB$oav
      )
    
    write.table(
      x,
      file = filename,
      sep = ",",
      row.names = FALSE,
      qmethod = "double"
    )
  }
  
  dev.off()
}
