# ChemGridDataProcessingFunctions.R
#
# 2006 Jan Wildenhain
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

library(DBI)
library(RMySQL)

kDB_HOST <- Sys.getenv("R_DB_HOST")
kDB_NAME <- Sys.getenv("R_DB_NAME")
kDB_USER <- Sys.getenv("R_DB_USER")
kDB_PWD  <- Sys.getenv("R_DB_PWD")

kOrderBias <- "exp_info.plate_number, exp_info.plate_column, exp_info.plate_row"

kMyDBHandShake <- dbConnect(MySQL(), user = kDB_USER, dbname = kDB_NAME, host = kDB_HOST, password = kDB_PWD)

kZscoreBioactiveThreshold <- -3  # minimum Z Score to be classified as active
kZscoreToxicThreshold <- -20  # Z Score to be considered highly active
kPValueThreshold <- 0.001  # P value cutoff to be classified as active 


# MySQL storage function
distributeControlsToTables <- function(lib, s, sp) {
  # Positive and negative controls are stored in different tables in
  # ChemGRID. Positions where given at upload and used for screen summary
  # statistics
  # Args: (Unique SQL key)
  #   lib: library (compound library string identifier)
  #   s: screen identifier (strain string identfier)
  #   sp: species identifier (NCBI species integer identifier)
  stmt <-
    paste0(
      "INSERT ignore into ctrl_neg (sp, lib, s, pn, pc, pr, r, v)  select organism, library, expid, plate_number, plate_column, plate_row,  replicate_nr, plate_value from exp_info where expid = '",
      s,
      "' and library = '",
      lib,
      "' and organism= '",
      sp,
      "' and ( plate_column in (1, 2, 11) or plate_column = 12 and ( plate_row = 'A' or plate_row = 'C' or plate_row = 'E' or plate_row = 'G' ) )")
  dbGetQuery(kMyDBHandShake, stmt)
  
  stmt <-
    paste0(
      "INSERT ignore INTO ctrl_pos (sp, lib, s, pn, pc, pr, r, v)
      SELECT organism, library, expid, plate_number, plate_column, plate_row, replicate_nr, plate_value
      FROM exp_info WHERE expid = '",
      s,
      "' and library = '",
      lib,
      "' and organism= '",
      sp,
      "'
      AND (
      plate_column =12
      AND (
      plate_row = 'B'
      OR plate_row = 'D'
      OR plate_row = 'F'
      OR plate_row = 'H'
      )
      )")
  dbGetQuery(kMyDBHandShake, stmt)
}

# visualisation
showReplicateScatterPlot <- function(dbset, plcount, modlabel, viewcol) {
    # Plot two screens as scatterplot to find outliers
    # Args:
    #   dbset: Data object (data.frame) with raw and normalised data
    #   plcount: Plate count (index)
    #   modlabel: Set title for screen in plot
    #   view col: Data column to display in plot
    if (viewcol == "z_score") {
      dbset$av <- dbset$z_score
    }
    
    replicateA <- paste0("Replicate ", (plcount - 1), "")
    replicateB <- paste0("Replicate ", plcount, "")
    title <-
      paste0(
        "OD Scores for ",
        noquote(dbset$strain),
        " ",
        noquote(dbset$library),
        " (",
        noquote(dbset$species),
        ")")
    plot(
      dbset$av[, (plcount - 1)],
      dbset$av[, plcount],
      xlab = replicateA,
      ylab = replicateB,
      main = title)
    
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
            (dbset$z_score[i, (plcount - 1)] < kZscoreBioactiveThreshold &
             dbset$z_score[i, (plcount - 1)] >= kZscoreToxicThreshold) &
            (dbset$z_score[i, plcount] < kZscoreBioactiveThreshold &
             dbset$z_score[i, plcount] >= kZscoreToxicThreshold)
          ) &
          (dbset$z_pvalue[i, (plcount - 1)] < kPValueThreshold &
           dbset$z_pvalue[i, plcount] < kPValueThreshold)) {
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
      else if (dbset$z_score[i, (plcount - 1)] < kZscoreToxicThreshold &
               dbset$z_score[i, plcount] < kZscoreToxicThreshold) {
        points(dbset$av[i, (plcount - 1)], dbset$av[i, plcount], pch = 20, col =
                 "yellow")
        text(dbset$av[i, (plcount - 1)],
             dbset$av[i, plcount],
             cex = fontsize,
             pos = mypos ,
             plotlabel)
      }
      else if (dbset$z_score[i, (plcount - 1)] > abs(kZscoreBioactiveThreshold) &
               dbset$z_score[i, plcount] > abs(kZscoreBioactiveThreshold)) {
        points(dbset$av[i, (plcount - 1)], dbset$av[i, plcount], pch = 20, col =
                 "green")
        text(dbset$av[i, (plcount - 1)],
             dbset$av[i, plcount],
             cex = fontsize,
             pos = mypos ,
             plotlabel)
      }
      else if ((dbset$z_score[i, (plcount - 1)] < kZscoreToxicThreshold |
                dbset$z_score[i, plcount] < kZscoreToxicThreshold) &
               dbset$z_score[i, (plcount - 1)] < kZscoreBioactiveThreshold &
               dbset$z_score[i, plcount] < kZscoreBioactiveThreshold) {
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

# visualisation
PlotNormalisedScores <- function(s, pn, lib, sp) {
  # Perform scatter plot
  # Args:
  #   s: strain identifier
  #   pn: plate number identifier
  #   lib: library identifier
  #   sp: species identifier
  #
  title <- paste("OD Scores ", noquote(s), sep = "")
  sqlstring <-
    paste0(
      "SELECT Z_norm.supplier_obj_id, value, z_score, outlier, p_value FROM `Z_norm`,`may_info` WHERE strain='",
      noquote(s),
      "' and Z_norm.plate_number = ",
      pn,
      " and species='",
      noquote(sp),
      "' and Z_norm.plate_number=may_info.plate_number and Z_norm.plate_row=may_info.plate_row and Z_norm.plate_column=may_info.plate_column and db='",
      noquote(lib),
      "' order by Z_norm.plate_number, Z_norm.plate_row, Z_norm.plate_column"
    )
  rs          <- dbSendQuery(kMyDBHandShake, sqlstring)
  dbset       <- fetch(rs, n = -1)
  num_rowcols <- dim(dbset)
  
  plot(1:num_rowcols[1], dbset$value,
    xlab = "Datapoints",
    ylab = "OD Values",
    main = title
  )
  
  for (i in 1:num_rowcols[1]) {
    if (dbset$outlier[i] < 1 &
        dbset$z_score[i] < kZscoreBioactiveThreshold &
        dbset$z_score[i] >= kZscoreToxicThreshold  &
        dbset$p_value[i] < kPValueThreshold) {
      points(i, dbset$value[i], pch = 20, col = "blue")
      text(i,
           dbset$value[i],
           cex = 0.5,
           pos = 1 ,
           dbset$supplier_obj_id[i])
    }
    if (dbset$outlier[i]) {
      points(i, dbset$value[i], pch = 20, col = "red")
      text(i,
           dbset$value[i],
           cex = 0.5,
           pos = 1 ,
           dbset$supplier_obj_id[i])
    }
    if (dbset$z_score[i] < kZscoreToxicThreshold) {
      points(i, dbset$value[i], pch = 20, col = "yellow")
      text(i,
           dbset$value[i],
           cex = 0.5,
           pos = 1 ,
           dbset$supplier_obj_id[i])
    }
    if (dbset$z_score[i] > abs(kZscoreBioactiveThreshold)) {
      points(i, dbset$value[i], pch = 20, col = "green")
      text(i,
           dbset$value[i],
           cex = 0.5,
           pos = 2 ,
           dbset$supplier_obj_id[i])
    }
  }
  
}

# normalisation
GetCorrectionValues <- function(l, s, sp, r, ctrl) {
  # calculate the means, medians for the control
  # Args:
  #   l: library identifier
  #   s: strain identifier
  #   sp: species identifier
  #   r: replicate identifier
  #   ctrl: control type
  # Return:
  #   ministat: data.frame with plate statistic
  ch_plate_md <- c()
  ch_plate_m  <- c()
  ch_plate_sd <- c()
  p_nx        <- c()
  p_nsd       <- c()
  c_dim       <- c()
  per_enh     <- c()
  raw_enh     <- c()
  
  if (ctrl) {
    table <- paste("ctrl_pos")
  } else {
    table <- paste("ctrl_neg")
  }
  gwhere <-
    paste0("sp='",
           noquote(sp),
           "' and lib='",
           noquote(l),
           "' and s = '",
           noquote(s),
           "'")
  
  sqlstring <-
    paste("SELECT max(pn) FROM `",
          noquote(table),
          "` WHERE ",
          noquote(gwhere),
          sep = "")
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  max.rep.num <- fetch(rs, n = -1)
  starter <- 0
  # SELECT MAX(plate_number) FROM `ctrl_XXX` # eg 4, 8 16, 625, ...
  for (i in 1:max.rep.num[1, 1]) {
    sqlstring <-
      paste(
        "SELECT v FROM `",
        noquote(table),
        "` WHERE ",
        noquote(gwhere),
        " and r=",
        r,
        " and pn=",
        i,
        " ORDER BY pc, pr",
        sep = ""
      )
    rs <- dbSendQuery(kMyDBHandShake, sqlstring)
    tmp_wtr1 <- fetch(rs, n = -1)
    check <- dim(tmp_wtr1)
    c_dim[i]       <- check[1]
    ch_plate_md[i] <- median(as.numeric(tmp_wtr1[, 1]))
    ch_plate_m[i]  <- mean(as.numeric(tmp_wtr1[, 1]))
    ch_plate_sd[i] <- sd(as.numeric(tmp_wtr1[, 1]))
    for (j in 1:check[1]) {
      # normalize through the plate median
      per_enh[starter + j] <-
        as.numeric(tmp_wtr1[j, 1]) / ch_plate_md[i]
      raw_enh[starter + j] <-
        as.numeric(tmp_wtr1[j, 1])
    }
    # simple lowess to get rid of plate bias
    #dol <- lowess(1:check[1],per_enh[(starter+1):(starter+check[1])], f=8/10)
    #warning(check[1],"|",starter,"\n")
    #per_enh[(starter+1):(starter+check[1])] <- per_enh[(starter+1):(starter+check[1])]/dol$y
    p_nx[i]  <-
      mean(per_enh[(starter + 1):(starter + check[1])])
    p_nsd[i] <-
      sd(per_enh[(starter + 1):(starter + check[1])])
    
    starter <- starter + check[1]
  }
  ministat <-
    list(
      norm = per_enh,
      raw = raw_enh,
      p_med = ch_plate_md,
      p_x = ch_plate_m,
      p_sd = ch_plate_sd,
      p_nx = p_nx,
      p_nsd = p_nsd,
      cmpds_per_p = c_dim
    )
  
  return(ministat)
}

# normalisation
TestLowessAdjustment <- function(lib, s, sp, pn, r, pv) {
  # Ignore plates with high variation in LOWESS procedure
  # Args:
  #   lib: Library
  #   s  : Strain experiment id 
  #   sp : Species
  #   r  : Replicate
  #   pv : Plate values of drugs in this plate
  # Return:
  #   Boolean pass on performing LOWESS 1 or not 0 

  if (sd(pv) < 0.33) { 
    return(1) 
  } else { 
    return(0) 
  }
}

PlateActivityAdjustment<- function(pv) {
  # Do quantile evaluation if plate is highly active
  # prior to LOWESS normalisation. If it is replace 
  # highly active values with the lower or upper quantile.
  # Args:
  #   pv : Plate values of drugs in this plate
  # Return:
  #   Adjusted plate activity values

  # get the lower and upper quartile and median of the negative controls
    res <- boxplot(pv, plot = FALSE)
    uppercut <- res$stats[5]
    lowercut <- res$stats[1]
    middle <- res$stats[3]
  # replace all acitves with the average for the average of the platevalues 
    pv[pv > uppercut] <- middle
    pv[pv < lowercut] <- middle
    return(pv)
}

get_top_density <- function(pv) {
  # Get the lower and upper quartile and median from raw data
  # and remove strong outliers
  # Args:
  #   pv : Plate values of drugs in this plate
    res <- boxplot(pv, plot = FALSE)
    uppercut <- res$stats[5]
    lowercut <- res$stats[2]
    return(pv[pv < uppercut & pv > lowercut])
}

get_top_density_median<- function(pv) {
# pv:plate values of drugs in this plate
# replace all acitves with the average for the average of the platevalues 
    return(median(get_top_density(pv)))
}

get_top_density_IQR<- function(pv) {
# pv:plate values of drugs in this plate
# replace all acitves with the average for the average of the platevalues 
    return(IQR(get_top_density(pv)))
}

plate_non_activity_background_signal <- function(lib,s,sp,pn,r,pv) {
  # Compare plate data with neg. control information
  # Args:
  # lib :library
  # s   :strain experiment id 
  # sp  :species
  # r   :replicate
  # pv  :plate values of drugs in this plate
  # Returns:
  #    Median from the non active fraction of the data

    gwhere <- paste("sp='",noquote(sp),"' and lib='",noquote(lib),"' and s = '",noquote(s),"'",sep="")
    sqlstr <- paste("SELECT v FROM `ctrl_neg` WHERE ",noquote(gwhere)," and r=",r," and pn=",pn," ORDER BY pc, pr",sep="")
    rs     <- dbSendQuery(kMyDBHandShake, sqlstr)
    tmp    <- fetch(rs, n = -1)
    if ( (median(pv) / median(tmp$v)) < 0.67 ) {
       return(median(tmp$v))  # if too active take the median of the negative control
    } else {
       return(get_top_density_median(pv))  # otherwise use data
    }
}

experiment_non_activity_background_signal <- function(lib,s,sp,pv) {
  # If whole screen too active negative controls are used normalised
  # to 1, otherwise remove outliers and normalise based on median
  # Args: 
  # lib : library
  # s   : strain experiment id 
  # sp  : species
  # pv  : plate values of drugs in this plate
  # Returns:
  #   Normalised data frame
  
    if ( (median(pv)) < 0.67 ) {
       return(1)
    } else {
       return(get_top_density_median(pv))
    }
}

get_norm_values <- function(library,strain, species, replicate, negnormvec, posnormvec) {
  # Calculate normalized drug OD values #
  # Args:
  # library    : library ID
  # strain     : strain experiment ID 
  # replicate  : replicate ID 
  # posnormvec : vector with positive controls 
  # negnormvec : vector with negative controls
  # Returns:
  #   Normalised data frame
  
  c_dim       <- c()
  per_enh     <- c()
  raw_enh     <- c()
  ch_plate_md <- c()
  ch_plate_m  <- c()
  ch_plate_sd <- c()
  p_nsd       <- c()
  sqlstring   <-
    paste0(
      "SELECT  max(plate_number) FROM `exp_info` WHERE organism='",
      noquote(species),
      "' and library='",
      noquote(library),
      "' and expid = '",
      noquote(strain),
      "'")
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  max.rep.num <- fetch(rs, n = -1)
  
  starter <- 0
  for (i in 1:max.rep.num[1, 1]) {
    sqlstring <-
      paste0(
        "SELECT plate_value FROM `exp_info`, `may_info` WHERE organism='",
        noquote(species),
        "' and library='",
        noquote(library),
        "' and expid = '",
        noquote(strain),
        "' and replicate_nr=",
        replicate,
        " and exp_info.plate_number=",
        i,
        " and exp_info.plate_number = may_info.plate_number and exp_info.plate_column = may_info.plate_column and exp_info.plate_row = may_info.plate_row and db = library ORDER BY ",
        noquote(kOrderBias)
      )
    rs <- dbSendQuery(kMyDBHandShake, sqlstring)
    tmp_raw <- fetch(rs, n = -1)
    check <- dim(tmp_raw)
    c_dim[i] <- check[1]
    # normalize towards the plate median
    ch_plate_md[i] <- median(as.numeric(tmp_raw[, 1]))
    ch_plate_m[i]  <- mean(as.numeric(tmp_raw[, 1]))
    ch_plate_sd[i] <- sd(as.numeric(tmp_raw[, 1]))
    for (j in 1:check[1]) {
        # normalize through the plate median
        per_enh[starter + j] <-
           as.numeric(tmp_raw[j, 1]) / plate_non_activity_background_signal(library,
                                                                         strain,
                                                                         species,
                                                                         i,
                                                                         replicate,
                                                                         as.numeric(tmp_raw[, 1]))
        raw_enh[starter + j] <- as.numeric(tmp_raw[j, 1])
    }
    # adjust for plate activity
    if (TestLowessAdjustment(library, strain, species, i, replicate, per_enh[(starter +
                                                                                1):(starter + check[1])])) {
        prelow <- PlateActivityAdjustment(per_enh[(starter + 1):(starter + check[1])])
        # simple lowess to get rid of plate bias
        dol <- lowess(1:check[1], prelow, f = 1 / 3)
        per_enh[(starter + 1):(starter + check[1])] <-
           per_enh[(starter + 1):(starter + check[1])] / dol$y
    }
    p_nsd[i] <-
      sd(per_enh[(starter + 1):(starter + check[1])])
    starter <- starter + check[1]
    
  }
  # create output
  ministat <-
    list(
      norm      = per_enh,
      raw       = raw_enh,
      p_med     = ch_plate_md,
      p_x       = ch_plate_m,
      p_sd      = ch_plate_sd,
      p_nsd     = p_nsd,
      cmpds_per_p = c_dim
    )
  
  return(ministat)
}


do_lowess <- function(data) {
  # Normalisation for plate effects using LOWESS
  # Args:
  #   Plate data in row, column order
  # Returns:
  #   Lowess adjusted data
  
  # for 4 plates f=1 / 10
  dolowessa   <- lowess(1:data$datapoints, data$av[, 1], f = 1 / 6)
  dolowessb   <- lowess(1:data$datapoints, data$av[, 2], f = 1 / 6)
  data$av[, 1] <- data$av[, 1] / dolowessa$y;  # ratio normalisation
  data$av[, 2] <- data$av[, 2] / dolowessb$y;
  return(data)
}


##########################################################################
# calculate sd per replicate and throw all out which do not fit the test #
##########################################################################

library(outliers)

find_nonreplicates <- function(mat, ps, ctrl_x, ctrl_sd) {
  # calculate sd per replicate and mark MAD outliers and Z-factor between neg ctrl and
  # molecule (sd estimation from replicates)
  # Args:
  #   mat: data plate matrix with row and columns
  #   ps: row index
  #   ctrl_x: neg control mean
  #   ctrl_sd: neg control sd
  # Returns: 
  #   Data matrix with 
  quo      <- c()
  z.factor <- c()
       for (i in 1:ps) {
                 quo[i]     <- sd(mat[i, ])
                 z.factor[i] <- 1 - ( (3 * quo[i] + 3 * ctrl_sd) / (abs(mean(mat[i, ]) - ctrl_x)) )
       }
  # do MAD statistics outlier based on the variance between plates
  out      <- scores(quo, type = "mad", prob = 0)
  ministat <- cbind(out, quo, z.factor)
  return(ministat)
}

undermine_find_outliers <- function(base) {
  # keep outliers with high variability BUT always BIOACTIVE OR TOXIC  
  # adjustment to keep highly active highly variable replicates
  # Args:
  #   Plate screen data frame (data structure)
  # Returns:
  #   Adjusted values 
  
   # use 5 standard deviations from mean distribution
   keep<-(apply(abs(base$z_score), 1, min) > abs(5))  
   for(i in 1:base$datapoints) {
      if (keep[i]) {
         base$oav[i, 1] <- 0;
      }
      # dirty work around single experiments (no replicates)
      if (is.na(base$oav[i, 1])) {
         base$oav[i, 1]<- 0
      }
      if (is.na(base$oav[i, 3])) {
        base$oav[i, 3] < -0
      }
   }
   return(base)
}

exp_normalize <- function(library, strain, species) {
  # read all the data from the MySQLdatabase and perform plate data normalisation
  # Args:
  #   library: library ID
  #   strain : strain ID
  #   species: species ID
  # Return:
  #   Data structure with raw and normalised data
  sqlstring <-
    paste(
      "SELECT supplier_obj_id , exp_info . plate_number , exp_info . plate_row , exp_info . plate_column FROM `exp_info` , `may_info` WHERE organism='",
      noquote(species),
      "' and library='",
      noquote(library),
      "' and expid = '",
      noquote(strain),
      "' and replicate_nr = 1 and exp_info . plate_number = may_info . plate_number and exp_info . plate_column = may_info . plate_column and exp_info . plate_row = may_info . plate_row and db = library ORDER BY ",
      noquote(kOrderBias),
      sep = ""
    )
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  db_drug_id <- fetch(rs, n = -1)
  
  sqlstring <-
    paste(
      "SELECT  count(replicate_nr) FROM `exp_info` WHERE organism='",
      noquote(species),
      "' and library='",
      noquote(library),
      "' and expid = '",
      noquote(strain),
      "' GROUP BY replicate_nr",
      sep = ""
    )
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  max.rep.num <- fetch(rs, n = -1)
  blah <- dim(max.rep.num) # get replicate number
  
  sqlstring <-
    paste(
      "SELECT  max(plate_number) FROM `exp_info` WHERE organism='",
      noquote(species),
      "' and library='",
      noquote(library),
      "' and expid = '",
      noquote(strain),
      "'",
      sep = ""
    )
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  max_p_num <- fetch(rs, n = -1) # get maximum plate number
  
  sqlstring <-
    paste(
      "SELECT  count(*) FROM `ctrl_neg` WHERE sp='",
      noquote(species),
      "' and lib='",
      noquote(library),
      "' and s = '",
      noquote(strain),
      "' and r = 1",
      sep = ""
    )
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  max_ctrl_neg <- fetch(rs, n = -1) # get maximum plate number
  
  sqlstring <-
    paste(
      "SELECT  count(*) FROM `ctrl_pos` WHERE sp='",
      noquote(species),
      "' and lib='",
      noquote(library),
      "' and s = '",
      noquote(strain),
      "' and r = 1",
      sep = ""
    )
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  max_ctrl_pos <- fetch(rs, n = -1) # get maximum plate number
  
  sqlstring <-
    paste(
      "SELECT  count(*) as datapoints FROM `exp_info`, `may_info` WHERE organism='",
      noquote(species),
      "' and library='",
      noquote(library),
      "' and expid = '",
      noquote(strain),
      "' and db = library and exp_info . plate_number = may_info . plate_number and exp_info . plate_column = may_info . plate_column and exp_info . plate_row = may_info . plate_row and replicate_nr = 1",
      sep = ""
    )
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  max_dp <- fetch(rs, n = -1) # get datapoints per library
  
  pc <-
    array(dim = c(max_ctrl_pos[1, 1], blah[1]))  # get controls (positives)
  nc <-
    array(dim = c(max_ctrl_neg[1, 1], blah[1]))  # get controls (negative)
  pc_mean <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get median for plate OD
  pc_median <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get median for plate OD
  pc_sd <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get sd for plate OD
  nc_median <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get median for plate OD
  nc_mean <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get median for plate OD
  nc_sd <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get sd for plate OD
  nc_norm_mean <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get median for plate OD
  nc_norm_sd <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get sd for plate OD
  
  z_factor <- array(dim = c(max_p_num[1, 1], blah[1]))
  plate_mean <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get median for plate OD
  plate_median <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get median for plate OD
  plate_sd <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get sd for plate OD
  plate_norm_sd <-
    array(dim = c(max_p_num[1, 1], blah[1]))  # get sd for plate OD
  av <- array(dim = c(max_dp[1, 1], blah[1]))  # normalize plate values
  rav <- array(dim = c(max_dp[1, 1], blah[1]))  # raw plate values
  oav <- array(dim = c(max_dp[1, 1], 3)) # outlier matrix
  
  # go through all replicates
  for (i in 1:blah[1]) {
    # REMEMBER 1 = positive ctrl 0 = negative ctrl
    tmp <- GetCorrectionValues(library, strain, species, i, 1)
    pc_median[, i] <- tmp$p_med
    pc_mean[, i]   <- tmp$p_x
    pc_sd[, i]     <- tmp$p_sd
    pc[, i]        <- tmp$raw #norm/median(tmp$norm)
    tmp            <- GetCorrectionValues(library, strain, species, i, 0)
    nc_median[, i] <- tmp$p_med
    nc_mean[, i]   <- tmp$p_x
    nc[, i]        <- tmp$raw #norm/median(tmp$norm)
    nc_sd[, i]     <- tmp$p_sd
    nc_norm_mean[, i] <- tmp$p_nx
    nc_norm_sd[, i]   <- tmp$p_nsd
    tmp <- get_norm_values(library, strain, species, i, nc[, i], pc[, i])
    plate_median[, i] <- tmp$p_med
    plate_sd[, i]     <- tmp$p_sd
    plate_mean[, i]   <- tmp$p_x
    av[, i]  <-
      tmp$norm / experiment_non_activity_background_signal(library, strain, species, tmp$norm) # normalize between experiments
    rav[, i] <- tmp$raw # raw data
    plate_norm_sd[, i] <- tmp$p_nsd
    # do ZFACTOR
    z_factor[, i] <-
      1 - (((3 * pc_sd[, i]) + (3 * nc_sd[, i])) / (abs(pc_mean[, i] - nc_mean[, i])))
    z_factor <- replace(z_factor, z_factor == 'NaN', 1)
  }
  
  # find non replicates in the datasets
  te <- c()
  ti <- c()
  starter <- 1
  
  for (i in 1:max_p_num[1, 1]) {
    sqlstring <-
      paste(
        "SELECT  count(*) as platepoints FROM `exp_info`, `may_info` WHERE organism='",
        noquote(species),
        "' and library='",
        noquote(library),
        "' and expid = '",
        noquote(strain),
        "' and replicate_nr = 1 and exp_info . plate_number = ",
        i,
        " and db = library and exp_info . plate_number = may_info . plate_number and exp_info . plate_column = may_info . plate_column and exp_info . plate_row = may_info . plate_row",
        sep = ""
      )
    rs <- dbSendQuery(kMyDBHandShake, sqlstring)
    max_p <- fetch(rs, n = -1) # get datapoints per library
    
    te[i] <- starter
    ti[i] <- (starter + max_p[1, 1] - 1)
    oav[te[i]:ti[i], 1:3] <-
      find_nonreplicates(av[te[i]:ti[i], ], max_p[1, 1], mean(nc_norm_mean[i, ]), mean(nc_norm_sd[i, ]))
    
    starter <- starter + max_p[1, 1]
    #plot(testdata[,2])
    
  }
  
  base <-
    list(
      db_identifier = db_drug_id,
      av = av,
      rav = rav,
      oav = oav,
      nc = nc,
      pc = pc,
      replicates = blah[1],
      datapoints = max_dp[1, 1],
      outlier = sum(oav[, 1]),
      strain = strain,
      library = library,
      species = species,
      plate_raw_median = plate_median,
      plate_raw_sd = plate_sd,
      pc_median = pc_median,
      pc_sd = pc_sd,
      nc_median = nc_median,
      nc_sd = nc_sd,
      plates = max_p_num[1, 1],
      neg_controls = max_ctrl_neg,
      pos_controls = max_ctrl_pos,
      z_factor = z_factor,
      pc_mean = pc_mean,
      pc_mean = nc_mean,
      plate_raw_mean = plate_mean,
      plate_norm_sd = plate_norm_sd,
      nc_norm_sd = nc_norm_sd
    )
  
  return(base)
}


get_zscores <- function(base) {
  # do z_score statistics
  # Args:
  #   Takes base data structure (see above)
  # Returns:
  #   Extended base data structure
  exp_results <- base$av  # wtbind, mtbind or wtav, mtav
  exp_ctrl_p  <- base$pc  # get z score cutoff for toxicity
  
  p      <- ncol(exp_results)
  p_ctrl <- ncol(exp_ctrl_p)
  
  #result$estimators <- 'quartiles'
  # quartile-based estimations
  if (mean(base$plate_raw_median / base$nc_median) < 0.7) {
    
    tmpset <- base$nc / median(base$nc)
    m.est  <- apply(tmpset, 2, median)
    iqr    <- apply(tmpset, 2, IQR)
    s.est  <- iqr / (qnorm(0.75) - qnorm(0.25))
    
  } else {
    
    m.est <- apply(exp_results, 2, get_top_density_median)
    iqr   <- apply(exp_results, 2, get_top_density_IQR)
    s.est <- iqr / (qnorm(0.75) - qnorm(0.25))
    
  }
  # normalize data using quartile-based estimators of central tendency and dispersion
  
  z <- exp_results
  P.value <- exp_results
  E.value <- exp_results
  
  for (i in 1:p) {
    # calculate z-score
    # 0.003 is average reader error
    z[, i] <- (z [, i] - m.est[i]) / (s.est[i] + 0.003)  
    # format(z,nsmall=4)
    P.value[, i] <- pnorm(abs(z[, i]), lower.tail = F)
  }
  
  if (mean(exp_ctrl_p) > 1) {
    
    z_ctrl <- 1 / exp_ctrl_p
  
  } else {
    
    z_ctrl <- exp_ctrl_p

  }
  for (i in 1:p) {

        z_ctrl[, i] <- (z_ctrl[, i] - m.est[i]) / (s.est[i] + 0.003)  
    
  }
  
  base$z_ctrl_p <- z_ctrl
  # estimate Z-Score for toxic controls
  base$z_score  <- z
  base$z_pvalue <- P.value
  
  return(base)
}


get_sam_score <- function(baseA, baseB) {
  # do SAM statistics
  # Args:
  #   baseA: data structure for dataset A 
  #   baseB: data structure for dataset B
  # Returns:
  #   baseA: with added stats (used for visualisation A vs. B)
  
  p_value   <- c()
  t_stat    <- c()
  sam       <- c()
  std_error <- c()
  
  cA <- IQR(baseA$av)
  cB <- IQR(baseB$av)
  
  for (i in 1:baseA$datapoints) {
    
    t_test       <- t.test(baseA$av[i, ], baseB$av[i, ])
    p_value[i]   <- t_test$p.value
    t_stat[i]    <-  t_test$statistic
    std_error[i] <- sqrt(var(baseA$av[i, ]) / baseA$replicate + var(baseB$av[i, ] /
                                                                     baseB$replicate))
    sam[i]       <- (mean(baseA$av[i, ]) - mean(baseB$av[i, ])) / (std_error[i] +
                                                              ((cA + cB) / 2))
    
  }
  
  baseA$sam_score     <- sam
  baseA$sam_pvalue    <- p_value
  baseA$sam_std_error <- std_error
  
  return(baseA)
}

platestat2db <- function(base) {
  # DB inserts for plate data statistics
  # Args:
  #  base: data structure with raw and normalised data and statistics

   for (j in 1:base$plates) {
      for (i in 1:base$replicates) {
          sql.str <- paste("INSERT INTO `res_plate` (`sp`, `lib`, `s`, `pn`, `r`, `z`, `v_pn`, `v_ctrl_n`, `v_ctrl_p`, `v_pn_no`, `v_ctrl_n_no`, `m_pn`, `m_ctrl_n`, `m_ctrl_p`) VALUES ('", base$species,"','",base$library,"','",base$strain,"',",j,",",i,",",base$z_factor[j,i],",",base$plate_raw_sd[j,i],",",base$nc_sd[j,i],",",base$pc_sd[j,i],",",base$plate_norm_sd[j,i],",",base$nc_norm_sd[j,i],",",base$plate_raw_median[j,i],",",base$nc_median[j,i],",",base$pc_median[j,i],") ON DUPLICATE KEY UPDATE z=",base$z_factor[j,i]," , v_pn=",base$plate_raw_sd[j,i],", v_ctrl_n=",base$nc_sd[j,i],", v_ctrl_p=",base$pc_sd[j,i],", v_pn_no=",base$plate_norm_sd[j,i],", v_ctrl_n_no=",base$nc_norm_sd[j,i],", m_pn=",base$plate_raw_median[j,i],", m_ctrl_n=",base$nc_median[j,i],", m_ctrl_p=",base$pc_median[j,i],sep="")
          dbGetQuery(kMyDBHandShake, sql.str);
      }
   }
}



expstat2db <- function(base) {
  # void DB inserts into screen statistics table
  # Args:
  #  base: data structure with raw and normalised data and statistics
  #
  plate_z <-
    1 -  ((3 * mean(sd(base$pc))) + (3 * mean(sd(base$nc)))) / (abs((mean(sd(
      base$pc
    ))) - mean(median(base$nc))))
  tox     <- mean(median(base$z_ctrl_p))
  sql.str <-
    paste0(
      "INSERT INTO `res_exp` (`sp`, `lib`, `s`, `z`, `tox`, `v_exp`, `m_exp_ctrl_n`, `v_exp_ctrl_n`, `m_exp_ctrl_p`, `v_exp_ctrl_p`) VALUES ('",
      base$species,
      "','",
      base$library,
      "','",
      base$strain,
      "',",
      plate_z,
      ",",
      tox,
      ",",
      mean(sd(base$rav)),
      ",",
      mean(median(base$nc)),
      ",",
      mean(sd(base$nc)),
      ",",
      mean(median(base$pc)),
      ",",
      mean(sd(base$pc)),
      " ) ON DUPLICATE KEY UPDATE z=",
      plate_z,
      ", tox=",
      tox,
      ", v_exp=",
      mean(sd(base$rav)),
      ", m_exp_ctrl_n=",
      mean(median(base$nc)),
      ", v_exp_ctrl_n=",
      mean(sd(base$nc)),
      ", m_exp_ctrl_p=",
      mean(median(base$pc)),
      ", v_exp_ctrl_p=",
      mean(sd(base$pc))
    )
  
  dbGetQuery(kMyDBHandShake, sql.str)
}



insert_stat2db <- function(base) {
  # DB inserts into normalised data point statistics table
  # Args:
  #  base: data structure with raw and normalised data and statistics
  for (i in 1:base$datapoints) {
    if (is.na(base$oav[i, 1])) {
      base$oav[i, 1] <- 1
    }
    
    sql <-
      paste0(
        "insert ignore into `Z_norm` (supplier_obj_id, plate_number, plate_row, plate_column, strain, species, value, outlier, v, z_score, v_z, p_value, v_p, z_factor, library) values ('",
        base$db_identifier[i, 1],
        "', ",
        base$db_identifier[i, 2],
        ",'",
        base$db_identifier[i, 3],
        "',",
        base$db_identifier[i, 4],
        ",'",
        base$strain,
        "','",
        base$species,
        "',",
        median(base$av[i, ]),
        ",",
        base$oav[i, 1],
        ",",
        base$oav[i, 2],
        ",",
        median(base$z_score[i, ]),
        ",",
        sd(base$z_score[i, ]),
        ",",
        median(base$z_pvalue[i, ]),
        ",",
        sd(base$rav[i, ]),
        ",",
        base$oav[i, 3],
        ", '",
        base$library,
        "')"
      )
    dbGetQuery(kMyDBHandShake, sql)
    
  }
}

insert_sam2db <- function(baseA, baseB) {
  # void DB inserts sam difference write
  # Args:
  #  baseA: data structure A with raw and normalised data and statistics
  #  baseB, data structure B
  for (i in 1:baseA$datapoints) {
    sqlstring <-
      paste0(
        "insert into `sam_info` (supplier_obj_id, plate_number, plate_row, plate_column, strain_a, strain_b, SAM, p_value) values ('",
        baseA$db_identifier[i, 1],
        "', ",
        baseA$db_identifier[i, 2],
        ",'",
        baseA$db_identifier[i, 3],
        "',",
        baseA$db_identifier[i, 4],
        ",'",
        baseA$strain,
        "','",
        baseB$strain,
        "',",
        baseA$sam_score[i],
        ",",
        baseA$sam_pvalue[i],
        ")"
      )
    dbGetQuery(kMyDBHandShake, sqlstring)
    
    
  }
  
}

update_stat2db <- function(base) {
  # DB updates if dataset with keys is already there  
  # Args:
  #   data structure
  for (i in 1:base$datapoints) {
    sqlstring <-
      paste0(
        "update `Z_norm` set z_factor = ",
        base$oav[i, 3],
        ", outlier = ",
        base$oav[i, 1],
        ", v = ",
        base$oav[i, 2],
        ", value = ",
        median(base$av[i, ]) ,
        ", v_z = ",
        sd(base$z_score[i, ]) ,
        ", z_score = ",
        median(base$z_score[i, ]) ,
        ", v_p = ",
        sd(base$rav[i, ]) ,
        ", p_value = ",
        median(base$z_pvalue[i, ]) ,
        " where supplier_obj_id='",
        base$db_identifier[i, 1] ,
        "' and strain='",
        base$strain,
        "' and library='",
        base$library,
        "' and species='",
        base$species,
        "'"
      )
    dbGetQuery(kMyDBHandShake, sqlstring)
  }
}

check_data <- function(lib, s, species) {
  # DB entry check, test if already entered if yes update plate data
  # Args:
  #  lib    : library ID
  #  species: species ID
  #  s      : strain ID
  # Returns:
  #     Boolean true or false
  sqlstring <-
    paste0(
      "select distinct strain from `Z_norm` where species='",
      noquote(species),
      "' and library = '",
      noquote(lib),
      "' and strain='",
      noquote(s),
      "'"
    )
  rs   <- dbSendQuery(kMyDBHandShake, sqlstring)
  data <- fetch(rs, n = -1) # get datapoints per library
  
  bool <- TRUE
  if (nrow(data))
    bool <- FALSE
  
  return(bool)
}

turn_around_plates <- function(sp, lib, s, pn, rn) {
  # in DB turn around plates
  # useful if experimenters accidently turned plates (rotate)
  # Args:
  #   sp <- 4932  # species ID 
  #   pn <- 9  # plate number ID
  #   s <- "Tep1"  # strain ID
  #   lib <- "Lopac"  # library ID
  #   rn <- 2  # replicate ID
  internal   <- paste0("")
  # internal <- paste0("and plate_column = 1 and plate_column = 12")
  # internal <- paste0("and plate_column <> 1 and plate_column <> 12")
  # read data
  sqlstring <-
    paste0(
      "select plate_value from `exp_info` where organism = '",
      noquote(sp),
      "' and library = '",
      noquote(lib),
      "' and expid='",
      noquote(s),
      "' and plate_number = ",
      pn,
      " and replicate_nr = ",
      rn ,
      " ",
      noquote(internal),
      " order by  plate_row ASC, plate_column ASC"
    )
  rs        <- dbSendQuery(kMyDBHandShake, sqlstring)
  data      <- fetch(rs, n = -1) # get datapoints per library
  data[, 1] <- as.numeric(data[, 1])
  dp        <- dim(data)
  # reorder data
  sqlstring <-
    paste0(
      "select plate_row, plate_column from `exp_info` where organism = '",
      noquote(sp),
      "' and library = '",
      noquote(lib),
      "' and expid='",
      noquote(s),
      "' and plate_number = ",
      pn,
      " and replicate_nr = ",
      rn ,
      " ",
      noquote(internal),
      " order by  plate_row DESC, plate_column DESC"
    )
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  reorder <- fetch(rs, n = -1) # get datapoints per library
  # update plate with new order
  for (i in 1:dp[1]) {
    sqlstring <-
      paste0(
        "update `exp_info` set plate_value = ",
        data[i, 1],
        " where organism ='",
        sp ,
        "' and library='",
        lib,
        "' and expid='",
        s,
        "' and plate_number = ",
        pn,
        " and replicate_nr =",
        rn,
        " and plate_row = '",
        reorder[i, 1],
        "' and plate_column = ",
        reorder[i, 2]
      )
    dbGetQuery(kMyDBHandShake, sqlstring)
    
  }
  
}


switch_plates_exp <- function(pnA, pnB, s, rn, db, org) {
  # DB rotate two plates if experimenters accidently 
  # switched plates (wrong order) 
  # Args:
  #   pnA <- 4  # switched plate A ID
  #   pnB <- 3  # plate B ID
  #     s <- "MT2481"  # strain ID
  #    db <- "Lopac"  # library ID
  #    rn <- 1  # replicate ID
  #   org <- 4932 # strain ID

   sqlstring <- paste0("update exp_info set plate_number = 9999 where plate_number = ",pnB," and expid ='",noquote(s),"' and library = '",noquote(db),"' and replicate_nr = ", rn ," and organism='",noquote(org),"'");
   dbGetQuery(kMyDBHandShake, sqlstring);
   sqlstring <- paste0("update exp_info set plate_number = ",pnB," where plate_number = ",pnA," and expid ='",noquote(s),"' and library = '",noquote(db),"' and replicate_nr = ", rn ," and organism='",noquote(org),"'");
   dbGetQuery(kMyDBHandShake, sqlstring);
   sqlstring <- paste0("update exp_info set plate_number = ",pnA," where plate_number = 9999 and expid ='",noquote(s),"' and library = '",noquote(db),"' and replicate_nr = ", rn ," and organism='",noquote(org),"'");
   dbGetQuery(kMyDBHandShake, sqlstring);

}

switch_plates_norm <- function(pnA, pnB, s, db, org) {
  # DB rotate two plates if experimenters accidently 
  # switched plates (wrong order) in normalised stats table
  # Args:
  #   pnA <- 4  # switched plate A ID
  #   pnB <- 3  # plate B ID
  #     s <- "MT2481"  # strain ID
  #    db <- "Lopac"  # library ID
  #   org <- 4932 # strain ID
  
   sqlstring <- paste("update Z_norm set plate_number = 9999 where plate_number = ",pnB," and strain ='",noquote(s),"' and library = '",noquote(db),"' and species='",noquote(org),"'",sep="");
   dbGetQuery(kMyDBHandShake, sqlstring);
   sqlstring <- paste("update Z_norm set plate_number = ",pnB," where plate_number = ",pnA," and strain ='",noquote(s),"' and library = '",noquote(db),"' and species='",noquote(org),"'",sep="");
   dbGetQuery(kMyDBHandShake, sqlstring);
   sqlstring <- paste("update Z_norm set plate_number = ",pnA," where plate_number = 9999 and strain ='",noquote(s),"' and library = '",noquote(db),"' and  species='",noquote(org),"'",sep="");
   dbGetQuery(kMyDBHandShake, sqlstring);
   sqlstring <- paste("ALTER TABLE `Z_norm` DROP PRIMARY KEY");
   dbGetQuery(kMyDBHandShake, sqlstring);
   sqlstring <- paste("update Z_norm, may_info set Z_norm.supplier_obj_id = may_info.supplier_obj_id where may_info.plate_number = ",pnA," and may_info.plate_number = Z_norm.plate_number and may_info.plate_column = Z_norm.plate_column and may_info.plate_row = Z_norm.plate_row and library = '",noquote(db),"' and library = db and strain ='",noquote(s),"' and species='",noquote(org),"'",sep="");
   dbGetQuery(kMyDBHandShake, sqlstring);
   sqlstring <- paste("update Z_norm, may_info set Z_norm.supplier_obj_id = may_info.supplier_obj_id where may_info.plate_number = ",pnB," and may_info.plate_number = Z_norm.plate_number and may_info.plate_column = Z_norm.plate_column and may_info.plate_row = Z_norm.plate_row and library = '",noquote(db),"' and library = db and strain ='",noquote(s),"' and species='",noquote(org),"'",sep="");
   dbGetQuery(kMyDBHandShake, sqlstring);
   sqlstring <- paste("ALTER TABLE `Z_norm` ADD PRIMARY KEY ( `species` , `library` , `strain` , `supplier_obj_id` )");
   dbGetQuery(kMyDBHandShake, sqlstring);

}

myImagePlot <- function(x, ...) {
  # ImagePlot produces heatmap of the plate values
  # Min Args:
  #   x    : data.frame with row and column names
  # List Args : ## if property list contains parameters:
  #   title   : title of plot
  #   xLabels : x axis label
  #   yLabels : y axis label
  #   title   : title of plot
  #   zlim    : min and max to scale heatmap colours e.g. (-1,5)
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
  ColorLevels <- seq(min, max, length = length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x):1
  yLabels <- yLabels[reverse]
  x <- x[reverse, ]
  
  # Data Map
  par(mar = c(1, 2, 2.5, 2))
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
    cex.axis = 0.7
  )
  axis(
    LEFT <-
      2,
    at = 1:length(yLabels),
    labels = yLabels,
    las = HORIZONTAL <- 1,
    cex.axis = 0.7
  )
  
  # Color Scale
  par(mar = c(3, 2, 2.5, 2))
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
    cex.axis = 0.7
  )
  
  layout(1)
}

create_cluster_table <- function(lib, sp) {
  # Hierachical Clustering of Libraries Compounds vs. Strains
  # Args:
  #  lib: library ID
  #  sp: species ID
  kMyDBHandShake <- dbConnect(MySQL(), user = kDB_USER, dbname = kDB_NAME, host = kDB_HOST, password = kDB_PWD)
  
  if (sp == "%") {
    spid <- "all"
  } else {
    spid <- sp
  }

  #  CODE to create specific combinations cluster
  #  subselection <- " join HM_columns on HM_columns.strain = Z_norm.strain and hours between 20 and 50 "
  subselection <- " "
  sqlstring <-
    paste0(
      "select z_score,v_z, Z_norm.strain, supplier_obj_id, outlier, p_value, tox from Z_norm JOIN res_exp ON Z_norm.strain = s and species = sp and library = lib ",
      subselection,
      " where library ='",
      noquote(lib),
      "' and species like '",
      noquote(sp),
      "' order by Z_norm.strain, supplier_obj_id"
    )
  rs <- dbSendQuery(kMyDBHandShake, sqlstring)
  ds <- fetch(rs, n = -1) # get all datapoints
  #a <- matrix(data[,1], 79,11,byrow=FALSE)
  num_ds <- dim(ds)
  
  sqlstring <-
    paste(
      "select supplier_obj_id from may_info where db ='",
      noquote(lib),
      "' order by supplier_obj_id",
      sep = ""
    )
  rs       <- dbSendQuery(kMyDBHandShake, sqlstring)
  dr       <- fetch(rs, n = -1) # get chemicals per library
  num_rows <- dim(dr)
  
  sqlstring <-
    paste(
      "select distinct Z_norm.strain from Z_norm ",
      subselection,
      " where library ='",
      noquote(lib),
      "' and species like '",
      noquote(sp),
      "' order by Z_norm.strain",
      sep = ""
    )
  rs       <- dbSendQuery(kMyDBHandShake, sqlstring)
  dc       <- fetch(rs, n = -1) # get the strains per library
  num_cols <- dim(dc)
  
  # Hypergeometric test to calculate activity enrichment for a compound across strains
  hyperg <- function(x, c1) {
    a <- length(x[x < c1])
    b <- length(x[x > abs(c1)])
    k <- length(x)
    N <- ceiling(3450 * 0.92)
    m <- ceiling(3450 * 0.08)
    return(c(dhyper((a + b),
                    m = m,
                    n = N,
                    k = (k - 5)
    ), a, b, k))
  }
  # calculate variability for each compound, likeliness of beeing active again, mean, median,
  xzd <-
    matrix(
      ds$z_score,
      num_rows[1],
      num_cols[1],
      byrow = FALSE,
      dimnames = list(dr$supplier_obj_id, dc$strain)
    )
  hgeometric <- apply(xzd, 1, hyperg, c1 = kZscoreBioactiveThreshold)
  a <- data.frame(
    "xsd"     = apply(xzd, 1, sd),
    "xmean"   = apply(xzd, 1, mean),
    "xmedian" = apply(xzd, 1, median),
    "xmin"    = apply(xzd, 1, min),
    "xmax"    = apply(xzd, 1, max),
    "hg"      = hgeometric[1],
    "sensitive" = hgeometric[2],
    "resistent" = hgeometric[3],
    "nonactive" = (hgeometric[4] - (hgeometric[2] + hgeometric[3])),
    "screened"  = hgeometric[4],
    "ratio"     = ((hgeometric[2] + hgeometric[3]) / hgeometric[4])
  )
  # compute binary representation and store in DB for fingerprint calculations
  dbGetQuery(kMyDBHandShake, (
    paste0(
      "DROP TABLE IF EXISTS compounds_stats_",
      noquote(lib),
      "_",
      noquote(spid)
    )
  ))
  
  chemical_identifiers <- rownames(a)
  a          <- data.frame(cid = chemical_identifiers, a)
  colattr    <- sapply(a, dbDataType, dbObj = kMyDBHandShake)
  colattr[1] <- "char(11)"
  dbWriteTable(
    kMyDBHandShake,
    paste("compounds_stats_", noquote(lib), "_", noquote(spid), sep = ""),
    field.types = colattr,
    row.names = FALSE,
    a
  )
  
  # set outlier zero
  for (i in 1:num_ds[1]) {
    # take standard deviation into account
    # if (ds$z_score[i] < 0) { ds$z_score[i] <- ds$z_score[i]+ds$v_z[i] }
    # else { ds$z_score[i]<-ds$z_score[i]-ds$v_z[i] }
    # classify all the data based on tox and bioactivity thresholds
    if (ds$outlier[i] == 0 &&
        ds$z_score[i] <= (kZscoreToxicThreshold + ds$tox[i]) &&
        ds$p_value[i] <= kPValueThreshold) {
      ds$z_score[i] <- 12
    } # yellow (highly active)
    else if (ds$outlier[i] == 0 &&
             ds$z_score[i] <= (kZscoreBioactiveThreshold) &&
             ds$p_value[i] <= kPValueThreshold) {
      ds$z_score[i] <-
        round(abs(ds$z_score[i]) * 10 / (abs(ds$tox[i]) - abs(kZscoreToxicThreshold)) /
                2) + 1
    } #  bioactive
    else if (ds$outlier[i] == 0 &&
             ds$z_score[i] >= abs(kZscoreBioactiveThreshold) &&
             ds$p_value[i] <= kPValueThreshold) {
      ds$z_score[i] <- 11
    } # green enhancers
    else if (ds$outlier[i] &&
             abs(ds$z_score[i]) > abs(kZscoreBioactiveThreshold) &&
             ds$p_value[i] <= kPValueThreshold) {
      ds$z_score[i] <- 13
    } # red non-replicate outliers
    else if ((
      ds$z_score[i] > kZscoreBioactiveThreshold &&
      ds$z_score[i] < abs(kZscoreBioactiveThreshold)
    ) ||
    ds$p_value[i] > kPValueThreshold)  {
      ds$z_score[i] <- 0
    } # zeros none active
  }
  
  yx <-
    matrix(
      ds$z_score,
      num_rows[1],
      num_cols[1],
      byrow = FALSE,
      dimnames = list(dr$supplier_obj_id, dc$strain)
    )
  
  # remove all rows where none of the values are above 
  # lower threshold -3 and below upper threshold 3.
  rm(ds)
  rm(dr)
  rm(rs)
  
  xx <- yx[apply(abs(yx) > 0, 1, any), ]
  rm(yx)
  # Cluster binary and average linkage
  d         <- dist(xx, method = "binary") 
  hc        <- hclust(d, "ave") # ave, single
  hc$labels <- rep("", length(hc$labels))
  ord       <- hc$order
  dg        <- as.dendrogram(hc)
  
  # create a dendogram and create a 2D cluster plot genes vs. drugs
  transp.d      <- dist(t(xx), method = "binary")
  hc.col        <- hclust(transp.d, "complete")    # clustering
  hc.col$labels <- rep("", length(hc.col$labels))  # remove the labels
  dendro.col    <- as.dendrogram(hc.col)           # dendrogram
  ord.col       <- hc.col$order                    # order
  
  rnames  <- data.frame(rownames(xx[ord, ord.col]))
  entries <- data.frame(xx[ord, ord.col])
  
  binentries <- xx[ord, ord.col] > 0
  binentries[FALSE] = 0
  binbox <- c()
  for (i in 1:nrow(binentries)) {
    binbox[i] <- paste(as.character(binentries[i, ]), collapse = "")
  }
  # write to DB
  cnames <- colnames(xx[ord, ord.col])
  output <- cbind(rnames, entries)
  binbox <- cbind(rnames, binbox)
  colnames(binbox) <- c("supplier_obj_id", "fp")
  bincolattr    <- sapply(binbox, dbDataType, dbObj = kMyDBHandShake)
  bincolattr[1] <- "char(11)"
  bincolattr[2] <- "char(255)"
  colattr <- sapply(output, dbDataType, dbObj = kMyDBHandShake)
  
  colattr[1] <- "char(11)"
  for (i in 1:num_cols[1]) {
    colattr[i + 1] <- "tinyint(4)"
  }
  
  # Web heatmap; write to DB heatmap table for browser
  kMyDBHandShake <- dbConnect(MySQL(), user = kDB_USER, dbname = kDB_NAME, host = kDB_HOST, password = kDB_PWD)
  
  dbGetQuery(kMyDBHandShake, (
    paste0(
      "DROP TABLE IF EXISTS compounds_",
      noquote(lib),
      "_",
      noquote(spid)
    )
  ))
  dbWriteTable(
    kMyDBHandShake,
    paste("compounds_", noquote(lib), "_", noquote(spid), sep = ""),
    field.types = colattr,
    row.names = FALSE,
    output
  )
  # rownames.xx.ord..ord.col..  ## previously version rownames_xx_ord__ord_col__
  dbGetQuery(kMyDBHandShake, (
    paste0(
      "ALTER TABLE `compounds_",
      noquote(lib),
      "_",
      noquote(spid),
      "` CHANGE `rownames.xx.ord..ord.col..` `supplier_obj_id` CHAR( 11 ) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL "
    )
  ))
  
  # compute binary representation and store in DB for fingerprint calculations
  dbGetQuery(kMyDBHandShake, (
    paste(
      "DROP TABLE IF EXISTS compounds_binary_",
      noquote(lib),
      "_",
      noquote(spid),
      sep = ""
    )
  ))
  dbWriteTable(
    kMyDBHandShake,
    paste("compounds_binary_", noquote(lib), "_", noquote(spid), sep = ""),
    field.types = bincolattr,
    row.names = FALSE,
    binbox
  )
  
}

insert_all_sam2db <- function(l, s) {
  # store sam statistics in database
  # Args:
  #   l <- "Lopac"  # library ID
  #   s <- "ynl224c"  # strain ID
    sql  <- paste0("SELECT distinct expid  FROM `exp_info` where expid <>'",noquote(s),"' and library = '",noquote(l),"'")
    rs   <- dbSendQuery(kMyDBHandShake, sql)
    data <- fetch(rs, n = -1) # get datapoints per library
    tmplopacA <- exp_normalize(l, s)
    for (i in 1:nrow(data)) {
       tmplopacB <- exp_normalize(l, data[i, 1])
       tmplopacB <- get_sam_score(tmplopacB, tmplopacA)
       insert_sam2db(tmplopacB, tmplopacA)
    }
}
