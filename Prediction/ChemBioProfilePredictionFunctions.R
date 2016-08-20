# ChemBioProfilePredictionFunctions.R
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
# Functions used to process data from Pipeline Pilot Multi Class Naive Bayes
# likelihood values into a csv table for further processing
# 
#
# Functions
#
performanceStats <- function(predAct) {
  # Calculate prediciton performance statistics
  #
  # Args:
  #    predAct is a two column dataframe of prediction, actual outcomes
  # Returns:
  #    matrix with prediction summary statistics
  preds = predAct[, 1]
  trues = predAct[, 2]
  xTab <- table(preds, trues)
  clss <- as.character(sort(unique(preds)))
  r <- matrix(NA, ncol = 7, nrow = 1, 
              dimnames = list(c(), c('Acc',
                                    paste("P", clss[1], sep='_'), 
                                    paste("R", clss[1], sep='_'), 
                                    paste("F", clss[1], sep='_'), 
                                    paste("P", clss[2], sep='_'), 
                                    paste("R", clss[2], sep='_'), 
                                    paste("F", clss[2], sep='_'))))
  r[1, 1] <- sum(xTab[1, 1], xTab[2, 2]) / sum(xTab)  # Accuracy
  r[1, 2] <- xTab[1, 1] / sum(xTab[, 1])  # Miss Precision
  r[1, 3] <- xTab[1, 1]/sum(xTab[1,])  # Miss Recall
  r[1, 4] <- (2 * r[1, 2] * r[1, 3]) / sum(r[1, 2], r[1, 3])  # Miss F
  r[1, 5] <- xTab[2, 2] / sum(xTab[, 2])  # Hit Precision
  r[1, 6] <- xTab[2, 2] / sum(xTab[2, ])  # Hit Recall
  r[1, 7] <- (2 * r[1, 5] * r[1, 6]) / sum(r[1, 5], r[1, 6])  # Hit F
  return(r)
}

createCompoundPredictionList <- function(minidat, id, lstart, lend, vstart, vend, rstr1, rstr2) {
  # Parse single entry from data matrix
  #
  # Args:
  #   minidat: PP csv export
  #   lstart: label start index
  #   lend: label end index
  #   vstart: value start index
  #   vend: value end index
  #   rstr1: string expression 1 to detect gene symbol tag
  #   rstr2: string expression 2 to detect value tag
  #
  # Returns:
  #  data.frame with symbol and corresponding value
  #
  idx1 <- lend - lstart
  idx2 <- vend - vstart
  data <- data.frame(sym = rep("none", idx1), score = rep(0, idx2))
  data$sym <- as.character(data$sym) 
  
  for (i in 1:idx1) { 
    data$sym[i]   <- as.character(gsub(rstr1, "", gsub(rstr2, "", minidat[id, (lstart + i )])))
    data$score[i] <- gsub(rstr1, "", gsub(rstr2, "", minidat[id, (vstart + i)]))
  }
  
  return(data)
}

convertPP2RListObj <- function(dat, dat.uID, prefixstr, suffixstr, db.table.name, my.con) {
  # Parse PP csv and build condensed csv
  #
  # Args:
  #   dat: PP csv export
  #   dat.uID: unique data id
  #   prefixstr: string expression 1 to detect gene symbol tag
  #   suffixstr: string expression 2 to detect value tag
  #   db.table.name: store results in MySQL data table name
  #   my.con: database connection (RMySQL)
  #
  # Returns:
  #   is void with full data stored in MySQL in db.table.name
  
  rstr1 <- prefixstr  # prefix 
  rstr2 <- suffixstr  # suffix
  
  tags.getAllNames <- names(dat)
  if(length(i <- grep(".AllNames", tags.getAllNames)))
    i
  min(i)
  max(i)
  max(i) - min(i)
  
  lstart <- min(i) # label starts
  lend <- max(i) # label end
  
  tags.getAllValues <- names(dat)
  if(length(i <- grep(".AllValues", tags.getAllValues)))
    i
  min(i)
  max(i)
  max(i) - min(i)
  
  vstart <- min(i) # values start
  vend <- max(i) # values end #
  
  # verification
  head(dat)[lstart - 2:lstart + 2]
  head(dat)[vstart - 2:vstart + 2]
  head(dat)[vend - 2:vend]
  
  namevector   <- c()
  cid_names    <- c()
  compoundlist <- list()
  
  i <- 1
  datid <- i

  namevector[i] <- as.character(dat[datid, dat.uID])
  cid_names[i] <- as.character(dat[datid, dat.uID])     
  print(paste(i, ":", cid_names[i]))

  compoundlist[[i]] <- createCompoundPredictionList(dat, datid, lstart, lend, vstart, vend, rstr1, rstr2)
  print(head(compoundlist[[i]]))

  tmp.data <- data.frame(compoundlist[[i]])
  tmp.cid  <- rep(cid_names[i], nrow(tmp.data))
  db.table <- data.frame(cid=tmp.cid, sym=tmp.data$sym, lh=tmp.data$score)
  dbWriteTable(my.con, db.table.name, db.table, overwrite = F, append = T, row.names = F)
  

  for (i in 2:nrow(dat) ) {
    
    datid <- i
    
    if (length(datid)) {
      namevector[i] <- as.character(dat[datid, dat.uID])
      cid_names[i]  <- as.character(dat[datid, dat.uID])
      
      print(paste(i, ":", cid_names[i]))
      compoundlist[[i]] <- createCompoundPredictionList(dat, datid, lstart, lend, vstart, vend, rstr1, rstr2)
      print(head(compoundlist[[i]]))

      tmp.data <- data.frame(compoundlist[[i]])
      tmp.cid  <- rep(cid_names[i], nrow(tmp.data))
      db.table <- data.frame(cid = tmp.cid, sym = tmp.data$sym, lh = tmp.data$score)
      dbWriteTable(mycon, db.table.name, db.table, overwrite= F, append = T, row.names = F)
 
    }
  }
}

PP2RListObj2datamatrix <- function(testList) {
  # Parse PP list object and return it as a data matrix
  #
  # Args:
  #   # testList list object with gene symbol and NB likelihoods
  # Returns:
  #   # transformed random forest data matrix
  
  imax <- nrow(testList[[1]][[1]])
  jmax <- length(testList[[1]])
  
  my.gene.pos.vector <- list()
  my.gene.pos <- 1
  
  my.data.matrix <- matrix(0, nrow = imax, ncol = jmax)
  dim(my.data.matrix)
  colnames(my.data.matrix) <- testList[[2]]
  
  for (j in 1:jmax) {
    for (i in 1:imax) {
  #    for (k in 1:length(testList[[1]][[j]]$score)) {
        if (is.null(my.gene.pos.vector[[testList[[1]][[j]]$sym[i]]])) { 
        
          cat(paste("new:", testList[[1]][[j]]$sym[i], "\n")) 
          my.gene.pos.vector[[testList[[1]][[j]]$sym[i]]] <- my.gene.pos
          my.data.matrix[my.gene.pos, j] <- testList[[1]][[j]]$score[i]
          my.gene.pos <- my.gene.pos + 1
        
        } else {
          
          pos.tmp <- my.gene.pos.vector[[testList[[1]][[j]]$sym[i]]]
          cat(paste("pos:", testList[[1]][[j]]$sym[i], ":", pos.tmp, "\n")) 
          my.data.matrix[pos.tmp, j] <- testList[[1]][[j]]$score[i]
          
        }
  #    }
      #takeOrder <- order(testList[[1]][[j]]$sym)
      #resortScore <- testList[[1]][[j]]$score[takeOrder]
      #my.data.matrix[i, j] <- resortScore[i]
    }
  }
  
  dim(my.data.matrix)
  head(my.data.matrix)
  t(my.data.matrix)
  head(t(my.data.matrix))
  
  as.data.frame(my.gene.pos.vector)[1, ]
  rownames(my.data.matrix) <- names(as.data.frame(my.gene.pos.vector))
  rownames(my.data.matrix)
  dim(my.data.matrix)
  head(my.data.matrix)
  t(my.data.matrix)
  head(t(my.data.matrix))
  
  return(my.data.matrix)
}


adjustedRandomForestDataMatrix <- function(d,idat.matrix.adj) {
  # Transform data for random forest analysis
  #
  # Args:
  #   d: data
  #   idat.matrix.adj: adjusted data matrix
  #   prefixstr: string expression 1 to detect gene symbol tag
  dat.matrix.adj <- t(idat.matrix.adj)
  get.strains <- toupper(colnames(dat.matrix.adj))
  get.strains
  which(get.strains == "REF2")  # println test
  
  colnames(d)[15]
  get.transition <- c()
  
  for (j in 2:length(colnames(d))) {
    if(length(i <- grep(colnames(d)[j], get.strains)))
      get.transition[j] <- i
  }
  get.transition
  max.dat4rf  <- ncol(d)
  max.datcmpds <- nrow(dat.matrix.adj)
  dat4rf       <- matrix(0, nrow = max.dat4rf, ncol = max.datcmpds)
  dim(dat4rf)
  
  length(dat.matrix.adj[,get.transition[i]])
  for (i in 2:length(get.transition)) {
    dat4rf[i, ] <- (dat.matrix.adj[, get.transition[i]])
  }
  
  get.transcolnames   <- toupper(colnames(dat.matrix.adj)[get.transition])
  get.transcolnames
  rownames(dat4rf)    <- get.transcolnames
  rownames(dat4rf)[1] <- "cid"
  get.transrownames   <- rownames(dat.matrix.adj)
  get.transrownames
  colnames(dat4rf)    <- get.transrownames
  dat4rf[1, ]          <- get.transrownames
  dat4rf[1:4, ]
  
  return(t(dat4rf)) 
}

