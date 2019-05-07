#!/bin/env Rscript

## get options from cmd
###############################################################################
args<-commandArgs(TRUE)
if(length(args) < 6) {
  cat(length(args))
  cat("\n")
  cat("   Function: visualize number_of_molecule/total_size_of_molecule against molecule_cov/molecule_len/molecule_readNum.\n\n")
  cat("   USAGE: plotMoleculeCovLen.R mCOV.txt mLEN.txt mPAR.txt mREN.txt mRDI visCovLen \n")
  cat("                               mCOV.txt \t file: molecule_base_cov ~ molcule_num \n")
  cat("                               mLEN.txt \t file: molecule_len      ~ total_molecule_len/molecule_num \n")
  cat("                               mPAR.txt \t file: molecule_num      ~ partition/barcode_num \n")
  cat("                               mREN.txt \t file: read_num          ~ molecule_num \n")  
  cat("                               mRDI.txt \t file: read_distance     ~ case_num \n")  
  cat("                               visCovLen \t output prefix string \n\n")
}else
{
  ## check args
  covfile  = args[[1]]
  lenfile  = args[[2]]
  parfile  = args[[3]]
  readNfile= args[[4]]
  readDfile= args[[5]]
  prefix   = args[[6]]
  if(file.exists(covfile)!=TRUE)
  {
    stop(paste("Cannot find molecule-stat file: ", covfile, ". Please check it! Program exited.\n",sep=""))
  }
  if(file.exists(lenfile)!=TRUE)
  {
    stop(paste("Cannot find molecule-stat file: ", lenfile, ". Please check it! Program exited.\n",sep=""))
  } 
  if(file.exists(parfile)!=TRUE)
  {
    stop(paste("Cannot find molecule-stat file: ", parfile, ". Please check it! Program exited.\n",sep=""))
  }
  if(file.exists(readNfile)!=TRUE)
  {
    stop(paste("Cannot find molecule-stat file: ", readNfile, ". Please check it! Program exited.\n",sep=""))
  }
  if(file.exists(readDfile)!=TRUE)
  {
    stop(paste("Cannot find molecule-stat file: ", readDfile, ". Please check it! Program exited.\n",sep=""))
  }  
  
  if(1)
  {
    # colors <- c("red", "darkcyan", "blue","darkorange","cyan3", "azure3")
    statbase <- basename(covfile)
    pdf(paste(prefix, '_molecule_stat.pdf', sep=""), width = 8, height = 9)
    par(mfrow=c(4,2))
    par(mai = c(0.6, 0.8, 0.5, 0.5)); # margin: bottom, left, top, right
    
    ########################################## A. cov info: two plot ################################################
    stat     <- read.table(covfile)
    if(length(stat$V1) > 0)
    {
      #1.1 molecule base coverage ~ molecule number
      ymax     <- max(stat$V2[1:min(length(stat$V1), 50)])
      xmax     <- which.max(stat$V2[1:min(length(stat$V1), 50)])
      xleft    <- stat$V1[1]
      xrange   <- xmax*2
      if(xrange > length(stat$V1))
      {
        xrange <- length(stat$V1)
      }
      xright <- 1
      if(stat$V1[xmax]>0.5) xright <- 3*stat$V1[xmax]
      plot(stat$V1, 
           stat$V2, 
           col="red",
           cex=1,
           cex.lab = 1.2, 
           cex.axis = 1.2,    
           type = "b",
           xlim=c(0, xright), 
           ylim=c(1, 1.2*ymax), 
           xlab="Molecule base coverage (x)", 
           ylab="Number of molecules")
      abline(v=stat$V1[xmax], col="red", lty="dotted")              
      title(paste("A.1 Molecule base coverage", sep=""))   
      legend("topright", 
             pch=c("*","*","*","*"),
             col="red",
             text.col = "red",
             c(paste("peak    : ", formatC(stat$V1[xmax], width=4, digits = 4), "x",sep=""), 
               paste("mean   : ",  formatC(round(sum(stat$V1[1:xrange]* stat$V2[1:xrange])/sum(stat$V2[1:xrange]), digits = 4), width=4), "x", sep=""), 
               paste("median: ",   formatC(round(median(rep(stat$V1[1:xrange], stat$V2[1:xrange]/100)), digits = 4), width=4), "x", sep=""), 
               paste("max     : ", formatC(max(stat$V1), width=4), "x", sep="")),
             cex = 1,
             box.col="NA")
      #1.2 molecule bases in total ~ molecule base coverage
      total_bases_of_molecules_with_that_cov <- stat$V2*stat$V1/1000;
      ymax <-  max(total_bases_of_molecules_with_that_cov)
      xmax <- which.max(total_bases_of_molecules_with_that_cov[1:min(length(stat$V1), 50)])
      xright <- 1
      if(stat$V1[xmax]>0.5) xright <- 3*stat$V1[xmax]     
      plot(stat$V1, 
           total_bases_of_molecules_with_that_cov, 
           col="red",
           cex=1,
           cex.lab = 1.2, 
           cex.axis = 1.2,    
           type = "b",
           xlim=c(xleft, 1.0), 
           ylim=c(1, 1.2*ymax), 
           xlab="Molecule base coverage (x)",
           ylab="Molecule_base_coverage *\n Number_mole (Kb)")
      abline(v=stat$V1[xmax], col="red", lty="dotted")     
      title(paste("A.2 Total molecule bases", sep=""))
      legend("topright", 
             pch=c("*","*","*","*"),
             col="red",
             text.col = "red",
             c(paste("peak    : ", formatC(stat$V1[xmax], width=4, digits = 4), "x",sep=""), 
               paste("mean   : ",  formatC(round(sum(stat$V1[1:2*xmax]* stat$V2[1:2*xmax])/sum(stat$V2[1:2*xmax]), digits = 4), width=4), "x", sep=""), 
               paste("median: ",   formatC(round(median(rep(stat$V1[1:2*xmax], stat$V2[1:2*xmax]/100)), digits = 4), width=4), "x", sep=""), 
               paste("max     : ", formatC(max(stat$V1), width=4), "x", sep="")),
             cex = 1,
             box.col="NA")
      rm(total_bases_of_molecules_with_that_cov)
    }
    ############################################ B. len info: two plots  #############################################
    statbase <- basename(lenfile)
    stat     <- read.table(lenfile)
    if(length(stat$V1) > 0)
    {
      ymax     <- max(stat$V2[1:min(50, length(stat$V1))])/1000000
      xmax     <- which.max(stat$V2[1:min(50, length(stat$V1))])
      xleft    <- stat$V1[1]
      xright   <- 150
      if(stat$V1[xmax] > 80) xright <- 3*stat$V1[xmax]
      # 2.1 Molcule length (1kb bin) ~ Molecule number
      plot(stat$V1, 
           stat$V3, 
           col="blue",
           cex = 1,
           cex.lab = 1.2, 
           cex.axis = 1.2,    
           type = "b",
           xlim=c(xleft, xright), 
           ylim=c(1, 1.2*max(stat$V3[2:min(xright, length(stat$V3))])), 
           xlab="Molecule length (kb)", 
           ylab="Number of molecules") 
      abline(v=stat$V1[xmax], col="blue", lty="dotted")
      legend("topright", 
             pch="*",
             text.col = "blue",
             col="blue",
             paste("total molecules: ", format(sum(stat$V3), format="d", big.mark=","), "",sep=""),
             cex=1,
             box.col="NA")
      title(paste("B.1 Molcule length (1kb bin)", sep=""))      
      # 2.2 molecule length ~ total molecule length
      plot(stat$V1, 
           stat$V2/1000000000, 
           col="blue",
           cex = 1,
           cex.lab = 1.2, 
           cex.axis = 1.2,
           type = "b",
           xlim=c(xleft, xright), 
           ylim=c(0, 1.2*max(stat$V2)/1000000000), 
           xlab="Molecule length (kb)", 
           ylab="Molecule_length *\n Number_mole (Gb)")
      abline(v=stat$V1[xmax], col="blue", lty="dotted")      
      title(paste("B.2 Total molecule length (1kb bin)", sep=""))   
      len <- min(2*xmax, length(stat$V1))
      legend("topright", 
             pch=c("*","*","*","*"),
             text.col = "blue",
             col = "blue",
             c(paste("peak    : ", formatC(stat$V1[xmax], width=4), " kb",sep=""), 
               paste("mean   : ",  formatC(round(mean(  rep(stat$V1[1:len], stat$V2[1:len]/1000000))), width=4), " kb", sep=""), 
               paste("median: ",   formatC(round(median(rep(stat$V1[1:len], stat$V2[1:len]/1000000))), width=4), " kb", sep=""), 
               paste("max     : ", formatC(max(stat$V1), width=4), " kb", sep="")),
             cex = 1,
             box.col="NA")
    }
    ################################################# 3. read number per molecule  ###################################
    rawstat2     <- read.table(readNfile)
    if(length(rawstat2$V1) > 0)
    {
      # 3.1 read number per molecule ~ number of molecules
      stat     <- rawstat2
      ymax     <- max(stat$V2[1:min(50, length(stat$V1))])
      xmax     <- which.max(stat$V2[1:min(50, length(stat$V1))])
      xleft    <- stat$V1[1]
      xright   <- 100 
      if(stat[xmax,1] > 30) xright <- 5*stat[xmax,1]
      #if(xmax<=3) xright <- 20
      plot(stat[,1], 
           stat[,2], 
           col="darkcyan",
           cex = 1,
           cex.lab = 1.2, 
           cex.axis = 1.2,
           type = "b",
           xlim=c(0, xright), 
           ylim=c(1, 1.2*max(stat[,2])), 
           xlab="Number of reads (of one molecule)", 
           ylab="Number of molecules")
      abline(v=stat[xmax,1], col="darkcyan", lty="dotted")     
      legend("topright", 
             pch=c("*", "*"),
             text.col = "darkcyan",
             col = "darkcyan",
             c(paste("total reads: ", format(sum(stat[,1]*stat[,2]), format="d", big.mark=","), "",sep=""),
               paste("peak    : ", formatC(stat[xmax,1], width=4),sep="")
             ),
             cex=1,
             box.col="NA")   
      title(paste("C.1 Molecule read number", sep=""))
      ################ 3.2 read number per molecule ~ total number of reads of molecules
      cntcol = 2;
      # recover missed cnt, so that we have continuous read number as 1,2,3...
      rawstat <- cbind(1:rawstat2[length(rawstat2[,1]), 1], rep(0, rawstat2[length(rawstat2[,1]), 1]), rep(0, rawstat2[length(rawstat2[,1]), 1]))
      i  <- 0
      while( i < length(rawstat[,1]))
      {
        i               <- i + 1
        rawstat[rawstat2[i, 1], 1] <- rawstat2[i, 1]
        rawstat[rawstat2[i, 1], cntcol] <- rawstat2[i, cntcol] # to change to column 2
      }
      
      # smoothing: combine read number like (1,2), (3,4), (5,6),...
      len = length(rawstat[,1])
      if(length(rawstat[,1])%%2 == 1) len = len - 1;
      y1 <- rawstat[seq(0,len, 2),]
      y2 <- rawstat[seq(1,len, 2),]
      
      stat     <- matrix(0, nrow = len/2, ncol = 2)
      stat[,1] <- round((y1[,1]*y1[,cntcol]+y2[,1]*y2[,cntcol])/(y1[,cntcol]+y2[,cntcol]+0.001), digits = 2)
      stat[,2] <- round(y1[,cntcol]+y2[,cntcol]+0.001)
      
      xmax     <- stat[which.max(stat[3:(len/2),1]*stat[3:(len/2),2])+2, 1]
      xleft    <- stat[1,1]
      xright   <- 100
      if(xmax > 30) xright <- 5*xmax
      #if(xright == 0 || which.max(stat[3:(len/2),1]*stat[3:(len/2),2])+2<=5) xright <- 20
      plot(stat[stat[,1]>0,1], 
           stat[stat[,1]>0,1]*stat[stat[,1]>0,2], 
           col="darkcyan",
           cex = 1,
           cex.lab = 1.2, 
           cex.axis = 1.2,
           type = "b",
           xlim=c(xleft, xright), 
           ylim=c(1, 1.2*max(stat[,1]*stat[,2])), 
           xlab="Number of reads (of one molecule)", 
           ylab="Number_of_reads *\n Number_mole")
      abline(v=xmax, col="darkcyan", lty="dotted")     
      legend("topright", 
             pch=c("*", "*"),
             text.col = "darkcyan",
             col = "darkcyan",
             c(paste("total reads: ", format(sum(stat[,1]*stat[,2]), format="d", big.mark=","), "",sep=""),
               paste("peak    : ", formatC(xmax, width=4),sep="")
             ),
             cex=1,
             box.col="NA")         
    }
    # this is required by DrLin recombis, for automatic read number cutoff setting.
    # note output is total number of reads from all molecules with a certain read number.
    write.table(cbind(stat[stat[,1]>0, 1], round(stat[stat[,1]>0, 1]*stat[stat[,1]>0, 2])), 
                file=paste("smoothed_", readNfile, sep=""),
                col.names = F,
                row.names = F)
    title(paste("C.2 Total molecule read number", sep=""))
    ################################################# 4. read distance distribution   ################################
    stat     <- read.table(readDfile)
    if(length(stat$V1) > 0)
    {
       yyyy    <- abs(stat$V2/1000000*stat$V1)
       ymax    <- max(yyyy[min(200,length(yyyy)):min(5000, length(yyyy))]) # physical distance by adjacent reads
       xmax    <- which.max(yyyy[min(200,length(yyyy)):min(5000, length(yyyy))])+199
       xmaxval <- stat[xmax, 1]
       #cat(paste("xmax=",xmax,"\n"))
       #cat(paste("xmaxval=",xmaxval, "\n"))
       plot(stat$V1, 
            yyyy, 
            col="darkorange",
            cex = 0.3,
            cex.lab = 1.2, 
            cex.axis = 1.2,
            type = "b",
            xlim=c(-5000, 25000), 
            ylim=c(1, 1.2*ymax), 
            xlab="Distances (of R1-R1 and R2-R2; bp) ", 
            ylab="Physical coverage (Mb)")
       abline(v=xmaxval, col="darkorange", lty="dotted") 
       legend("topright", 
              pch="*",
              text.col = "darkorange",
              col = "darkorange",
              paste("peak: ", xmaxval, " bp",sep=""),
              cex=1,
              box.col="NA")
       title(paste("D. Total molecule read distance", sep=""))
    }
    ################################################# 5. molecule per partition info  ################################
    stat     <- read.table(parfile)
    if(length(stat$V1) > 0)
    {
      ymax     <- max(stat$V1*stat$V2)
      xmax     <- which.max(stat$V1*stat$V2)
      xleft    <- stat$V1[1]
      xright   <- min(stat$V1[xmax]*3, length(stat$V1))
      if(xright > length(stat$V1)) xright <- length(stat$V1)
      plot(stat$V1, 
           stat$V1*stat$V2, 
           col="darkgoldenrod",
           cex = 1,
           cex.lab = 1.2, 
           cex.axis = 1.2,
           type = "b",
           xlim=c(xleft, xright), 
           ylim=c(1, 1.2*max(stat$V1*stat$V2)), 
           xlab="Number of molecules", # in one partition
           ylab="No. of molecules *\n No. of partition")
      abline(v=stat$V1[xmax], col="azure3", lty="dotted")     
      legend("topright", 
             pch="*",
             text.col = "darkgoldenrod",
             col = "darkgoldenrod",
             paste("total partitions: ", format(sum(stat$V2), format="d", big.mark=","), "",sep=""),
             cex=1,
             box.col="NA")      
      title(paste("E. Molecule number in partitions", sep=""))
    }    
    # close pdf
    dev.off()
  }
}
