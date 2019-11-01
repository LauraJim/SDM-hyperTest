# Laura Jimenez 
# First version: May, 2018

### FUNCTIONS ----------

# ARGUMENTS:
# 'sp.name' -> character chain wiht species name, used for plot legends
# 'output.mod' -> raster with raw output from SDM
# 'occ.pnts' -> csv file with 3 columns and as many rows as presences of the species.
#               First column contains the name of the species and columns 2:3 contain the lon,lat coordinates.
# 'null.mod' -> indicate the distribution used as null model, possible values are "binomial" and "hypergeom".
# 'conlev' -> probability that indicates the confidence level to be drawn around the null model.
#             When 'flag=F', 'conlev' must be equal to 0.
# 'flag' -> if FALSE, values are calculated and returned but not plotted

## Difference between 'accum.occ' and 'accum.occ1': 'output.mod' and 'occ.pnts' format is different, and,
## the second one don't produce plots in environmental space.

accum.occ1 <- function(sp.name,output.mod,occ.pnts,null.mod,conlev=0,flag=T){
  # Raster to point conversion of MaxEnt output
    mxnt.pns <- rasterToPoints(output.mod)
  # Sort the points according to the suitability calculated by MaxEnt (column 3 from mxnt.pns)
    iord <- order(mxnt.pns[,3],decreasing = T) # highly suitable cells first
    nmxnt <- length(iord)
    mxnt.ord <- mxnt.pns[iord,]
  # Extract the raw outputs for the occurrences
    mxnt.occ <- na.omit(extract(output.mod,occ.pnts[,2:3]))
    # Number of occurrences
    nocc <- length(mxnt.occ)
    # Prevalence = number of occurrences / number of cells in G
    preval <- nocc / (nmxnt + nocc)
  # Use vector of Maxent values for the occurrences to classify all the cells
    # acording to the partition of the interval (0,1) defined by that vector
    brks <- sort(unique(mxnt.occ))
    if(brks[length(brks)] != 1){ # 1 is not included
      if(brks[1] != 0){ # 0 is not included
        mxnt.cls <- tabulate(cut(mxnt.ord[,3],breaks=c(0,brks,1),labels=F,right=T))
      } else{ # 0 is included
        mxnt.cls <- tabulate(cut(mxnt.ord[,3],breaks=c(brks,1),labels=F,right=T))
      }
    } else{ # 1 is included
      if(brks[1] != 0){ # 0 is not included
        mxnt.cls <- tabulate(cut(mxnt.ord[,3],breaks=c(0,brks),labels=F,right=T))
      } else{ # 0 is included
        mxnt.cls <- tabulate(cut(mxnt.ord[,3],breaks=brks,labels=F,right=T))
      }
    }
    # calculate the accumulated number of cells in the subregions
    mxnt.acc <- c(0,cumsum(rev(mxnt.cls)))
  # Count the number of occurrence points in each subregion
    counts <- vector(mode="numeric")
    occs <- sort(mxnt.occ,decreasing = T)
    dupl <- duplicated(occs)*1 # ==1 if value is duplicated, ==0 otherwise
    for(i in 1:nocc){
      if(dupl[i]==0){ # a value can be replicated two or more times
        counts <- c(counts,1)
      } else {
        nn <- length(counts)
        counts[nn] <- counts[nn] + 1
      }
    }
    # calculate the accumulated number of occurrences in the subregions
    occ.acc <- c(0,cumsum(counts),rep(nocc,length(mxnt.acc)-length(counts)-1))
    # select values that contain important imformation
    last <- sum(occ.acc < nocc)
    if(mxnt.acc[last] == nmxnt){
      ntt <- 1:(last+1)
    } else{
      ntt <- 1:length(mxnt.acc)
    }
  # Print the important values
    print(paste("Number of cells in MaxEnt output:",nmxnt),quote=F)
    print(paste("Number of occurrence points:",nocc),quote=F)
    print(paste("Probability of selecting an occurrence point:",round(preval,4)),quote=F)
  # Calculate values of the intervals around the random line using the null model
    if(conlev > 0){
      conlev1 <- (1 - conlev) / 2
      if(null.mod == "binomial"){
        infs <- qbinom(conlev1,mxnt.acc[ntt],preval)
        sups <- qbinom(conlev1,mxnt.acc[ntt],preval,lower.tail = F)
      }
      if(null.mod == "hypergeom") {
        infs <- qhyper(conlev1,m=nocc,n=nmxnt,k=mxnt.acc[ntt])
        sups <- qhyper(conlev1,m=nocc,n=nmxnt,k=mxnt.acc[ntt],lower.tail = F)
      }
    }
    # Now make all the plots
    if(flag==T){
      # before making the plots, we will use shades of gray to identfy the different subareas
      nsub <- length(mxnt.acc)
      cols <- gray((0:nsub/nsub)) #zero indicates black, and one indicates white
      ci <- vector("numeric")
      for (i in 1:nsub){
        # black indicates highly suitable, and white indicates highly unsuitable
        if(i==nsub){
          c <- rep(cols[nsub],nmxnt-length(ci))
          ci <-c(ci,c)
        } else{
          c <- rep(cols[i],mxnt.acc[i+1]-mxnt.acc[i])
          ci <-c(ci,c)
        }
      }
      # we will use the world map from 'maptools'
      data("wrld_simpl", package="maptools")
      # but, before plotting we need to crop the world map using output.mod create the clipping polygon
      mxnt.ext <- bbox(SpatialPoints(mxnt.pns[,1:2]))
      ###
      # Plot 1: subregions in geographic space
      x11()
      plot(wrld_simpl,xlim=mxnt.ext[1,],ylim=mxnt.ext[2,],col="wheat1",axes=T,bg="azure2",main="Subregions in Geographical Space")
      # add points with corresponding gray shades
      points(mxnt.ord[,1:2],pch=15,col=ci,cex=0.5)
      # add occurrences
      points(occ.pnts[,2:3],pch=19,col="red")
      ###
      # Plot 2: comparison among counts under random selection hypothesis
      x11()
      plot(mxnt.acc[ntt],mxnt.acc[ntt]*preval,type="b",col="red",xlab="Number of cells", ylab="Occurrences",
           main="Accumulation of occurrences",xlim=c(0,mxnt.acc[length(ntt)]),ylim=c(0,nocc),lwd=2)
      # confidence intervals from hypergeom/binomial distribution
      if(conlev > 0){
        #lines(mxnt.acc[ntt],infs,type="b",col="skyblue3",lwd=2)
        #lines(mxnt.acc[ntt],sups,type="b",col="skyblue3",lwd=2)
        segments(mxnt.acc[ntt],infs,mxnt.acc[ntt],sups,col = "gray")
        points(mxnt.acc[ntt],infs,pch=19,col="grey25")
        points(mxnt.acc[ntt],sups,pch=19,col="grey25")
        if(null.mod == "binomial") legmod <- paste("Binomial CI, p =",conlev)
        if(null.mod == "hypergeom") legmod <- paste("Hypergeometric CI, p =",conlev)
      }
      # under non-random selection hypothesis
      lines(mxnt.acc[ntt],occ.acc[ntt],type="s",col="blue",lwd=2)
      if(max(ntt)<=50){
        text(mxnt.acc[ntt],occ.acc[ntt],labels=occ.acc,pos=2)
      } else {
        rind <- seq(1,length(ntt),by=200) #%#
        text(mxnt.acc[rind],occ.acc[rind],labels=occ.acc[rind],pos=2)
      }
      legend("bottomright",legend=c(sp.name,"Random counts",legmod,"SDM counts"),lwd=2,col=c("white","red","skyblue3","blue"),bty="n")
    }
    resul <- cbind(occ.acc,mxnt.acc,round((occ.acc/nocc)*100,2),round((mxnt.acc/nmxnt)*100,2))
    colnames(resul) <- c("No.occurrences","No.cells","%Gained Occ","%Area")
    return(resul)
}
### End
