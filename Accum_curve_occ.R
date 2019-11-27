# Laura Jimenez 
# First version: May, 2018

### FUNCTIONS ----------

# ARGUMENTS:
# 'sp.name' -> character chain wiht species name, used for plot legends
# 'output.mod' -> csv file with the following columns: long, lat,
#                 output values from SDM, environmental variables
#                 and as much rows as cells in the area of study
# 'occ.pnts' -> csv file with same colums than 'output.mod' and
#               as much rows as occurrence points
# 'null.mod' -> indicate the distribution used as null model, possible values are
#               "binomial" and "hypergeom".
# 'conlev' -> probability that indicates the confidence level to be drawn around
#             the null model. When 'flag=F', 'conlev' must be equal to 0.
# 'bios' -> vector with two numbers that indicate which environmental variables to plot,
#             if bios==0, the scatterplot pairs will be drawn

# Required packages: maptools, fields, dismo

## Difference between 'accum.occ' and 'accum.occ1': 'output.mod' format is different, and,
## the second one don't produce plots in environmental space, while the first one always produces plots

accum.occ <- function(sp.name,output.mod,occ.pnts,null.mod,conlev=0,bios=0){
  # Sort the points according to the suitability calculated by the ENM model (column 3 from mxnt.pns)
    iord <- order(output.mod[,3],decreasing = T) # highly suitable cells first
    nmxnt <- length(iord)
    mxnt.ord <- output.mod[iord,]
  # Raw outputs for the occurrences
    mxnt.occ <- occ.pnts[,3]
  # Number of occurrences
    nocc <- length(mxnt.occ)
  # Prevalence = number of occurrences / number of cells in G
    preval <- nocc / (nmxnt + nocc)
  # Use vector of output values for the occurrences to classify all the cells in the M region
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
  # Print relevant values
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
  # Before making the plots, we will use shades of gray to identfy the different subareas
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
  # but, before plotting we need to crop the world map using output.mod
  # create the clipping polygon
    mxnt.ext <- bbox(SpatialPoints(output.mod[,1:2]))
  ### Start plotting
  # Plot 1: subregions in geographic space
    x11()
    #if( (ncol(output.mod)-3)==2 ) { par(mfrow=c(1,2)) }
    plot(wrld_simpl,xlim=mxnt.ext[1,],ylim=mxnt.ext[2,],col="wheat1",axes=T,bg="azure2",main="Subregions in Geographical Space")
    # add points with corresponding gray shades
    points(mxnt.ord[,1:2],pch=15,col=ci,cex=0.5)
     # add occurrences
    points(occ.pnts[,1:2],pch=19,col="red")
    # add a color bar to distinguish between low and high suitabilities
    # coordinates for 'text' function need to be changed when using different values for 'bios'
    fields::colorbar.plot(mxnt.ext[1,1]+1,mxnt.ext[2,1],
                          adj.y=0,adj.x = 0,strip=0:nsub,col=ci,strip.length = 0.2,strip.width = 0.05)
    text(mxnt.ext[1,1]+4.5,mxnt.ext[2,1]+2,"Suitability of points in M")
    text(mxnt.ext[1,1]+1.5,mxnt.ext[2,1]-0.5,"High")
    text(mxnt.ext[1,1]+7.5,mxnt.ext[2,1]-0.5,"Low")
    legend("bottomleft",legend=sp.name,pch=19,col="red",bty="n")
  ###
  # Plot 2: subregions in environmental space
    cn <- names(output.mod)
    if(length(bios)==1){
      if((ncol(output.mod)-3)>2){ # more than two environmental variables
        cis <- c(rev(ci),rep(2,nocc)) # grey shades for subregions + red for occurrence points
        nc <- 4:ncol(output.mod)
        # build a matrix with background and presence points, sorted for visualization
        plotm <- rbind(as.matrix(mxnt.ord[nmxnt:1,nc]),as.matrix(occ.pnts[,nc]))
        pch.occ <- ifelse(nocc<=50,19,20)
        pcht <- c(rep(18,length(ci)),rep(pch.occ,nocc))
        mypanel <- function(x,y,col=cis,pch=pcht,bgcolor="steelblue4"){
        # function used to color the background in the panels from function pairs
          ll <- par("usr")
          rect(ll[1], ll[3], ll[2], ll[4], col=bgcolor)
          points(x,y,col=col,pch=pch,bg=bgcolor)
        }
        x11()
        pairs(plotm,panel=mypanel,main="Subregions in Environmental Space")
      }
      else{ # there are only two environmental variables
        x11()
        plot(output.mod[,4],output.mod[,5],type="n",xlab=cn[bios[1]],ylab=cn[bios[2]],
             main="Subregions in Environmental Space")
        u <- par("usr")
        rect(u[1], u[3], u[2], u[4], col = "steelblue4", border = "red")
        points(mxnt.ord[nmxnt:1,4:5],pch=15,col=rev(ci),cex=0.5)
        # add occurrences
        points(occ.pnts[,4:5],pch=19,col="red")
        legend("topleft",legend=c(sp.name,"Occurence points","Suitability of points in M"),pch=c(19,19,15),col=c("white","red","grey"),bty="n")
      }
    } else{ # the user selected only two environmental variables
      x11()
      plot(output.mod[,bios[1]],output.mod[,bios[2]],type="n",xlab=cn[bios[1]],ylab=cn[bios[2]],
           main="Subregions in Environmental Space")
      u <- par("usr")
      rect(u[1], u[3], u[2], u[4], col = "steelblue4", border = "red")
      # add points from low to high suitability, works better for the visualization
      points(mxnt.ord[nmxnt:1,bios],pch=15,col=rev(ci),cex=0.5)
      # add occurrences
      points(occ.pnts[,bios],pch=19,col="red")
      # add a color bar to distinguish between low and high suitabilities
        # coordinates for 'text' function need to be changed when using different values for 'bios'
      fields::colorbar.plot(min(output.mod[,bios[1]])+10,max(output.mod[,bios[2]])-110,
                            adj.y=0,adj.x = 0,strip=0:nsub,col=ci,strip.length = 0.2,strip.width = 0.05)
      text(min(output.mod[,bios[1]])+30,max(output.mod[,bios[2]])-45,"Suitability of points in M")
      text(min(output.mod[,bios[1]])+9,max(output.mod[,bios[2]])-130,"High")
      text(min(output.mod[,bios[1]])+50,max(output.mod[,bios[2]])-130,"Low")
      legend("topleft",legend=sp.name,pch=19,col="red",bty="n")
    }
  ###
  # Plot 3: comparison among counts under random selection hypothesis
    x11()
    plot(mxnt.acc[ntt],mxnt.acc[ntt]*preval,type="b",col="red",xlab="Number of cells", ylab="Occurrences",
        main="Accumulation of occurrences",xlim=c(0,mxnt.acc[length(ntt)]),ylim=c(0,nocc),lwd=2)
  # confidence intervals from hypergeom/binomial distribution
    if(conlev > 0){
      #lines(mxnt.acc[ntt],infs,type="b",col="skyblue3",lwd=2)
      #lines(mxnt.acc[ntt],sups,type="b",col="skyblue3",lwd=2)
      segments(mxnt.acc[ntt],infs,mxnt.acc[ntt],sups,col = "gray")
      points(mxnt.acc[ntt],infs,pch=20,col="grey25")
      points(mxnt.acc[ntt],sups,pch=20,col="grey25")
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
  resul <- cbind(occ.acc,mxnt.acc,round((occ.acc/nocc)*100,2),round((mxnt.acc/nmxnt)*100,2))
  colnames(resul) <- c("No.occurrences","No.cells","%Gained Occ","%Area")
  return(resul)
}
### End
