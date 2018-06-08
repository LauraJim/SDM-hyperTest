# Jorge Soberón & Laura Jiménez
# May, 2018

### FUNCTIONS ----------

# ARGUMENTS:
# 'out.maxent' -> raster with raw output from MaxEnt
# 'occ.pnts' -> csv file with 3 columns and as many rows as presences of
#  the species. First column contains the name of the species and columns 2:3
#  contain the lon,lat coordinates.
# 'null.mod' -> indicate the distribution used as null model, possible values are
# "binomial" and "hypergeom".
# 'conlev' -> probability that indicates the confidence level to be drawn around
# the null model. When 'flag=F', 'conlev' must be equal to 0.
# 'flag' -> indicates if values must be calculated and plotted,
#  flag=F means that values are caluculated and returned but not plotted

bin.thresholds <- function(out.maxent,occ.pnts,null.mod,conlev=0,flag=T){
  # Raster to point conversion of MaxEnt output
    mxnt.pns <- rasterToPoints(out.maxent)
  # Sort the points according to the suitability calculated by MaxEnt (column 3 from mxnt.pns)
    iord <- order(mxnt.pns[,3],decreasing = T) # highly suitable cells first
    nmxnt <- length(iord)
    mxnt.ord <- mxnt.pns[iord,]
  # Extract the raw outputs for the occurrences
    mxnt.occ <- na.omit(extract(out.maxent,occ.pnts[,2:3]))
    # Number of occurrences
    nocc <- length(mxnt.occ)
    # Prevalence = number of occurrences / number of cells in G
    preval <- nocc / (nmxnt + nocc)
  # Use vector of Maxent values for the occurrences to classify all the cells
    # acording to the partition of the interval (0,1) defined by that vector
    brks <- sort(unique(mxnt.occ))
    mxnt.cls <- tabulate(cut(mxnt.ord[,3],breaks=c(0,brks,1),labels=F,right=T))
    # calculate the accumulated number of cells in the sugregions
    mxnt.acc <- c(0,cumsum(rev(mxnt.cls)))
  # Count the number of occurrence points in each subregion
    counts <- vector(mode="numeric")
    occs <- sort(mxnt.occ,decreasing = T)
    dupl <- duplicated(occs)*1 # ==1 if value is duplicated, ==0 else
    for(i in 1:nocc){
      if(dupl[i]==0){ # a value can be replicated two or more times
        counts <- c(counts,1)
      } else {
        nn <- length(counts)
        counts[nn] <- counts[nn] + 1
      }
    }
    # calculate the accumulated number of occurrences in the sugregions
    occ.acc <- c(0,cumsum(counts),rep(nocc,length(mxnt.acc)-length(counts)-1))
  # Print the important values
    if(nocc <= 200){
      print(paste("Number of cells in MaxEnt output:",nmxnt),quote=F)
      print(paste("Number of occurrence points:",nocc),quote=F)
      print(paste("Probability of selecting an occurrence point:",round(preval,4)),quote=F)
      print("Number of points from MaxEnt in each subregion:",quote=F)
      print(mxnt.acc)
      print("Number of occurrence points in each subregion:",quote=F)
      print(occ.acc)
    }
  # Calculate values of the intervals around the random line using the null model
    if(conlev > 0){
      conlev1 <- (1 - conlev) / 2
      if(null.mod == "binomial"){
        infs <- qbinom(conlev1,mxnt.acc,preval)
        sups <- qbinom(conlev1,mxnt.acc,preval,lower.tail = F)
      } 
      if(null.mod == "hypergeom") {
        infs <- qhyper(conlev1,m=nocc,n=nmxnt,k=mxnt.acc)
        sups <- qhyper(conlev1,m=nocc,n=nmxnt,k=mxnt.acc,lower.tail = F)
      }
    }
    # Now make all the plots
    if(flag==T){
      # Plot 1: subregions
      nsub <- length(mxnt.acc)
      cols <- gray(0:nsub/nsub) #zero indicates black, and one indicates white
      # we use the world map from 'maptools'
      data("wrld_simpl", package="maptools")
      # before plotting we need to crop the map using out.maxent
      # create the clipping polygon
      mxnt.ext <- extent(out.maxent)
      x11()
      plot(wrld_simpl,xlim=mxnt.ext[1:2],ylim=mxnt.ext[3:4],col="wheat1",
           axes=T,bg="azure2",main="Subregions defined by accumulation of occurrences")
      for (i in 1:length(mxnt.acc[-1])) {
        yi <- mxnt.ord[(mxnt.acc[i]+1):mxnt.acc[i+1],1:2]
        # use shades of gray for the different subareas
        # black indicates highly suitable, and white indicates highly unsuitable
        points(yi,pch=15,col=cols[i],cex=0.5)
      }
      # add occurrences
      points(occ.pnts[,2:3],pch=8,col="blue")
      # Plot 2: comparison among counts
      x11()
      # under random selection hypothesis
      plot(mxnt.acc,mxnt.acc*preval,type="b",col="blue",xlab="Number of cells",
           ylab="Occurrences",main="Accumulation of occurrences vs. random classification",
           xlim=c(0,nmxnt),ylim=c(0,nocc),lwd=2)
      # confidence intervals from binomial distribution
      if(conlev > 0){
        lines(mxnt.acc,infs,type="b",col="skyblue3",lwd=2)
        lines(mxnt.acc,sups,type="b",col="skyblue3",lwd=2)
        if(null.mod == "binomial") legmod <- paste("Binomial CI, alpha =",1-conlev)
        if(null.mod == "hypergeom") legmod <- paste("Hypergeometric CI, alpha =",1-conlev)
      }
      # under non-random selection hypothesis
      lines(mxnt.acc,occ.acc,type="b",col="red",lwd=2)
      if(nocc<=50){
        text(mxnt.acc,occ.acc,labels=occ.acc,pos=2)
      } else {
        rind <- seq(1,nocc,by=10)
        text(mxnt.acc[rind],occ.acc[rind],labels=occ.acc[rind],pos=2)
      }
      legend("bottomright",legend=c("Random counts",legmod,"Maxent counts"),lwd=2,col=c("blue","skyblue3","red"),bty="n")
    }
    return(list(mxnt.acc,occ.acc))
}

### MAIN ------------
library(dismo)
# NOTE: packages 'maptools' is also needed

# Read data
  # occurrence points
  lanP <- read.csv("Z:\\Laura\\AUC_Bien_Hecha\\DataPoints\\Lantana20.csv",header=T)
  garP <- read.csv("Z:\\Laura\\AUC_Bien_Hecha\\DataPoints\\P_garamas.csv",header=T)
  nimP <- read.csv("Z:\\Laura\\AUC_Bien_Hecha\\DataPoints\\Catasticta_nimbice.csv",header=T)
  head(nimP)

  # MaxEnt output with small area
  lanS <- raster("Z:\\Laura\\AUC_Bien_Hecha\\MaxentSmall\\Lantana.asc")
  # Another species
  garS <- raster("Z:\\Laura\\AUC_Bien_Hecha\\MaxentSmallPgaramas\\Pgaramas.asc")
  nimS0.5 <- raster("Z:\\Laura\\AUC_Bien_Hecha\\MaxentCatastictaSmall_Beta.5\\Catasticta_nimbice.asc")
  nimS1 <- raster("Z:\\Laura\\AUC_Bien_Hecha\\MaxentCatastictaSmall_Beta1\\Catasticta_nimbice.asc")
  nimB1 <- raster("Z:\\Laura\\AUC_Bien_Hecha\\MaxentCatastictaBig_Beta1\\Catasticta_nimbice.asc")
  nimS2 <- raster("Z:\\Laura\\AUC_Bien_Hecha\\MaxentCatastictaSmall_Beta2\\Catasticta_nimbice.asc")
  x11()
  plot(nimS1)
  points(nimP[,2:3])

    # MaxEnt output with big area
  lanB <- raster("Z:\\Laura\\AUC_Bien_Hecha\\MaxentBig\\Lantana.asc")
  garB <- raster("Z:\\Laura\\AUC_Bien_Hecha\\MaxentBigPgaramas\\Pgaramas.asc")
  x11()
  plot(lanB)

  # Small area
lanS.counts <- bin.thresholds (out.maxent=lanS,occ.pnts=lanP,null.mod="binomial",conlev=0.95,flag=T)

# Big area
lanB.counts <-bin.thresholds (out.maxent=lanB,occ.pnts=lanP,null.mod="binomial",conlev=0.95,flag=T)

# Small area p garamas with more thresholds
garS.counts <-bin.thresholds (out.maxent=garS,occ.pnts=garP,null.mod="binomial",conlev=0.95,flag=T)

# Big area p garamas with more thresholds
garB.counts <-bin.thresholds (out.maxent=garB,occ.pnts=garP,null.mod="binomial",conlev=0.95,flag=T)

# Catastica_nimbice, small area
nimSBeta0.5 <-bin.thresholds (out.maxent=nimS0.5,occ.pnts=nimP,null.mod="binomial",conlev=0.95,flag=T)
nimSBeta1 <-bin.thresholds (out.maxent=nimS1,occ.pnts=nimP,null.mod="binomial",conlev=0.95,flag=T)
nimBBeta1 <-bin.thresholds (out.maxent=nimB1,occ.pnts=nimP,null.mod="binomial",conlev=0.95,flag=T)
nimSBeta2 <-bin.thresholds (out.maxent=nimS2,occ.pnts=nimP,null.mod="binomial",conlev=0.95,flag=T)

# k<- 0:40#lanS.counts[[1]][4]
# plot(k,dbinom(k,garS.counts[[1]][82],20/7486))
# inf <- qbinom(0.025,garS.counts[[1]][82],20/7486)
# sup <- qbinom(0.025,garS.counts[[1]][82],20/7486,lower.tail = F)
# sum(dbinom(inf:sup,garS.counts[[1]][82],20/7486))

x<- 0:80
plot(x,dhyper(x,192,7486,garS.counts[[1]][80]))
