# Laura Jiménez
# May, 2019

### FUNCTIONS ----------
# Description: compare the occurrence-accumulation curves of different SDM/ENMs for a single species

# ARGUMENTS:
# 'mods' -> list with as much elements as models to be compared, each element must be the resulting matrix of values from
#           the function 'accum.occ' applied to the same occurrence data (first two columns only)
# 'nocc' -> number of occurrence points
# 'ncells' -> number of cells in M
# 'sp.name' -> character chain with the name of the species under study
# 'mods.names' -> character vector with the names of the SDMs to be compared, must have the same length as 'mods'
# 'alpha' -> values between 0 and 1 representing the probability of the CI for the null model 

comp.accplot <- function(mods,nocc,ncells,xrange=c(0,ncells),sp.name,mods.names,alpha){
  # number of models to compare
  nmods <- length(mods)
  if(length(mods.names==nmods)){
    # calculate the prevalence of the species in M
    preval <- nocc / ncells
    # Plot the curves of the models to be compared
    x11()
    plot(0:nocc,0:nocc,type="n",xlim=xrange,ylim=c(0,nocc),
         #main=substitute(expr = italic(tl),env = list(tl=sp.name)),
         xlab="Area (number of cells)",ylab="Number of Occurrences")
    # add the curves for each model and calculate the max value in x-axis
    colmod <- colorRampPalette(c("orange","royalblue"))(nmods)
    pt <- 15:(15+nmods)
    for (i in 1:nmods) {
      lines(mods[[i]][,2],mods[[i]][,1],type="o",pch=pt[i],lwd=2,col=colmod[i])
      #xmax <- c(xmax,modelB[[i]][(nrow(modelB[[i]])-1),2])
    }
    # add the line of random classification and confidence bands
    pnts <- floor(seq(0,xrange[2],length=nocc)) #ifelse(xrange[2]==ncells,ncells,max(xmax[2:nmods]))
    a1 <- (1-alpha)/2
    infs <- qhyper(a1,m=nocc,n=ncells-nocc,k=pnts)
    sups <- qhyper(1-a1,m=nocc,n=ncells-nocc,k=pnts)
    lines(0:ncells,0:ncells*preval,col="red",lwd=2)
    for (j in 1:nocc) {
      segments(pnts[j],infs[j],pnts[j],sups[j],col = "gray")
      points(pnts[j],infs[j],pch=20,col="grey25")
      points(pnts[j],sups[j],pch=20,col="grey25")
    }
    # add a legend to identify the different lines
    legend("bottomright",legend=c(mods.names,"Random counts","Hypergeometric-CI"),
    #       pch=c(pt,NA,NA),col=c(colmod,"red","gray"),lwd=3)
  } else{
    print("Warning! 'mods' and 'mods.names' should have the same length")
  }
}

### A second version of the function that allows to cut the x-axis
comp.accplot1 <- function(mods,nocc,ncells,sp.name,mods.names,alpha,xgap){
  # number of models to compare
  nmods <- length(mods)
  if(length(mods.names==nmods)){
    # calculate the prevalence of the species in M
    preval <- nocc / ncells
    # create a sequence of x-values
    pnts <- floor(seq(0,ncells,length=nocc))
    # Plot the curves of the models to be compared
    x11()
    # Start by plotting the curves of each model
    colmod <- colorRampPalette(c("orange","royalblue"))(nmods)
    pt <- 15:(15+nmods)
    # first model
    plotrix::gap.plot(mods[[1]][,2],mods[[1]][,1],type="o",pch=pt[1],col=colmod[1],lwd=2,
                      gap=xgap,gap.axis="x",breakcol="gray",xlim=c(0,ncells),ylim=c(0,nocc),
                      main=substitute(expr = italic(tl),env = list(tl=sp.name)),
                      xlab="Area (number of cells)",ylab="Number of Occurrences")
    # add the curves of the remaining models
    for (i in 2:nmods) {
      plotrix::gap.plot(mods[[i]][,2],mods[[i]][,1],type="o",pch=pt[i],lwd=2,col=colmod[i],
                        gap=xgap,gap.axis="x",breakcol="gray",xlim=c(0,ncells),ylim=c(0,nocc),add=T)
    }
    
    # add the confidence bands
    a1 <- (1-alpha)/2
    infs <- qhyper(a1,m=nocc,n=ncells-nocc,k=pnts)
    sups <- qhyper(1-a1,m=nocc,n=ncells-nocc,k=pnts)
    for (j in 1:nocc) {
      segments(pnts[j],infs[j],pnts[j],sups[j],col = "gray")
      points(pnts[j],infs[j],pch=20,col="grey25")
      points(pnts[j],sups[j],pch=20,col="grey25")
    }
    #plotting the line of random classification
    lines(pnts,pnts*preval,col="red",lwd=2)
    # add a legend to identify the different lines
    legend("bottomright",legend=c(mods.names,"Random counts","Hypergeometric-CI"),col=c(colmod,"red","gray"),lwd=3)
  } else{
    print("Warning! 'mods' and 'mods.names' should have the same length")
  }
}

## END
