rectbands <- function(x,waveseq,ymin,ymax,width=0.5,plotcol,alfa){
  #This function adds highlighted, semi-transparent rectangular sections to plots
  #This is useful for highlighting regions/bands of a spectrum. The function assumes that 
  #there is a matrix 'X' with x-axis scale vector 'waveseq' that is plotted:
  #matplot(waveseq,t(X),type="l",...) 
  #inputs:
  #x: vector with indices of variables (columns of X) to highlight
  #waveseq: scale for x-axis
  #ymin,ymax: minimum and maximum bounds for rectangle(s)
  #width: Isolated variables are padded with a small amount of space to the left and right to make them more visible when plotted 
  #example for width -> if x=[1 2 4 6 7], the indices 1&2 and 6&7 are both considered groups, whereas 4 is isolated.
  #In this case, the rectangular regions will be at waveseq[1]:waveseq[2], (waveseq[4]-width):(waveseq[4]+width), and waveseq[6]:waveseq[7]
  #plotcol: color of rectangular region(s)
  #alfa: color transparency level in [0,1]. 0=fully transparent, 1=no transparency
  #R code by Cannon Giglio, 2016-2018, Brown Research Group, University of Delaware
  require(scales)
  xseq <- c(1:length(waveseq))
  xlen <- length(x)
  if(xlen<=2){
    if(xlen==1){
      rect((waveseq[x[1]]-width),ymin,(waveseq[x[1]]+width),ymax,border=scales:::alpha(plotcol,alfa),col=scales:::alpha(plotcol,alfa))
    }
    else{
      if((x[1]+1)==x[2]){
        rect((waveseq[x[1]]-width),ymin,(waveseq[x[2]]+width),ymax,border=scales:::alpha(plotcol,alfa),col=scales:::alpha(plotcol,alfa))
      }
      else{
        rect((waveseq[x[1]]-width),ymin,(waveseq[x[1]]+width),ymax,border=scales:::alpha(plotcol,alfa),col=scales:::alpha(plotcol,alfa))
        rect((waveseq[x[2]]-width),ymin,(waveseq[x[2]]+width),ymax,border=scales:::alpha(plotcol,alfa),col=scales:::alpha(plotcol,alfa))
      }
    }
  }
  else{
    xl <- rep(0, length(x))
    lj <- 1
    if((x[1]+1)==x[2]){
      xl[lj] <- x[1]
      lj <- lj+1
    }
    for(i in 2:(xlen-1)){
      if(x[i+1]==(x[i]+1) && x[i-1]!=(x[i]-1)){
        xl[lj] <- x[i]
        lj <- lj+1
      }
    }
    xr <- rep(0,length(x))
    rj <- 1
    for(i in 2:(xlen-1)){
      if(x[i-1]==(x[i]-1) && x[i+1]!=(x[i]+1)){
        xr[rj] <- x[i]
        rj <- rj + 1
      }
    }
    if(x[xlen-1]==(x[xlen]-1)){xr[rj] <- x[xlen]}
    xsolo <- rep(0,length(x))
    xs <- 1
    if((x[1]+1)!=x[2] &&xl[1]!=2){
      xsolo[xs] <- x[1]
      xs <- xs+1
    }
    for(i in 2:(xlen-1)){
      if(x[i-1]!=(x[i]-1) && x[i+1]!=(x[i]+1)){
        xsolo[xs] <- x[i]
        xs <- xs+1
      }
    }
    if((x[xlen]-1)!=x[xlen-1]){xsolo[xs] <- x[xlen]}
    xl <- waveseq[xl[which(xl!=0)]]
    xr <- waveseq[xr[which(xr!=0)]]
    xsolo <- waveseq[xsolo[which(xsolo!=0)]]
    xvec.l <- sort(c(xl,xsolo-width))
    xvec.r <- sort(c(xr,xsolo+width))
    wavemat <- rbind(xvec.l,xvec.r)
    for(i in 1:length(xvec.l)){
      rect((wavemat[1,i]),ymin,(wavemat[2,i]),ymax,border=scales:::alpha(plotcol,alfa),col=scales:::alpha(plotcol,alfa))
    }
  }
}