# Script for processing X,Y coordinates from Poplar tension wood bending to extract curvature

bezierCurve <- function(x, y, n=10)
{
  outx <- NULL
  outy <- NULL
  
  i <- 1
  for (t in seq(0, 1, length.out=n))
  {
    b <- bez(x, y, t)
    outx[i] <- b$x
    outy[i] <- b$y
    
    i <- i+1
  }
  
  #return (list(x=outx, y=outy))
  return(cbind(outx,outy))
}

bez <- function(x, y, t)
{
  outx <- 0
  outy <- 0
  n <- length(x)-1
  for (i in 0:n)
  {
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
  }
  
  return (list(x=outx, y=outy))
}


# grep_key=".txt"
# hook=1800
# j=1
# output_name="out"


StemCurve<-function(directory,grep_key=".txt",hook=1800,len=NA)
{
  require(sp)
  
  all_files<-list.files(directory,full.names=TRUE)
  subset<-grep(grep_key,all_files,perl=TRUE)
  
  files<-all_files[subset]
  n<-length(files)
  
  #check if all files are .txt
  check1<-length(grep(grep_key, files))
  
  if(n!=check1) {print("Not all files are .txt") } else
    
    output <- data.frame(matrix(vector(), n, 2, dimnames=list(c(), c("curvature", "height"))), stringsAsFactors=F)
  out<-data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c("file","x","y"))), stringsAsFactors=F)
  vec<-data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("file","xis","yis","angle"))), stringsAsFactors=F)
  
  for (j in 1:n)
  {
    wt_03.13<-NULL
    outi<-NULL
    veci<-NULL
    data <- read.table(files[j], sep='\t')
    data[,2]<-2500 - data[,2]
    
    #crop line to distance to hook
    #returns the segment lengths
    S <- SpatialPoints(data)
    seg<-LineLength(Line(S),sum=F)
    sum<-0
    for(k in 1:length(seg))
    {
      sum<-sum+seg[k]
      if(sum>=hook) break;
    }
    
    wt_03.13<-data[1:k,]
    outi<-cbind(files[j],data[1:k,])
    
    
    names(wt_03.13) <- c('x','y') # relabel columns for clarity
    names(outi)<-c("file","x","y")
    
    wt_03.13$x <- wt_03.13$x - wt_03.13$x[1]
    wt_03.13$y <- wt_03.13$y - wt_03.13$y[1] # y-axis is inverted in ImageJ by default
    head(wt_03.13,30)
    
    # calculate vectors
    xis=NULL
    for(i in 2:length(wt_03.13$x)) xis[i-1] <- wt_03.13$x[i] - wt_03.13$x[i-1]
    yis=NULL
    for(i in 2:length(wt_03.13$y)) yis[i-1] <- wt_03.13$y[i] - wt_03.13$y[i-1]
    vectors <- data.frame(cbind(xis, yis))
    #head(vectors,20)
    
    #http://www.euclideanspace.com/maths/algebra/vectors/angleBetween/index.htm
    # calculate angles
    vectors$angle=NULL
    for (i in 1:length(vectors$xis)-1) vectors$angle[i] <- (atan2(vectors$yis[i+1], vectors$xis[i+1]) - atan2(vectors$yis[i], vectors$xis[i]))*180/pi
    #head(vectors,20)
    #vectors
    veci<-cbind(files[j],vectors)
    names(veci)<-c("file","xis","yis","angle")
    names(vec)<-c("file","xis","yis","angle")
    vec<-rbind(vec,veci)
    
    #total change in angle
    curv <- sum(vectors$angle)
    #curv # 1.296412 radians
    if(is.na(len)) 
    {
      output[j,1]<-curv
      output[j,2]<-wt_03.13$y[nrow(wt_03.13)]/sum(seg) 
    } else 
    {
      #output[j,1]<-curv
      #output[j,2]<-wt_03.13$y[nrow(wt_03.13)]/len
      #getfinal xy corrdinate
      output[j,1]<-data[k,1]
      output[j,2]<-data[k,2]
    }
    
    row.names(output)[j]<-files[j]
    names(out)<-c("file","x","y")
    out<-rbind(out,outi)
  }
  
  
  return(list(stat=output,data=out,vectors=vec))
  
  
}




BelzerCurve <- function(directory,grep_key=".txt",hook=1800,npoints=20)
{
  require(sp)
  #setwd(directory)
  
  all_files<-list.files(directory,full.names=T)
  subset<-grep(grep_key,all_files,perl=TRUE)
  
  files<-all_files[subset]
  n<-length(files)
  
  #check if all files are .txt
  check1<-length(grep(grep_key, files))
  
  if(n!=check1) {print("Not all files are .txt") } else
    
    output <- data.frame(matrix(vector(), n, 2, dimnames=list(c(), c("curvature", "height"))), stringsAsFactors=F)
  out<-data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c("file","x","y"))), stringsAsFactors=F)
  
  
  for (j in 1:n)
  {
    wt_03.13<-NULL
    outi<-NULL
    
    data <- read.table(files[j], sep='\t')
    data[,2]<-2500 - data[,2]
    
    #crop line to distance to hook
    #returns the segment lengths
    S <- SpatialPoints(data)
    seg<-LineLength(Line(S),sum=F)
    sum<-0
    for(k in 1:length(seg))
    {
      sum<-sum+seg[k]
      if(sum>=hook) break;
    }
    
    
    
    #fit with belzer spline curve
    
    wt_03.13<-data.frame(bezierCurve(data[1:k,1],data[1:k,2],npoints))
    outi<-cbind(files[j],wt_03.13)
    names(wt_03.13) <- c('x','y') # relabel columns for clarity
    names(outi)<-c("file","x","y")
    
    wt_03.13$x <- wt_03.13$x - wt_03.13$x[1]
    wt_03.13$y <- wt_03.13$y - wt_03.13$y[1] # y-axis is inverted in ImageJ by default
    head(wt_03.13,30)
    
    # calculate vectors
    xis=NULL
    for(i in 2:length(wt_03.13$x)) xis[i-1] <- wt_03.13$x[i] - wt_03.13$x[i-1]
    yis=NULL
    for(i in 2:length(wt_03.13$y)) yis[i-1] <- wt_03.13$y[i] - wt_03.13$y[i-1]
    vectors <- data.frame(cbind(xis, yis))
    #head(vectors,20)
    
    # calculate angles
    vectors$angle=NULL
    for (i in 1:length(vectors$xis)-1) vectors$angle[i] <- atan2(vectors$yis[i+1], vectors$xis[i+1]) - atan2(vectors$yis[i], vectors$xis[i])
    #head(vectors,20)
    
    curv <- sum(vectors$angle)
    #curv # 1.296412 radians
    output[j,1]<-curv#/sum(seg)
    output[j,2]<-wt_03.13$y[nrow(wt_03.13)]/sum(seg)
    row.names(output)[j]<-files[j]
    names(out)<-c("file","x","y")
    out<-rbind(out,outi)
  }
  
  
  return(list(stat=output,data=out))
  
  #write.table(output,file=paste(directory,"/",output_name,sep=""),sep="\t")
  
}


