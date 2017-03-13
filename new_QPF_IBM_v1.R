
#******************************************************
#******************************************************
# The code is based on the paper:
# Paper information (submitted to the Annals of Applied Statistics)
# "A Spatio-Temporal Modeling Approach for Weather Radar Reflectivity Data and Its Applications"
# Author1: Xiao Liu (IBM Thomas J. Watson Research Center) 
# Author2: Vikneswaran Gopal (Department of Statistics and Applied Probability, National University of Singapore)
# Author3: Jayant Kalagnanam (IBM Thomas J. Watson Research Center) 
#******************************************************
#******************************************************

#******************************************************
#******************************************************
# Code information (available on Github: )
# Author1: Xiao Liu (IBM Thomas J. Watson Research Center) 
# Author2: Vikneswaran Gopal (Department of Statistics and Applied Probability, National University of Singapore)
#******************************************************
#******************************************************

#******************************************************
#******************************************************
# Contact information: Xiao Liu, liuxiao314923@gmail.com
#******************************************************
#******************************************************

#******************************************************
#******************************************************
# version 6
# Last modified: March 14 2017
#******************************************************
#******************************************************
rm(list = ls())    # clear workspace

#******************************************************
#******************************************************
# Control parameters
#******************************************************
#******************************************************
cluster <- FALSE  # set TRUE for parallel computing on cluster (using the R snow package)
                  # The code on GitHub only supports the computation on a single machine
if (!cluster){
  #wd <- "C:/Users/IBM_ADMIN/Desktop/IBM Research/Data and R Code/Weather Vik/liu_xiao/after Vik/QPF-I-IBM"
  wd <- "E:/old/IBM Back-up/IBM Research/Data and R Code/Weather Vik/liu_xiao/after Vik/QPF-I-IBM/QPF paper/code/GitHub"
}else{
  wd <- "/nfs/home/xiaoliu/Rcode/liu_xiao/QPF-I-IBM"
}
wd.data <- paste(wd,"/Monsoon2010",sep="") # specify the where the data are saved (not used for the version on GitHub)
wd.data = paste(wd,"/data",sep="") # the sample data is saved in the data folder
wd.figure = paste(wd,"/output",sep="") # the sample figures are saved in the data folder
Monsoon = "Monsoon2010"
setwd(wd)

#******************************************************
#******************************************************
# Packages needed
#******************************************************
#******************************************************
library(sp)
library(rgdal)
library(rgeos)
library(gstat)
library(RColorBrewer)
library(spacetime)
library(selextR) # a special packages created by Vik Gopal, which needs to be installed from "selextR_0.1.tar.gz"
library(timeSeries)
library(fields)
library(MASS)
library(snow)
library(mvtnorm)
library(matrixcalc)

#******************************************************
#******************************************************
# Load all sub-routines
#******************************************************
#******************************************************
source("SubRoutine.R")
#source("TREC_function.R")
source("TREC_function_P1.R")
source("TREC_function_P2.R")
source("new_QPF_IBM_v1_Kernel.R")

#******************************************************
#******************************************************
# Create processes if parallel computing is needed
#******************************************************
#******************************************************
n.process <- 2
cl <- makeCluster(n.process)

#******************************************************
#******************************************************
# Plotting info.
#******************************************************
#******************************************************
# set coordinates projection system: SVY21 is used for Singapore
SVY21 <- '+proj=tmerc +lat_0=1.366666666666667 +lon_0=103.8333333333333 +k=1 +x_0=28001.642 +y_0=38744.572 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
frequency <- 5 * 60 # Radar scan frequency, in sec
PLOT <- TRUE # if plots are generated during the code execution

singList <- list("sp.polygons", sgBd)   #
myList <- list("sp.polygons", myBd)
inList <- list("sp.polygons", inBd)
sp.layout.list <- list(singList, myList, inList)

colPalette <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'YlOrRd'))(100), .85)
colPalette2 <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'Blues'))(100), .85)
#colPalette2 <- colorRampPalette(brewer.pal(n=9, 'Greens'))(100)
colPalette1 <- rev(rainbow(32, start=5/6, end=3/6, alpha=0.8))
at.points <- seq(from=-35, to=75, length=33)

######################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################

#******************************************************
#******************************************************
# Read in the data
#******************************************************
#******************************************************
#setwd(wd.data)
filenames = list.files(path = wd.data)
n.storm = length(filenames) # number of storms in the folder
performance.list.A =performance.list.B= list()


for (storm.id in 1:n.storm){
  
  setwd(wd.data)
  if (!cluster){
    #scan <- readRDS("Monsoon2010/scan_list_23")  # load a storm event
    scan <- readRDS(filenames[storm.id])  # load a storm event
  }else{
    #scan <- readRDS("Monsoon2010/scan_list_23")
  }
  setwd(wd)
  
  if (FALSE){ # the following few lines are not executed. But the users may use these code to obtain more details of the radar scan and generate the plot
    class(scan)  # class: SelexScanList, a user-defined class by the package SelextR
    head(summary(scan),n=60)   # summary of radar scan data
    n.scan <- nrow(summary(scan))  # number of scans available for this storm event
    Time <- (summary(scan)[,1])[seq(1,n.scan,5)]  # times when radar scans are taken
    plot(scan[[81]], boundaries=TRUE)  # plot a scan
    str(scan[[1]], max=2)  
    class(scan[[1]]$data) # Format: "SpatialGridDataFrame"
  }
  
  
  #******************************************************
  #******************************************************
  # Process the timestamps
  #******************************************************
  #******************************************************
  # Select radar scans 5 mins apart; only applies to 1km height
  Time <- unique(summary(scan)[summary(scan)[,2]==1,1] )  # get the times when scans are taken; 1km above sea level only.
  # after loading the Time, need to check if they are 5min apart
  # drop the first and last
  Time <- strptime(Time,"%Y-%m-%d %H:%M:%S")
  Time <- as.POSIXct(Time, tz = "GMT") + 3600*8  # get the time in GMT+8, where Singapore is located
  mark <- array(0,dim=c(length(Time),1))
  for (i.time in 1:(length(Time)-1)){
    temp <- Time[i.time] + 60*5
    temp2 <- !is.na(match(temp,Time))
    if (temp2){
      mark[i.time] <- 1
    }
  }
  if (mark[(length(Time)-1)] == 1){
    mark[length(Time)] = 1
  }
  count <- array(0,dim=c(length(Time),1))
  for (i.time in 1: (length(Time))){
    sum <- 0
    check <- 1
    for (i.count in i.time:(length(Time))){
      if (mark[i.count]==0){
        check <- 0
      }
      sum <- sum + ( mark[i.count] * check) 
    }
    count[i.time] <- sum
  }
  
  
  Time.select <- Time[which(count == max(count)):( which(count == max(count))+max(count)-1)]
  scan.id.select <- which( (summary(scan)[,2]==1) & !is.na(match(summary(scan)[,1],as.character(Time.select-3600*8))) )
  n.scan <- length(scan.id.select)
  n.scan.end <- n.scan - 6
  
  # only run the code when the following conditions are satisfied
  C1 = (max(count) > 10)   
  C2 = ((n.scan-1) > 18)   # the number of available radar scan is at least 18 for a storm
  
  
  # the length of scan.select should be long enough; 6 scans are used for modeling construction, and 6 are used for prediction. 
  if (C1&C2){
    
    scan.select <- scan.id.select[2:n.scan.end ] # the code will run on selected scans. 
    
    # after obtaining the scan.select, specifiy two parts:
    # start and end of scan for model construction (7 scans are needed)
    #scan.select.id.start = 37 # in the numerical example, the value is set to 37
    #scan.select.id.end = scan.select.id.start+6 # in the numerical example, the value is set to 37+6
    scan.select.id.start = round(length(scan.select)/2)-6 # in the numerical example, the value is set to 37
    scan.select.id.end = scan.select.id.start+6 # in the numerical example, the value is set to 37+6
    
    #******************************************************
    #******************************************************
    # Basic Grid Settings
    #******************************************************
    #******************************************************
    grd <- SpatialPoints(coordinates(scan[[1]]$data),proj4string=CRS(SVY21))
    grd.p <- SpatialPixels(grd,proj4string=CRS(SVY21))
    n.grd <- length(grd)
    
    # number of TREC centers, i.e., number of boxes
    n.TREC.ctr <- 93*93 
    
    # TREC array center
    TREC.ctr.index <- array(0/0, dim=c(93,93) )
    for (i in 1:93){
      start <- (i-1)*5+9
      TREC.ctr.index[i,] <- seq( (480*start + 10), (480*(start+1)-10), 5)
    }
    grd.TREC <- grd[TREC.ctr.index]
    
    if (PLOT){ # visualize the box, if needed
      #pdf(file = "../figures/PixelArray.pdf",   width=4, height=4)
      setwd(wd.figure); pdf(file = "PixelArray.pdf",   width=4, height=4)
      par(mar=c(4,4,1,2))
      par(las=TRUE)
      plot(sgBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]),col="azure2")
      plot(myBd,add=TRUE)
      plot(inBd,add=TRUE)
      plot(grd.TREC,add=TRUE,col="azure3",lwd=0.1)
      dev.off()
      setwd(wd)
    }
    
    # create an area of interest (the area which is close to Singapore)
    grd.ctr <- SpatialPoints(cbind(28000,38000),proj4string=CRS(SVY21))
    # plot(sgBd,axes=TRUE);plot(grd.ctr,add=TRUE,col='red')
    dist <- ( coordinates(grd.TREC)[,1]-coordinates(grd.ctr)[1,1] )^2 + ( coordinates(grd.TREC)[,2]-coordinates(grd.ctr)[1,2] )^2
    select.grd <- (dist < 50000^2)  # select grids which are 50km radius from the specified center
    grd.SG <- grd.TREC[select.grd] # select grids
    # plot(sgBd,axes=TRUE,xlim=c(-20000,80000),ylim=c(-10000,90000));plot(grd.SG,add=TRUE);plot(grd.SG,axes=TRUE)
    # length(grd.SG)
    
    # station (information of weather monitoring stations; this part of information is not needed)
    if (FALSE){
      stn.data <- read.csv("latLongStationInfo.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
      stn.data <- stn.data[,c(1:4)]
      stationLat <- char2dms(stn.data$lat, chd="d", chm="m", chs="s")
      stationLong <- char2dms(stn.data$long, chd="d", chm="m", chs="s")
      stationLocation.sp <- SpatialPoints(cbind(as.numeric(stationLong), 
                                                as.numeric(stationLat)), CRS('+proj=longlat +datum=WGS84'))
      stationLocationXY <- spTransform(stationLocation.sp,CRS(SVY21))
      stn.data$X <- coordinates(stationLocationXY)[,1]
      stn.data$Y <- coordinates(stationLocationXY)[,2]
    }
    
    
    # get all grd id associated with each box, and save it in a column
    # each box contain 19*19 pixels
    box.id <- array(0/0, dim=c(19*19,93*93))
    for (i in 1:n.TREC.ctr){
      for (j in 1:19){
        box.id[  c(((j-1)*19+1):(j*19)),i] <- c(TREC.ctr.index[i]-4329 +(j-1)*480):(TREC.ctr.index[i]-4329 +(j-1)*480 + 18)
      }
    } 
    if (PLOT){
      
      #pdf(file = "../figures/PixelArray2.pdf",   width=12, height=6)
      setwd(wd.figure); pdf(file = "PixelArray2.pdf",   width=12, height=6)
      par(mfrow=c(1,2))
      par(mar=c(4,4,1.5,2))
      par(las=TRUE)
      plot(sgBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]),col="azure2")
      plot(myBd,add=TRUE)
      plot(inBd,add=TRUE)
      plot(grd.TREC,add=TRUE,col="black",lwd=0.1)
      title("(a)")
      
      x.min = min(  c(coordinates(grd[box.id[,2]])[,1],coordinates(grd[box.id[,1]])[,1] ) )
      x.max = max(  c(coordinates(grd[box.id[,2]])[,1],coordinates(grd[box.id[,1]])[,1] ) )
      y.min = min(  c(coordinates(grd[box.id[,2]])[,2],coordinates(grd[box.id[,1]])[,2] ) )
      y.max = max(  c(coordinates(grd[box.id[,2]])[,2],coordinates(grd[box.id[,1]])[,2] ) )+2000
      plot(coordinates(grd[box.id[,1]]),ylim=c(y.min,y.max),col="black",axes=TRUE,xlab="",ylab="",pch=3)
      points(coordinates(grd[box.id[,2]]),pch=0,col="blue")
      points(coordinates(grd.TREC[1]),col="red",pch=20,cex=2)
      points(coordinates(grd.TREC[2]),col="red",pch=20,cex=2)
      legend("topright",legend=c("pixels from the first pixel array",
                                 "pixels from the second pixel array",
                                 "centers of the two pixel arrays"),
             pch=c(3,0,20),col=c("black","blue","red"),bty="n")
      title("(b)")
      
      dev.off()
      setwd(wd)
    }
    
    
    # get all grd id associated with each box, and save it in a column (non-overlapping)
    # each box has 5by5 pixels
    box.id.nv <- array(0/0, dim=c(5*5,93*93))  # nv: non-overlapping
    for (i in 1:n.TREC.ctr){
      for (j in 1:5){
        box.id.nv[  c(((j-1)*5+1):(j*5)),i] <- (TREC.ctr.index[i]+480*(j-3)-2):(TREC.ctr.index[i]+480*(j-3)+2)
      }
    } 
    
    #plot(grd[box.id.nv[,1]]);plot(grd.TREC[1],add=TRUE,col="red",pch=20)
    #plot(grd[box.id.nv[,2]],col="red");plot(grd.TREC[2],add=TRUE,col="red",pch=20)
    
    
    #******************************************************
    #******************************************************
    # The Main Code:
    #******************************************************
    #******************************************************
    g.value <- array(0/0,dim=c(n.TREC.ctr,6))  # g.value: the growth and decay, to be computed
    box.track.list=mean.track.list=list() 
    # box.track.list: a list that tracks the TREC boxes
    # mean.track.list: track the mean dBZ on each box
    
    G.procedure <- FALSE  
    # if the growth/decay is computed; initial value should be FALSE; the value will be set to TRUE after the 3rd iterations
    box.select.track <- array(0/0,dim=c(n.TREC.ctr, 7))
    box.mean.track <- array(0/0,dim=c(n.TREC.ctr, 7))
    
    # initial
    scan.count <- which(scan.select == scan.select[scan.select.id.start]) # specify the starting scan
    data0 <- (scan[[scan.id.select[scan.count]]]$data)@data # initial data
    data0.TREC <- array(0/0,dim=c(n.TREC.ctr,1))
    for (i in 1:n.TREC.ctr){
      temp.data <- data0[box.id[,i],1]
      weight <- exp( (temp.data-75)/30 )
      #weight <- rep(1,length(temp.data))
      na.rm <- !is.na(temp.data)
      temp.data <- temp.data[na.rm]
      weight <- weight[na.rm]
      data0.TREC[i,1] <- sum(   temp.data*weight, na.rm=TRUE) / sum(weight)
    }  # aggregate data for each TREC box
    
    box.mean.track.initial  <- data0.TREC # initial dBZ on TREC boxes 
    box.select.initial <- c(1:n.TREC.ctr) # initial TREC box locations (indexed)
    
    mean.track <- matrix(box.mean.track.initial,ncol=1)
    box.track <- matrix(box.select.initial,ncol=1)
    
    scan.i <- 0
    output2.list <- list()
    for (i.scan in scan.select[scan.select.id.start:scan.select.id.end]){ 
      scan.i <- scan.i + 1
      
      print(paste(Monsoon,"_",filenames[storm.id],"_",i.scan,sep=""))
      scan.count <- which(scan.select == i.scan)
      Time.target <- Time.select[scan.count + 7]
      data <- (scan[[scan.id.select[scan.count]]]$data)@data
      data2 <- (scan[[scan.id.select[scan.count+1]]]$data)@data
      
      if (PLOT) {
        setwd(wd.figure);  pdf(file = paste(filenames[storm.id],"_observed_",scan.i,".pdf",sep=""), width=6, height=6)
        data.spdf <- SpatialPixelsDataFrame(grd.p,data.frame(data),proj4string=CRS(SVY21))
        plot(sgBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]),col="blue")
        #plot(grd,add=TRUE)
        image(data.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[11],at.points[length(at.points)]))
        title(Time.select[scan.count])
        dev.off()
        setwd(wd)
      }
      if (PLOT) {
        setwd(wd.figure);  pdf(file = paste(filenames[storm.id], "_observed_",scan.i+1,".pdf",sep=""), width=6, height=6)
        data2.spdf <- SpatialPixelsDataFrame(grd.p,data.frame(data2),proj4string=CRS(SVY21))
        plot(sgBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]),col="blue")
        #plot(grd,add=TRUE)
        image(data2.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[11],at.points[length(at.points)]))
        title(Time.select[scan.count+1])
        dev.off()
        setwd(wd)
      }
      
      grd.TREC.coords <- coordinates(grd.TREC)
      
      # compute the velocity vectors
      timing1 <- proc.time()
      output1 <- parSapply(cl, c(1:8649),TREC.parallel1,
                           Data=data, Data2=data2, N.TREC.ctr=n.TREC.ctr, Box.id=box.id, Grd.TREC=grd.TREC.coords)
      timing2 <- proc.time()-timing1
      
      box_select <- output1[2,]
      Cross_R <- output1[1,]
      TREC.u <- coordinates(grd.TREC)[box_select,1] - coordinates(grd.TREC)[,1]
      TREC.v <- coordinates(grd.TREC)[box_select,2] - coordinates(grd.TREC)[,2]
      
      NARM <- FALSE  # this code is used to get rid of NA values, by setting all NA to zero
      if (NARM){
        TREC.u[is.na(TREC.u)] <- 0
        TREC.v[is.na(TREC.v)] <- 0	  
      }
      
      # compute the smoothed velocity field
      output2 <- parSapply(cl, c(1:8649),TREC.parallel2,
                           Grd.TREC=grd.TREC.coords, CRoss_R=Cross_R, TREC.U=TREC.u, TREC.V=TREC.v)
      timing2 <- proc.time()-timing1
      #print(timing2)
      
      output2.list[[scan.i]] <- output2
      TREC.u.s <-  output2[2,]
      TREC.v.s <-  output2[3,]
      box_select_revised <- output2[1,]
      
      # need to modify box_select_revised
      #box_unselected <- which( is.na(match(c(1:8649),box_select_revised)) )
      #Dist <- array()
      #j <- 0
      #for (i in 2:length(box_select_revised)){
      #  if ( sum(box_select_revised[i]==box_select_revised[1:(i-1)])>0 ){
      #    j <- j+1
      #    dist <- (coordinates(grd.TREC[box_select_revised[i]])[,1]-coordinates(grd.TREC[box_unselected])[,1])^2+(coordinates(grd.TREC[box_select_revised[i]])[,2]-coordinates(grd.TREC[box_unselected])[,2])^2
      #   case <- which(dist==min(dist,na.rm=TRUE))[1]
      #   Dist[j] <- sqrt(min(dist,na.rm=TRUE))
      #   box_select_revised[i] <- box_unselected[case]
      #    box_unselected <- box_unselected[-case]
      #  }
      #}
      
      # Form a sequence of tracked features 
      # Get the data on TREC box
      data2.TREC <- array(0/0,dim=c(n.TREC.ctr,1))
      for (i in 1:n.TREC.ctr){
        temp.data <- data2[box.id[,i],1]
        weight <- exp( (temp.data-75)/30 )
        #weight <- rep(1, length(temp.data))
        na.rm <- !is.na(temp.data)
        temp.data <- temp.data[na.rm]
        weight <- weight[na.rm]
        data2.TREC[i,1] <- sum(   temp.data*weight, na.rm=TRUE) / sum(weight)
      }
      
      #track.revised <- output2[1,][box.track[,scan.i]] # the issue here is that some vectors point to the same location
      #for (i in 2:length(track.revised)){
      #  if ( sum(track.revised[i]==track.revised[1:(i-1)])>0 ){track.revised[i]<-0/0}
      #}
      
      box.select.track[,scan.i] <-  box_select_revised[ box.track[,scan.i]  ]
      box.mean.track[,scan.i]  <- data2.TREC[ box.select.track[,scan.i] ]
      
      mean.track <- cbind(box.mean.track.initial, box.mean.track)
      box.track <- cbind(box.select.initial, box.select.track)
      
      box.track.list[[scan.i]] <- box.track
      mean.track.list[[scan.i]] <- mean.track
      
      
      if (FALSE){
        par(mfrow=c(1,2))
        
        data2.spdf <- SpatialPixelsDataFrame(grd.p,data.frame(data),proj4string=CRS(SVY21))
        plot(sgBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]),col="blue")
        #plot(grd,add=TRUE)
        image(data2.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[11],at.points[length(at.points)]))
        arrows(coordinates(grd.TREC)[,1],coordinates(grd.TREC)[,2],
               coordinates(grd.TREC)[,1]+TREC.u.s,
               coordinates(grd.TREC)[,2]+TREC.v.s,length=0.05,col="blue",lwd=1)
        
        data2.spdf <- SpatialPixelsDataFrame(grd.p,data.frame(data2),proj4string=CRS(SVY21))
        plot(sgBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]),col="blue")
        #plot(grd,add=TRUE)
        image(data2.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[1],at.points[length(at.points)]))
        arrows(coordinates(grd.TREC[box.track[,scan.i]])[,1],coordinates(grd.TREC[box.track[,scan.i]])[,2],
               coordinates(grd.TREC[box.track[,scan.i+1]])[,1],coordinates(grd.TREC[box.track[,scan.i+1]])[,2],
               length=0.05,col="blue",lwd=1)
      }
      if (PLOT){  # Plot the wind field
        #mypath <- file.path( "../figures/RadarFieldWind3.png" )
        #png(mypath,width = 10, height = 6, units = "in", res=300)	
        setwd(wd.figure); png(file = paste(filenames[storm.id], "_windfield_",scan.i,".png",sep=""), width=10, height=6, units = "in", res=300)
      
        par(mfrow=c(1,2))
        par(mar=c(2,2,4,1))
        
        data2.spdf <- SpatialPixelsDataFrame(grd.p,data.frame(data),proj4string=CRS(SVY21))
        plot(sgBd,axes=TRUE,xlim=c(-10000,60000),col="blue")
        #plot(grd,add=TRUE)
        image(data2.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[11],at.points[length(at.points)]))
        arrows(coordinates(grd.TREC)[,1],coordinates(grd.TREC)[,2],
               coordinates(grd.TREC)[,1]+TREC.u.s,
               coordinates(grd.TREC)[,2]+TREC.v.s,length=0.04,col="blue",lwd=1)
        title(paste(Time.select[scan.count],"(zoomed-in)",sep=" "))
        
        data2.spdf <- SpatialPixelsDataFrame(grd.p,data.frame(data),proj4string=CRS(SVY21))
        plot(sgBd,axes=TRUE,xlim=c(-10000,60000),col="blue")
        #plot(grd,add=TRUE)
        image(data2.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[11],at.points[length(at.points)]))
        arrows(coordinates(grd.TREC)[,1],coordinates(grd.TREC)[,2],
               coordinates(grd.TREC)[,1]+TREC.u.s,
               coordinates(grd.TREC)[,2]+TREC.v.s,length=0.04,col="blue",lwd=1)
        title(paste(Time.select[scan.count],"(zoomed-in)",sep=" "))
        dev.off()
        setwd(wd)
      }
      
      
      if (G.procedure){ # when scan.count is >= 2, we can start create the g.value matrix for growth and decay	
        i <- scan.i
        case <- !is.na(box.track[,i])
        grd.temp <- grd.TREC[box.track[case,i]]
        #temp1 <- box.track[,i-1]
        #temp2 <- box.track[,i+1] 
        g.value[,scan.i-1] <- ( mean.track[,i+1]-mean.track[,i-1] ) / 2
      }
      
      G.procedure <- TRUE
      
    } #for (i.scan in scan.select[scan.select.id.start:scan.select.id.end])
    
    
    
    print(paste(Monsoon,"_",filenames[storm.id], sep=""))
    print("proceed to g")
    
    # *********************************************
    # *********************************************
    # The following code is for the construction of the STCAR model for the growth and decay process, G
    # *********************************************
    # *********************************************
    
    i = 6
    case1 <- which( !is.na(box.track[,i+1]) )
    case2 <- which( (g.value[,i]!=0) )
    case <- intersect(case1, case2)
    grd.temp <- grd.TREC[box.track[,i+1][case]]
    value.spdf <- SpatialPixelsDataFrame(grd.temp,data.frame(g.value[case,i]),proj4string=CRS(SVY21))
    if (FALSE){
      plot(sgBd,axes=TRUE,xlim=c(bbox(sgBd)[1,1],bbox(sgBd)[1,2]),ylim=c(bbox(sgBd)[2,1],bbox(sgBd)[2,2]),col="blue")
      image(value.spdf,add=TRUE,col=colPalette1,zlim=c(-32,30))
      hist(g.value[ case2,i],30, col="grey")
      qqnorm(g.value[case2,i]);qqline(g.value[case2,i])
    }
    
    
    # create response, y, in the regression model
    n.seq <- 6
    y = mu = list()
    for (i in 1:n.seq){
      case1 <- which( !is.na(box.track[,i+1]) )
      case2 <- which( (g.value[,i]!=0) )
      case <- intersect(case1, case2)
      
      Grd.temp <- grd.TREC[box.track[,i+1][case]]
      Value <- g.value[case,i]
      
      tmp <- k.space.smooth1(grd=Grd.temp, value=Value, k.var=100, J=30 )  
      mu[[i]] <- tmp[[1]][,1]
      y[[i]] <- tmp[[1]][,2]
      #print(i)
    }
    
    # statistical inference for unknown parameters using the method described in Section 3.3 
    size = 1000
    Rou.hat = Sigma.hat = array()
    B.hat.list = Wd.list = Wd.inv.list = list()
    var.est.list = gt.list = yt.list = list()
    Yt = Mut = Select.grd =  array(0/0, dim=c(size,6))
    
    grd.ctr <-  SpatialPoints(cbind(28000,38000),proj4string=CRS(SVY21))
    select <- array()
    for (i in 1:8649){
      tmp <- (!is.na(g.value[i,1])) & (!is.na(box.track.list[[1+1]][i,2]))
      if (tmp){
        grd.temp <- grd.TREC[box.track.list[[1+1]][i,2]]
        dist[i] <- ( coordinates(grd.temp)[,1]-coordinates(grd.ctr)[1,1] )^2 + ( coordinates(grd.temp)[,2]-coordinates(grd.ctr)[1,2] )^2
        select[i] <- i
      }else{
        dist[i] <- Inf
      }
    }
    select.grd <- which( rank(dist,ties.method = c("random"))<=size)
    grd.temp <- grd.TREC[box.track.list[[1+1]][select.grd ,2]]
    
    grd.temp2 <- grd.TREC
    dist <- ( coordinates(grd.temp2)[,1]-coordinates(grd.ctr)[1,1] )^2 + ( coordinates(grd.temp2)[,2]-coordinates(grd.ctr)[1,2] )^2
    select.grd2 <- which( rank(dist,ties.method = c("random"))<=size)
    grd.temp2 <- grd.TREC[select.grd2]
    
    for (g.col in 1:ncol(g.value)){
      
      #grd.ctr <-  SpatialPoints(cbind(28000,38000),proj4string=CRS(SVY21))
      #case <- (!is.na(g.value[,g.col])) & (!is.na(box.track.list[[g.col+1]][,2]))
      mean.growth.temp <- g.value[select.grd,g.col]
      #grd.temp <- grd.TREC[box.track.list[[g.col+1]][case,2]]
      #plot(grd.temp,axes=TRUE);plot(grd.SG,add=TRUE,col="red")	
      #plot(sgBd,axes=TRUE);plot(grd.ctr,add=TRUE,col='red')
      #dist <- ( coordinates(grd.temp)[,1]-coordinates(grd.ctr)[1,1] )^2 + ( coordinates(grd.temp)[,2]-coordinates(grd.ctr)[1,2] )^2
      #select.grd <- (dist < 50000^2) 
      #select.grd <- which( rank(dist,ties.method = c("random"))<=size) 
      #grd.temp <- grd.temp[select.grd]
      # plot(sgBd,axes=TRUE,xlim=c(-20000,80000),ylim=c(-10000,90000));plot(grd.temp,add=TRUE)
      gt <- matrix( mean.growth.temp,ncol=1) 
      
      # impute the NA values
      gt.na <- which(gt==0)
      for (i in gt.na){
        gt[i] <- mean( gt[ which( abs(select.grd-select.grd[i])<=100 ) ] , na.rm=TRUE )
      }
      
      W <- array(0,dim=c(length(grd.temp),length(grd.temp)))
      for (row in 1:length(grd.temp)){
        grd.select <- grd.temp[row]
        temp <- (coordinates(grd.select)[,1] - coordinates(grd.temp)[,1])^2 + (coordinates(grd.select)[,2] - coordinates(grd.temp)[,2])^2 
        dist <- sqrt( temp )/1000 # KM
        #case <- (dist < 50)
        #W[row,case] <- 1
        #W[row,case] <- exp(-dist[case]/10)
        W[row,] <- exp(-dist/10)
        case <- (dist == 0)
        W[row,case] <- 0
      }
      Wd.inv <- diag( 1/rowSums(W), length(grd.temp), length(grd.temp))
      Wd <- diag( rowSums(W), length(grd.temp), length(grd.temp))
      Wd.list[[g.col]] <- Wd
      Wd.inv.list[[g.col]] <- Wd.inv
      
      beta <- Wd.inv %*% W
      #beta <- W
      #for (j in 1:nrow(W)){
      #  beta[j,] <- W[j,] * Wd.inv[j,j]
      #}
      eigen.beta <- eigen(beta, symmetric=TRUE, only.values = TRUE, EISPACK = FALSE)$values
      #is.positive.definite(beta)  
      
      # compute beta using LSE; in fact, beta should be estimated using MLE; here a simplied method is adopted
      tmp <- k.space.smooth1(grd=grd.temp, value=gt, k.var=100, J=30)  
      yt <- matrix( tmp[[1]][,2], ncol=1 )
      gamma <- tmp[[2]]
      
      judge <- 1
      jj <- 0 
      rou.hat = sigma.hat = delta.change = array()
      while (judge == 1){
        #print(jj)
        jj <- jj + 1
        likelihood <- function(rou){
          B <- diag(1,nrow(beta)) - rou * beta
          detB <- det(B)
          n <- nrow(beta)
          #sigma.square <- t(yt) %*% t(B) %*% B %*% yt/n
          sigma.square <- t(yt) %*% Wd %*% B %*% yt/n
          #output <- log(sigma.square) - 2/n * log(detB)
          #output <- log(sigma.square) - 2/n * determinant(Wd %*% B, logarithm = TRUE)$modulus
          output <- log(sigma.square) - 1/n * determinant(Wd %*% B, logarithm = TRUE)$modulus
          return(output)
        }
        reltol <- 1e-6
        trace <- 46
        #print("mle started")
        
        temp <- optim(0, likelihood, control=list(reltol,trace), method="L-BFGS-B", 
                      lower=1/min(eigen.beta), upper=1/max(eigen.beta))
        #temp <- optim(0, likelihood, method="L-BFGS-B", lower=-1, upper=1)
        rou.hat[jj] <- temp$par
        B.hat <- diag(1,nrow(beta)) - rou.hat[jj] * beta
        detB.hat <- det(B.hat)
        n <- nrow(beta)
        #sigma.hat[jj] <- sqrt(  t(yt) %*% t(B.hat) %*% B.hat %*% yt/n  )   
        sigma.hat[jj] <- sqrt(  t(yt) %*% Wd %*% B.hat %*% yt/n  )   
        
        SIGMA <-   solve(B.hat) %*% (sigma.hat[jj]^2 * Wd.inv)
        #is.positive.definite(round(SIGMA,2));
        tmp <- k.space.smooth2(grd=grd.temp, grd2=grd.temp2, value=gt, k.var=100, J=30, var_y=round(SIGMA,2))  
        
        yt <- matrix( tmp[[1]][,2], ncol=1 )
        mut <- matrix( tmp[[1]][,3], ncol=1 )
        gamma <- cbind(gamma, tmp[[2]])
        
        delta.change[jj] <- mean( (gamma[,jj+1]-gamma[,jj])^2 )
        if (jj == 1){
          judge <- 0	
        }
      }
      
      Rou.hat[g.col] <- rou.hat[length(rou.hat)]
      Sigma.hat[g.col] <- sigma.hat[length(sigma.hat)]
      B.hat.list[[g.col]] <- diag(1,nrow(beta)) - Rou.hat[g.col] * beta
      Yt[,g.col] <- yt
      Mut[,g.col] <- mut
    }
    
    
    # estimate r
    order.p = 4
    y.ar <- array(0,dim=c(2000,1))
    x.ar <- array(0/0,dim=c(2000,4))
    Cov.mat <- array(0, dim=c(2000,2000))
    ii <- 0
    for (i in 5:6){
      ii <- ii + 1
      
      Gt1 <- g.value[select.grd,i]
      Gt2 <- g.value[select.grd,i-1]
      Gt3 <- g.value[select.grd,i-2]
      Gt4 <- g.value[select.grd,i-3]
      Gt5 <- g.value[select.grd,i-4]
      
      cov.mat <-  Sigma.hat[i]^2 * Wd.inv.list[[i]]
      Cov.mat[((ii-1)*1000+1):(ii*1000),  ((ii-1)*1000+1):(ii*1000)] <- cov.mat
      y.ar[((ii-1)*1000+1):(ii*1000)] <- matrix( Gt1,ncol=1)
      x.ar[((ii-1)*1000+1):(ii*1000),] <- cbind(Gt2,Gt3,Gt4,Gt5)
      
    }
    
    lm.data <- data.frame( cbind(y.ar, x.ar) )
    colnames(lm.data) <- c("y",paste("X",c(1:ncol(x.ar)),sep=""))
    xnam <- paste("X", c(1:ncol(x.ar)),sep="")
    lm.gls.fit <- lm.gls(as.formula(paste("y~", paste(xnam, collapse= "+"),"-1",sep="")),
                         data=lm.data,W=Cov.mat,na.action=na.exclude,inverse=TRUE,method = "qr")
    r <- lm.gls.fit$coefficients
    #qqnorm(lm.gls.fit$residuals)
    
    #save.image(paste(Monsoon,"_",filenames[storm.id],"_","gvalue.RData",sep=""))
    setwd(wd.figure)
    save.image(paste(filenames[storm.id],"_","gvalue.RData",sep="")) # save the entire data set as .RData
    setwd(wd)
    #load(paste(Monsoon,"_",filenames[storm.id],"_","gvalue.RData",sep=""))
    

    # *********************************************
    # *********************************************
    # The following code is for the construction of the STCAR model for the growth and decay process, G
    # *********************************************
    # *********************************************
    print("proceed to prediction")
    # source the prediction subroutine
    source("new_QPF_IBM_v1_pred.R")
    
    # compute the performance in terms of MSE
    performance_A = data.frame( cbind(MSE1,MSE2,MSE1a,MSE2a) )  # raw MSE
    performance_B = data.frame( cbind(MSE1b,MSE2b,MSE1c,MSE2c) )  # smoothed MSE
    colnames(performance_A) = c("MSE1", "MSE2", "MSE1a", "MSE2a")
    colnames(performance_B) = c("MSE1", "MSE2", "MSE1a", "MSE2a")
    performance.list.A[[storm.id]] = performance_A
    performance.list.B[[storm.id]] = performance_B
    
    
  } # if (length(scan.select)>=18)
  
  setwd(wd.figure)
  save(performance.list.A, file=paste(Monsoon,"_test_A",".RData",sep="") )
  save(performance.list.B, file=paste(Monsoon,"_test_B",".RData",sep="") )
  setwd(wd)
  print(paste(Monsoon,"_",filenames[storm.id], sep=""))
  print("completed")
  
}






























