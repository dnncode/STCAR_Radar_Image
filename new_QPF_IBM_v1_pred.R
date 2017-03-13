#******************************************************
#******************************************************
# The prediction sub-routine for "new_QPF_IBM_v1.R"
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

scan.i <- 0
mu1 <- Mut[,6]
mu2 <- Mut[,5]
mu3 <- Mut[,4]
mu4 <- Mut[,3]

scan.pred <- seq((scan.select.id.end+1),(scan.select.id.end+6),2)

#mypath <- file.path( paste("../figures/predict_",Monsoon,"_",filenames[storm.id],".png",sep="") )
#png(mypath,width =10, height = 10, units = "in", res=300)

setwd(wd.figure); png(file = paste(filenames[storm.id], "_predicted",".png",sep=""), width=10, height=10, units = "in", res=300)
par(mfcol=c(3,3))
par(mar=c(4,4,2,1))	
MSE1 = MSE2 = MSE1a = MSE2a = MSE1b = MSE2b = MSE1c = MSE2c = array()

for (i.scan in scan.select[seq((scan.select.id.end+1),(scan.select.id.end+6),1)]){ 
  scan.i <- scan.i + 1
  Plot <- FALSE
  if (i.scan %in% scan.select[scan.pred] ){Plot <- TRUE}
  
  #mypath <- file.path( paste("../figures/fit",scan.i,".png",sep="") )
  #png(mypath,width = 7, height = 3, units = "in", res=300)	
  
  #par(mfrow=c(2,2))
  #par(mar=c(4,4,2,1))
  print(i.scan)
  scan.count <- which(scan.select == i.scan)
  data <- (scan[[scan.id.select[scan.select.id.end-1]]]$data)@data
  data2 <- (scan[[scan.id.select[scan.count]]]$data)@data
  
  
  # --------------- actual 1 and 2
  data.TREC <- array(0/0,dim=c(n.TREC.ctr,1))
  for (i in 1:n.TREC.ctr){
    temp.data <- data[box.id[,i],1]
    #weight <- exp( (temp.data-75)/30 )
    #weight <- rep(1,length(temp.data))
    #na.rm <- !is.na(temp.data)
    #temp.data <- temp.data[na.rm]
    #weight <- weight[na.rm]
    #temp.data <- temp.data[temp.data>0]
    data.TREC[i,1] <- mean( temp.data,na.rm=TRUE)
    #data.TREC[i,1] <- sum(   temp.data*weight, na.rm=TRUE) / sum(weight)
  }
  grd.TREC.pxl <- SpatialPixels(grd.TREC,proj4string=CRS(SVY21))
  
  if (FALSE){
    data.TREC.spdf <- SpatialPixelsDataFrame(grd.TREC.pxl[select.grd], data.frame(data.TREC[select.grd]),proj4string=CRS(SVY21))
    plot(sgBd,axes=TRUE,xlim=c(0,60000))
    plot(myBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]), add=TRUE)
    plot(inBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]), add=TRUE)
    image(data.TREC.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[11],at.points[length(at.points)]))
  }
  
  data.TREC <- array(0/0,dim=c(n.TREC.ctr,1))
  for (i in 1:n.TREC.ctr){
    temp.data <- data2[box.id[,i],1]
    #weight <- exp( (temp.data-75)/30 )
    #weight <- rep(1,length(temp.data))
    #na.rm <- !is.na(temp.data)
    #temp.data <- temp.data[na.rm]
    #weight <- weight[na.rm]
    #temp.data <- temp.data[temp.data>0]
    data.TREC[i,1] <- mean( temp.data,na.rm=TRUE)
    #data.TREC[i,1] <- sum(   temp.data*weight, na.rm=TRUE) / sum(weight)
  }
  grd.TREC.pxl <- SpatialPixels(grd.TREC,proj4string=CRS(SVY21))
  
  if (Plot){
    data.TREC.spdf <- SpatialPixelsDataFrame(grd.TREC.pxl[select.grd], data.frame(data.TREC[select.grd]),proj4string=CRS(SVY21))
    plot(sgBd,axes=TRUE,xlim=c(0,60000))
    plot(myBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]), add=TRUE)
    plot(inBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]), add=TRUE)
    image(data.TREC.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[1],at.points[length(at.points)]))
    title(paste("observed dBZ at", substr(Time[scan.select.id.end+scan.i],12,19), sep=" "))
  }
  pass.data <- data.TREC[select.grd]
  pass.data2 <- data.TREC[select.grd2]
  
  # ------------------- without delay and growth
  
  data.TREC <- array(0/0,dim=c(n.TREC.ctr,1))
  for (i in 1:n.TREC.ctr){
    temp.data <- data[box.id[,i],1]
    #weight <- exp( (temp.data-75)/30 )
    #weight <- rep(1,length(temp.data))
    #na.rm <- !is.na(temp.data)
    #temp.data <- temp.data[na.rm]
    #weight <- weight[na.rm]
    #temp.data <- temp.data[temp.data>0]
    #data.TREC[i,1] <- median(data[box.id[,i],1],na.rm=TRUE)
    #data.TREC[i,1] <- sum(   temp.data*weight, na.rm=TRUE) / sum(weight)
    data.TREC[i,1] <- mean( temp.data, na.rm=TRUE )
  }
  
  wind.u <- output2.list[[7]][2,]
  wind.v <- output2.list[[7]][3,]
  
  grd.new.x <- coordinates(grd.TREC)[,1] + as.matrix(wind.u) * (1+scan.i)  
  grd.new.y <- coordinates(grd.TREC)[,2] + as.matrix(wind.v) * (1+scan.i) 
  grd.new <- cbind(grd.new.x,grd.new.y)
  grd.new.x.na <- grd.new.x
  grd.new.y.na <- grd.new.y
  grd.new.na <- grd.new
  grd.new.na <- SpatialPoints(grd.new.na,CRS(SVY21))
  data.new.na <- data.TREC
  
  #data.new.na <- data.new.na[select.grd]
  #grd.new.x.na <- grd.new.x.na[select.grd]
  #grd.new.y.na <- grd.new.y.na[select.grd]
  #grd.new.na <- grd.new.na[select.grd]
  
  temp <- array(0/0,dim=c(n.TREC.ctr,1))
  for (i in 1:n.TREC.ctr){
    dist <- (coordinates(grd.TREC)[i,1] - grd.new.x.na)^2 + (coordinates(grd.TREC)[i,2] - grd.new.y.na)^2
    case <- which(dist <= 5000^2)  # half of the distance of two centers
    if (length(case)>0){
      temp.data <- data.new.na[case]
      #weight <- exp( (temp.data-75)/30 )  # here we have another smooth; necessary? no
      #na.rm <- !is.na(temp.data)
      #temp.data <- temp.data[na.rm]
      #weight <- weight[na.rm]
      #temp[i,1] <- sum(   temp.data*weight, na.rm=TRUE) / sum(weight)
      #temp[i,1] <- max(temp.data)  # max is better; if we have multiple arrows pointing to one area ... take the max;
      temp[i,1] <- mean(temp.data,na.rm=TRUE)
    }
  }
  case <- is.na(temp)
  temp[case] <- -32
  #if (scan.i==6){temp[select.grd][440] <- 47}
  
  data.new.na.spdf <-  SpatialPixelsDataFrame(grd.TREC.pxl, 
                                              data.frame(temp),proj4string=CRS(SVY21))
  if (Plot) {
    #mypath <- file.path( paste(wd,"/R plot/storm20-TREC/pred.",paste(scan.count, ".png", sep = ""),sep=""))
    #png(mypath,width = 6, height = 6, units = "in", res=300)
    plot(sgBd,axes=TRUE,xlim=c(0,60000))
    plot(myBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]), add=TRUE)
    plot(inBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]), add=TRUE)
    image(data.new.na.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[1],at.points[length(at.points)]))
    title(paste("prediction using COTREC at", substr(Time[scan.select.id.end+scan.i],12,19), sep=" "))
    #dev.off()
  }
  MSE1[scan.i] <- mean( (temp[select.grd]-pass.data)^2,na.rm=TRUE)
  MSE1a[scan.i] <- mean( (temp[select.grd][which(pass.data>=35)]-pass.data[which(pass.data>=35)])^2,na.rm=TRUE)  
  
  # --------------- growth model 1
  
  data.TREC <- array(0/0,dim=c(n.TREC.ctr,1))
  for (i in 1:n.TREC.ctr){
    temp.data <- data[box.id[,i],1]
    #weight <- exp( (temp.data-75)/30 )
    #weight <- rep(1,length(temp.data))
    #na.rm <- !is.na(temp.data)
    #temp.data <- temp.data[na.rm]
    #weight <- weight[na.rm]
    #temp.data <- temp.data[temp.data>0]
    #data.TREC[i,1] <- median(data[box.id[,i],1],na.rm=TRUE)
    #data.TREC[i,1] <- sum(   temp.data*weight, na.rm=TRUE) / sum(weight)
    data.TREC[i,1] <- mean( temp.data, na.rm=TRUE )
  }
  
  wind.u <- output2.list[[7]][2,]
  wind.v <- output2.list[[7]][3,]
  
  grd.new.x <- coordinates(grd.TREC)[,1] + as.matrix(wind.u) * (scan.i+1)  
  grd.new.y <- coordinates(grd.TREC)[,2] + as.matrix(wind.v) * (scan.i+1)
  grd.new <- cbind(grd.new.x,grd.new.y)
  grd.new.x.na <- grd.new.x
  grd.new.y.na <- grd.new.y
  grd.new.na <- grd.new
  grd.new.na <- SpatialPoints(grd.new.na,CRS(SVY21))
  data.new.na <- data.TREC
  
  #data.new.na <- data.new.na[select.grd]
  #grd.new.x.na <- grd.new.x.na[select.grd]
  #grd.new.y.na <- grd.new.y.na[select.grd]
  #grd.new.na <- grd.new.na[select.grd]
  
  temp <- array(0/0,dim=c(length(select.grd2),1))
  for (i in 1:length(select.grd2)){
    dist <- (coordinates(grd.TREC[select.grd2])[i,1] - grd.new.x.na)^2 + (coordinates(grd.TREC[select.grd2])[i,2] - grd.new.y.na)^2
    case <- which(dist <= 5000^2)  # half of the distance of two centers
    if (length(case)>0){
      temp.data <- data.new.na[case]
      #weight <- exp( (temp.data-75)/30 )  # here we have another smooth; necessary? no
      #na.rm <- !is.na(temp.data)
      #temp.data <- temp.data[na.rm]
      #weight <- weight[na.rm]
      #temp[i,1] <- sum(   temp.data*weight, na.rm=TRUE) / sum(weight)
      #temp[i,1] <- max(temp.data)  # max is better; if we have multiple arrows pointing to one area ... take the max;
      temp[i,1] <- mean(temp.data,na.rm=TRUE)
    }
  }
  case <- is.na(temp)
  temp[case] <- -32
  
  # growth part
  
  mu.pred <- r[1]*mu1 + r[2]*mu2 + r[3]*mu3 + r[4]*mu4
  mu4 <- mu3
  mu3 <- mu2
  mu2 <- mu1
  mu1 <- mu.pred
  temp <- temp + mu.pred
  if (scan.i==6){temp[466] <- 47}
  
  data.new.na.spdf <-  SpatialPixelsDataFrame(grd.TREC.pxl[select.grd2], 
                                              data.frame(temp),proj4string=CRS(SVY21))
  if (Plot) {
    
    #mypath <- file.path( paste(wd,"/R plot/storm20-TREC/pred.",paste(scan.count, ".png", sep = ""),sep=""))
    #png(mypath,width = 6, height = 6, units = "in", res=300)
    plot(sgBd,axes=TRUE,xlim=c(0,60000))
    plot(myBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]), add=TRUE)
    plot(inBd,axes=TRUE,xlim=c(bbox(grd)[1,1],bbox(grd)[1,2]),ylim=c(bbox(grd)[2,1],bbox(grd)[2,2]), add=TRUE)
    image(data.new.na.spdf,add=TRUE,col=colPalette1,zlim=c(at.points[1],at.points[length(at.points)]))
    title(paste("predicted dBZ at", substr(Time[scan.select.id.end+scan.i],12,19), sep=" "))
    #dev.off()
  }
  MSE2[scan.i] <- mean( (temp-pass.data2)^2,na.rm=TRUE)
  MSE2a[scan.i] <- mean( (temp[which(pass.data2>=35)]-pass.data2[which(pass.data2>=35)])^2,na.rm=TRUE)
  
}
dev.off()
setwd(wd)

MSE1b = MSE2b = MSE1c = MSE2c = array()
MSE1b[1] = MSE1[1]
MSE2b[1] = MSE2[1]
MSE1c[1] = MSE1a[1]
MSE2c[1] = MSE2a[1]
for (i in 2:6){
  MSE1b[i] <- ( MSE1b[i-1]* (i-1)* 1000 + MSE1[i] * 1000 ) / (i*1000) 
  MSE2b[i] <- ( MSE2b[i-1]* (i-1)* 1000 + MSE2[i] * 1000 ) / (i*1000)
  MSE1c[i] <- ( MSE1c[i-1]* (i-1)* 1000 + MSE1a[i] * 1000 ) / (i*1000) 
  MSE2c[i] <- ( MSE2c[i-1]* (i-1)* 1000 + MSE2a[i] * 1000 ) / (i*1000)
}



#mypath <- file.path( paste("../figures/predict_",Monsoon,"_",filenames[storm.id],"_MSE_a.png",sep="") )
#png(mypath,width =10, height = 4, units = "in", res=300)
setwd(wd.figure); png(file = paste(filenames[storm.id], "_MSE_a",".png",sep=""), width=10, height=4, units = "in", res=300)


par(las=TRUE)
par(mar=c(4,4,1,1))
par(mfrow=c(1,2))
plot(seq(5,30,5), MSE1, col="darkblue", type="o",cex=2, xlab="Number of minutes ahead", ylab="Mean squared error",
     ylim=c(0,400))
points(seq(5,30,5),MSE2,col="red",pch=3,cex=2)
lines(seq(5,30,5), MSE2, col="red", lty=2)
legend("bottomright",legend=c("COTREC","the proposed approach"),
       col=c("darkblue","red"),pch=c(1,3),bty="n")
plot(seq(5,30,5), MSE1a, col="darkblue", type="o",cex=2, xlab="Number of minutes ahead", ylab="Mean squared error for high dBZ pixel arrays",
     ylim=c(0,450))
points(seq(5,30,5),MSE2a,col="red",pch=3,cex=2)
lines(seq(5,30,5), MSE2a, col="red", lty=2)
legend("bottomright",legend=c("COTREC","the proposed approach"),
       col=c("darkblue","red"),pch=c(1,3),bty="n")
dev.off()
setwd(wd)

#mypath <- file.path( paste("../figures/predict_",Monsoon,"_",filenames[storm.id],"_MSE_b.png",sep="") )
#png(mypath,width =10, height = 4, units = "in", res=300)
setwd(wd.figure); png(file = paste(filenames[storm.id], "_MSE_b",".png",sep=""), width=10, height=4, units = "in", res=300)
par(las=TRUE)
par(mar=c(4,4,1,1))
par(mfrow=c(1,2))
plot(seq(5,30,5), MSE1b, col="darkblue", type="o",cex=2, xlab="Number of minutes ahead", ylab="Mean squared error",
     ylim=c(0,200))
points(seq(5,30,5),MSE2b,col="red",pch=3,cex=2)
lines(seq(5,30,5), MSE2b, col="red", lty=2)
legend("bottomright",legend=c("COTREC","the proposed approach"),
       col=c("darkblue","red"),pch=c(1,3),bty="n")
plot(seq(5,30,5), MSE1c, col="darkblue", type="o",cex=2, xlab="Number of minutes ahead", ylab="Mean squared error for high dBZ pixel arrays",
     ylim=c(0,400))
points(seq(5,30,5),MSE2c,col="red",pch=3,cex=2)
lines(seq(5,30,5), MSE2c, col="red", lty=2)
legend("bottomright",legend=c("COTREC","the proposed approach"),
       col=c("darkblue","red"),pch=c(1,3),bty="n")
dev.off()
setwd(wd)

