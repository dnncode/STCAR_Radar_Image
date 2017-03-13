#******************************************************
#******************************************************
# Spatial kernel smoothing 
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



k.space.smooth1 <- function(grd, value, k.var=100, J=30){

  

  grd.sp <- grd
  grd <- coordinates(grd)
  grd <- grd/1000
  #case <- which(value !=0)
  #grd <- grd[case,]
  #value <- value[case]
  #grd.sp <- grd.sp[case]

  # select center
  cl <- kmeans(grd, J)
  center <- cl$centers
  #plot(grd,col="grey");points( cl$centers[,1], cl$centers[,2], col="red")

  #get F and diag(pi)
  n <- nrow( grd )
  Fi <- cbind( rep(1,n), grd )
  pi <- diag( dmvnorm(grd, mean = center[1,], sigma = k.var*diag(2), log = FALSE) ) 
  M <- pi %*% Fi
  for (j in 2:J){
    pi <- diag( dmvnorm(grd, mean = center[j,], sigma = k.var*diag(2), log = FALSE) ) 
    M <- cbind(M, pi %*% Fi)
  }
  #

  #least squares
  data.fit <- data.frame( cbind(value, M) )
  colnames(data.fit) <- c("y", paste("x",seq(1:ncol(M)),sep=""))
  fit <- lm(y~.-1, data=data.fit)

  #fit
  value.fit <- predict(fit,data.fit)

  if (TRUE){
    Max <- max(c(value, value.fit), na.rm=TRUE)
    Min <- min(c(value, value.fit), na.rm=TRUE)
    par(mfrow=c(1,2))
    value.spdf <- SpatialPixelsDataFrame(grd.sp,data.frame(value),proj4string=CRS(SVY21))
    image(value.spdf,col=colPalette,zlim=c(Min,Max))
    value.fit.spdf <- SpatialPixelsDataFrame(grd.sp,data.frame(value.fit),proj4string=CRS(SVY21))
    image(value.fit.spdf,col=colPalette,zlim=c(Min,Max))
  }

  y <- (data.fit$y - value.fit)
  output <- cbind(value.fit, y)

  OUTPUT <- list()
  OUTPUT[[1]] <- output
  OUTPUT[[2]] <- fit$coefficients
  return(OUTPUT)
}

#k.space.smooth(grd=grd.temp, value=g.value[case,i-1], k.var=100, J=30)





k.space.smooth2 <- function(grd, grd2,value, k.var=100, J=30, var_y){

  PD = is.positive.definite(var_y)

  grd.sp <- grd
  grd.sp.2 <- grd2
  grd <- coordinates(grd)
  grd <- grd/1000
  grd2 <- coordinates(grd2)
  grd2 <- grd2/1000
  #case <- which(value !=0)
  #grd <- grd[case,]
  #value <- value[case]
  #grd.sp <- grd.sp[case]

  # select center
  cl <- kmeans(grd, J)
  center <- cl$centers
  #plot(grd,col="grey");points( cl$centers[,1], cl$centers[,2], col="red")
  #plot(grd2,col="grey");points( cl$centers[,1], cl$centers[,2], col="red")


  #get F and diag(pi)
  n <- nrow( grd )
  Fi <- cbind( rep(1,n), grd )
  Fi_2 <- cbind( rep(1,n), grd2 )
  pi <- diag( dmvnorm(grd, mean = center[1,], sigma = k.var*diag(2), log = FALSE) ) 
  pi_2 <- diag( dmvnorm(grd2, mean = center[1,], sigma = k.var*diag(2), log = FALSE) ) 
  M <- pi %*% Fi
  M_2 <- pi_2 %*% Fi_2
  for (j in 2:J){
    pi <- diag( dmvnorm(grd, mean = center[j,], sigma = k.var*diag(2), log = FALSE) ) 
    pi_2 <- diag( dmvnorm(grd2, mean = center[j,], sigma = k.var*diag(2), log = FALSE) ) 
    M <- cbind(M, pi %*% Fi)
    M_2 <- cbind(M_2, pi_2 %*% Fi_2)
  }
  #

  #least squares
  data.fit <- data.frame( cbind(value, M) )
  colnames(data.fit) <- c("y", paste("x",seq(1:ncol(M)),sep=""))
  if (PD){
    fit <- lm.gls(y~.-1, data=data.fit, W=var_y,
                  na.action=na.exclude, inverse=TRUE, method = "qr")
  }else{
    fit <- lm(y~.-1, data=data.fit,
                  na.action=na.exclude)
  }
  
  #fit$coefficients

  #fit
  value.fit <- fit$fitted.values
  value.fit2 <- M_2 %*% matrix(coefficients(fit),ncol=1) # this is the fitted mu for the entire grid

  if (TRUE){
    par(mfrow=c(2,2))
    par(mar=c(4,4,1,1))
    Max <- max(c(value, value.fit), na.rm=TRUE)
    Min <- min(c(value, value.fit), na.rm=TRUE)
    value.spdf <- SpatialPixelsDataFrame(grd.sp,data.frame(value),proj4string=CRS(SVY21))
    image(value.spdf,col=colPalette,zlim=c(Min,Max))
    value.fit.spdf <- SpatialPixelsDataFrame(grd.sp,data.frame(value.fit),proj4string=CRS(SVY21))
    image(value.fit.spdf,col=colPalette,zlim=c(Min,Max))
    value.fit2.spdf <- SpatialPixelsDataFrame(grd.sp.2,data.frame(value.fit2),proj4string=CRS(SVY21))
    image(value.fit2.spdf,col=colPalette,zlim=c(Min,Max))
  }

  y <- (data.fit$y - value.fit)
  output <- cbind(value.fit, y, value.fit2)

  OUTPUT <- list()
  OUTPUT[[1]] <- output
  OUTPUT[[2]] <- fit$coefficients
  return(OUTPUT)
}

#k.space.smooth(grd=grd.temp, value=g.value[case,i-1], k.var=100, J=30)




