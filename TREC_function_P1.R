#******************************************************
#******************************************************
# TREC function 1 (Parallel Computing Version)
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




###########################################################
# TREC function (START)
###########################################################

TREC.parallel1 <- function(i, Data, Data2, N.TREC.ctr, Box.id, Grd.TREC){

# Data: data.frame radar scan
# Data2: data.frame radar scan 2
# N.TREC.ctr: number of TREC boxes

	temp1 <- Data[Box.id[,i],1]
	search.R <- 10000 # search radius, 5min --> 10km, 120km/hour (big enough)
	center.select <- Grd.TREC[i,]
	dist <- sqrt(  ((center.select)[1]-(Grd.TREC)[,1])^2 + ((center.select)[2]-(Grd.TREC)[,2])^2  )	
	case <- which( (dist < search.R)&(dist>0) )

	temp3 <- array(0/0,dim=c(length(case),1))
	for (j in 1:length(case)){
		temp2 <- Data2[Box.id[,case[j]],1]
		sigma <- sqrt( sum( (temp1-mean(temp1))^2 )/360 * sum( (temp2-mean(temp2))^2 )/360)
		temp3[j,1] <- sum( (temp1-mean(temp1)) * (temp2-mean(temp2)) )/ (361) / sigma
	}

if (sum(!is.na(temp3))==0){
	Cross_R <- 0/0
	box_select <- i
}else{
	Cross_R <- max(temp3,na.rm=TRUE)
	box_select <- case[which(temp3==Cross_R)][1]
}

return( c(Cross_R,box_select) )
}


###########################################################
# TREC Parallel function (END)
###########################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################