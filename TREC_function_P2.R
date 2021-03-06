#******************************************************
#******************************************************
# TREC function 2 (Parallel Computing Version)
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

TREC.parallel2 <- function(i, Grd.TREC, CRoss_R, TREC.U, TREC.V){

# Data: data.frame radar scan
# Data2: data.frame radar scan 2
# Grd.p: grid point
# N.TREC.ctr: number of TREC boxes
	TREC.U <- matrix(TREC.U, nrow=1)
	TREC.V <- matrix(TREC.V, nrow=1)

	dist <- 10000
	temp <- ((Grd.TREC[i,1]) - (Grd.TREC)[,1])^2 + ((Grd.TREC[i,2]) - (Grd.TREC)[,2])^2
	case <- (sqrt(temp)<dist)
	w <- 1-exp(-CRoss_R[case]^2/0.2)	
	w_total <- sum(w,na.rm=TRUE)
	w <- w/w_total

	judge <- sum(!is.na(TREC.U[case]+TREC.V[case]))
	if ( judge>0 ){
		TREC.u.s <- sum( TREC.U[case] * w, na.rm=TRUE)
		TREC.v.s <- sum( TREC.V[case] * w, na.rm=TRUE)
		temp <- ((Grd.TREC[i,1]) + TREC.u.s - (Grd.TREC)[,1])^2 + ((Grd.TREC[i,2]) + TREC.v.s- (Grd.TREC)[,2])^2
		box_select_revised <- which(temp == min(temp,na.rm=TRUE))[1]
		#marker <- 0
	}else{
		TREC.u.s <- 0/0
		TREC.v.s <- 0/0
		box_select_revised <- 0/0
		#marker <- 1
	}

return(  c(box_select_revised,TREC.u.s,TREC.v.s) )

} # TREC.parallel2

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