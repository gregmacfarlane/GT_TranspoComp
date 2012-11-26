# Function to execute a gravity model using the convergence method
#		given in Meyer and Miller, chapter 5, problem 6.
gravityModel <- function(productions, attractions, cost.matrix, b, TOLERANCE) {
	
	# Input Checks
	# ==================================================================
	
	# Check that productions and attractions are the same size
	if(length(productions) != length(attractions))
		stop("The production and attraction vectors do not match.")
	
	# Make sure b > 0
	if(b <= 0)
		stop("The cost-sensitivity coefficient b must be positive.")
	
	# Tolerance > 1
	if(TOLERANCE <= 0)
		stop("The tolerance must be positive.")

	
	# Initialize convergence conditions
	# ==================================================================
	# run until difference between Trips, Trips1 fall below tolerance 
	DIFFERENCE = Inf;
	TRIPS  = matrix(rep(NA, length(productions)*length(attractions)),
									nrow=length(productions), ncol=length(attractions))
	TRIPS1 = matrix(rep(0, length(productions)*length(attractions)),
									nrow=length(productions), ncol=length(attractions))
	ASTAR = attractions
	
	
	# Program Loop
	# ==================================================================
	for(iteration in 1:100){ # maximum 100 iterations
		# check convergence
		if(DIFFERENCE < TOLERANCE) break
	
		# Calculate Trips from i to j
		for(i in 1:length(productions)){
			for(j in 1:length(attractions)){
				TRIPS[i,j] = ( productions[i] * ASTAR[j] * ( cost.matrix[i,j] )^-b ) /
												( sum(ASTAR * ( cost.matrix[i,] )^-b ));
			}
		}
		
		# Update attractions
		ASTAR1 = ASTAR * ( attractions / colSums( TRIPS ))
		ASTAR = ASTAR1
	
		# Calculate convergence condition
		DIFFERENCE = sum( abs( TRIPS1 - TRIPS ))
		TRIPS1 = TRIPS
	
	}
	#print(paste("Iterations: ", iteration, sep=""))
	return(TRIPS)
}


# Create a cost matrix based on distance and transfer assumptions

costMatrix <- function(distance, transfer) {
	# Graph of distances between stations
	d.matrix <- matrix(c(Inf,1,3,4,2,4,3,3,4,
											 1,Inf,2,3,1,3,2,2,3,
											 3,2,Inf,1,1,3,2,2,3,
											 4,3,1,Inf,2,4,3,3,4,
											 2,1,1,2,Inf,2,1,1,2,
											 4,3,3,4,2,Inf,1,3,4,
											 3,2,2,3,1,1,Inf,2,3,
											 3,2,2,3,1,3,2,Inf,1,
											 4,3,3,4,2,4,3,1,Inf),
										 nrow=9, ncol=9, TRUE)
	# Matrix indicating whether a tranfer penalty should be applied
	t.matrix <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1,
											 0, 0, 0, 0, 0, 1, 1, 1, 1,
											 0, 0, 0, 0, 0, 1, 1, 1, 1,
											 0, 0, 0, 0, 0, 1, 1, 1, 1,
											 0, 0, 0, 0, 0, 0, 0, 0, 0,
											 1, 1, 1, 1, 0, 0, 0, 0, 0,
											 1, 1, 1, 1, 0, 0, 0, 0, 0,
											 1, 1, 1, 1, 0, 0, 0, 0, 0,
											 1, 1, 1, 1, 0, 0, 0, 0, 0),
										 nrow=9,ncol=9, TRUE)
	# return weighted sum of both matrices
	return(distance*d.matrix + transfer*t.matrix)
}

# Function for random cost matrix
costRandomMatrix <- function(meandist = -0.5, sddist =0.5, 
														 meantransfer = 10, sdtransfer = 3,
														 mintransfer = 2, maxtransfer = 15,
														 meanmultiplier = 2.5, sdmultiplier = 1,
														 meanspeed=20.3, sdspeed = 5.296) {
	
	#===== set random parameters ==========================================================
	# transfer time
	random.transfer <- matrix(rtnorm(81,meantransfer, sdtransfer, 
																	 mintransfer, maxtransfer), nrow=9, ncol=9)

	
	# time perception multiplier
	time.multiplier <- rtnorm(1, meanmultiplier, sdmultiplier, 1.1)
	
	# distance between stations
	rdist <- rlnorm(8, meanlog= meandist, sdlog=sddist)
	
	# speed of the system (can be uniform because distances are random)
	randomspeed <- rnorm(1, meanspeed, sdspeed)
	
	
	#====== Graph of Distances  ==============
	d.matrix <- matrix(c(Inf,rdist[1],sum(rdist[1:3]),sum(rdist[1:4]),sum(rdist[1:2]),
											 	sum(rdist[c(1:2,5:6)]),sum(rdist[c(1:2,6)]),sum(rdist[c(1:2,7)]),
											 	sum(rdist[c(1:2,7:8)]),
											 rdist[1],Inf,sum(rdist[2:3]),sum(rdist[2:4]),rdist[2],
											 	sum(rdist[c(2,5:6)]),sum(rdist[c(2,6)]),sum(rdist[c(2,7)]),
											 	sum(rdist[c(2,7:8)]),
											 sum(rdist[1:3]),sum(rdist[2:3]), Inf, rdist[4],rdist[3],
											 	sum(rdist[c(3,5:6)]),sum(rdist[c(3,6)]),sum(rdist[c(3,7)]),
											 	sum(rdist[c(3,7:8)]),
											 sum(rdist[1:4]),sum(rdist[2:4]), rdist[4], Inf, sum(rdist[3:4]),
											 	sum(rdist[c(3:4,5:6)]),sum(rdist[c(3:4,6)]),sum(rdist[c(3:4,7)]),
											 	sum(rdist[c(3:4,7:8)]),
											 sum(rdist[1:2]),rdist[2],rdist[3],sum(rdist[3:4]),Inf,
											 	sum(rdist[5:6]),rdist[6],rdist[7],sum(rdist[7:8]),
											 sum(rdist[c(1:2,5:6)]),sum(rdist[c(2,5:6)]),sum(rdist[c(3,5:6)]),
											 	sum(rdist[c(3:4,5:6)]),sum(rdist[5:6]),Inf,rdist[5],sum(rdist[5:7]),
											 	sum(rdist[5:8]),
											 sum(rdist[c(1:2,6)]),sum(rdist[c(2,6)]),sum(rdist[c(3,6)]),
											 	sum(rdist[c(3:4,6)]),sum(rdist[6]),rdist[5],Inf,sum(rdist[6:7]),
											 	sum(rdist[6:8]),
											 sum(rdist[c(1:2,7)]),sum(rdist[c(2,7)]),sum(rdist[c(3,7)]),
											 	sum(rdist[c(3:4,7)]),sum(rdist[7]),sum(rdist[5:7]), sum(rdist[6:7]), 
											 	Inf, rdist[8],
											 sum(rdist[c(1:2,7:8)]),sum(rdist[c(2,7:8)]),sum(rdist[c(3,7:8)]),
											 	sum(rdist[c(3:4,7:8)]),sum(rdist[7:8]),sum(rdist[5:8]), sum(rdist[6:8]),
											 	rdist[8], Inf),
										 nrow=9, ncol=9, TRUE)
	
	traveltime <- d.matrix / randomspeed * 60 
								# miles/(miles per hour) * minutes per hour = minutes
	
	# ========= Transfer Time ================================
	t.matrix <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1,
											 0, 0, 0, 0, 0, 1, 1, 1, 1,
											 0, 0, 0, 0, 0, 1, 1, 1, 1,
											 0, 0, 0, 0, 0, 1, 1, 1, 1,
											 0, 0, 0, 0, 0, 0, 0, 0, 0,
											 1, 1, 1, 1, 0, 0, 0, 0, 0,
											 1, 1, 1, 1, 0, 0, 0, 0, 0,
											 1, 1, 1, 1, 0, 0, 0, 0, 0,
											 1, 1, 1, 1, 0, 0, 0, 0, 0),
										 nrow=9,ncol=9, TRUE)
	
	transfertime <- random.transfer * t.matrix
	
	
	# return sum of both matrices
	return(traveltime + time.multiplier * transfertime)
}


# Get the transfer volumes from a OD matrix
output.cleaner <- function(TRIPS)   {
	
	ws <- (TRIPS[1,6] + TRIPS[1,7] + TRIPS[2,6] + TRIPS[2,7])
	wn <- (TRIPS[1,8] + TRIPS[1,9] + TRIPS[2,8] + TRIPS[2,9])
	es <- (TRIPS[3,6] + TRIPS[3,7] + TRIPS[4,6] + TRIPS[4,7])
	en <- (TRIPS[3,8] + TRIPS[3,9] + TRIPS[4,8] + TRIPS[4,9])
	sw <- (TRIPS[6,1] + TRIPS[6,2] + TRIPS[7,1] + TRIPS[7,2])
	se <- (TRIPS[6,3] + TRIPS[6,4] + TRIPS[7,3] + TRIPS[7,4])
	nw <- (TRIPS[8,1] + TRIPS[8,2] + TRIPS[9,1] + TRIPS[9,2])
	ne <- (TRIPS[8,3] + TRIPS[8,4] + TRIPS[9,3] + TRIPS[9,4])
	
	transfers <- c(ws, wn, es, en, sw, se, nw, ne)
	
	return(transfers)
	
}