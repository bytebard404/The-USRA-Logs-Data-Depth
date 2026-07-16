irwPreSetup = function(dataFrame, target){
  df = dataFrame
  targetPoint = target

  #for now we are only finding the depth of a single point
  #we start by shifting all the points so the target point 
  #is equivalent to finding the depth of the origin
  
  #shift to origin
  shiftedMPG = df[,1] - targetPoint[1]
  shiftedWT = df[,2] - targetPoint[2]
  
  #we also need to calculate the true angle that the point forms with the x-axis
  #%% (2pi) allows to wrap the angle in the standard 0 to 2pi range
  thetaAngle = atan2(shiftedWT, shiftedMPG) %% (2*pi)
  
  #we calculate two angles for each observation
  # 1 - the switch from not being able to project to being able to project
  # 2 - vice versa
  inFrom_thetaMinusPiby2 = (thetaAngle - (pi/2)) %% (2*pi)
  outAt_thetaPlusPiby2   = (thetaAngle + (pi/2)) %% (2*pi)
  
  #For every angle in the 'in' list, we add (indicated by 1) one to our depth recorded so far
  #and for every angle in the 'out' list, we subtract (indicated by -1) one from our depth recorded
  
  actionForIn = rep(1, nrow(df))
  actionForOut = rep(-1, nrow(df))
  
  #create a data frame of all information and sort the angles for the sweep
  ListOfAngleSwitches = matrix(c(inFrom_thetaMinusPiby2, outAt_thetaPlusPiby2, 
                                 actionForIn, actionForOut), nrow = (2*nrow(df)), ncol = 2)
  ListOfAngleSwitches = ListOfAngleSwitches[order(ListOfAngleSwitches[,1]),]
  
  #calculating starting depth
  right_count_start = sum(shiftedMPG > 0)
  n = nrow(df)
  start_depth = min(right_count_start, n - right_count_start)

  #implementing the algorithm now
  #keeping a running tally of right_counts
  right_count = right_count_start + cumsum(ListOfAngleSwitches[,2])
  left_count = n - right_count
  
  #find the depth in each slice
  sweep_depths = pmin(left_count, right_count)
  univarDepth = c(start_depth, sweep_depths)
  
  #compute the size and weight of each slice
  temp1 = c(0,ListOfAngleSwitches[,1])
  temp2 = c(ListOfAngleSwitches[,1], (2*pi))
  slice_size = temp2 - temp1
  slice_weight = slice_size/(2*pi)
  
  #take the average to find the irw depth
  irwDepth = sum(slice_weight*univarDepth)
  return(irwDepth)
}

#-----------------------

irwPreSetupRandom = function(dataFrame, target){
  df = dataFrame
  targetPoint = target

  #for now we are only finding the depth of a single point
  #we start by shifting all the points so the target point 
  #is equivalent to finding the depth of the origin
  
  #shift to origin
  shiftedMPG = df[,1] - targetPoint[1]
  shiftedWT = df[,2] - targetPoint[2]
  
  #we also need to calculate the true angle that the point forms with the x-axis
  #%% (2pi) allows to wrap the angle in the standard 0 to 2pi range
  thetaAngle = atan2(shiftedWT, shiftedMPG) %% (2*pi)
  
  #we calculate two angles for each observation
  # 1 - the switch from not being able to project to being able to project
  # 2 - vice versa
  inFrom_thetaMinusPiby2 = (thetaAngle - (pi/2)) %% (2*pi)
  outAt_thetaPlusPiby2   = (thetaAngle + (pi/2)) %% (2*pi)
  
  #For every angle in the 'in' list, we add (indicated by 1) one to our depth recorded so far
  #and for every angle in the 'out' list, we subtract (indicated by -1) one from our depth recorded
  
  actionForIn = rep(1, nrow(df))
  actionForOut = rep(-1, nrow(df))
  
  #create a data frame of all information and sort the angles for the sweep
  ListOfAngleSwitches = matrix(c(inFrom_thetaMinusPiby2, outAt_thetaPlusPiby2, 
                                 actionForIn, actionForOut), nrow = (2*nrow(df)), ncol = 2)
  ListOfAngleSwitches = ListOfAngleSwitches[order(ListOfAngleSwitches[,1]),]
  
  #calculating starting depth
  right_count_start = sum(shiftedMPG > 0)
  n = nrow(df)
  start_depth = min(right_count_start, n - right_count_start)

  #implementing the algorithm now
  #keeping a running tally of right_counts
  right_count = right_count_start + cumsum(ListOfAngleSwitches[,2])
  left_count = n - right_count
  
  #find the depth in each slice
  sweep_depths = pmin(left_count, right_count)
  univarDepth = c(start_depth, sweep_depths)
  
  #upto this point we know what the depth in each slice is
  #in this iteration, we simply generate a random theta and
  #see what slice it lands in. 
  #remember we already know what the depth in that slice is
  
  angleSwitches = c(0,ListOfAngleSwitches[,1],2*pi)
  m = 5000
  thetas = runif(m, 0, 2*pi) #generate uniformly randomly
  univarDepthRand = numeric(m) #records the depth of each randomly generated theta
  
  #returns the index of the slice where each theta lands in
  j = findInterval(thetas, angleSwitches)
  univarDepthRand = univarDepth[j]
  
  PreSet_withRand_IrwDepth = mean(univarDepthRand)
  return(PreSet_withRand_IrwDepth)
}

#-------------------------------------

irwMonteCarlo = function(dataFrame, target){
  df = dataFrame
  df = as.matrix(df)
  targetPoint = target
  
  #number of directions
  m = 5000
  thetas = runif(m, 0, 2*pi) #generate uniformly randomly
  univarDepth = numeric(length(thetas)) #stores the depth associated with each randomly generated direction
  
  #function to calculate scalar projection
  #unit direction vector for each theta is given by <cos(theta), sin(theta)>
  calcProjection = function(point, theta){
    return((point[1]*cos(theta))+(point[2]*sin(theta)))
  }
  
  #for each theta, for each point, compute scalar projection of the target
  #and all the data points, find out the left and the right counts
  #record the univariate depth as the minimum of the two
  for(i in 1:length(thetas)){
    targetProjection = calcProjection(targetPoint, thetas[i])
    projection = calcProjection(df[,], thetas[i])
    lessThan = sum(projection < targetProjection)
    greaterThan = sum(projection >= targetProjection)
    univarDepth[i] = min(lessThan, greaterThan)
  }
  
  monteCarloIrwDepth = mean(univarDepth)
  return(monteCarloIrwDepth)
}