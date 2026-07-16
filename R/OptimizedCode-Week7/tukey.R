tukeyDepthOpt = function(dataFrame, target){
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
  #implementing the algorithm now
  #iterate through the list of switch angles, 
  #adding or subtracting always keeping track of the minimum seen
  
  minDepth = 0
  #number of points projecting on the right initially
  right_count_start = sum(shiftedMPG > 0)
  #keeping a running tally of right_counts
  right_count = right_count_start + cumsum(ListOfAngleSwitches[,2])
  n = nrow(df)
  left_count = n - right_count
  start_depth = min(right_count_start, n - right_count_start)
  minDepth = min(start_depth, left_count, right_count)
  
  return(minDepth)
}