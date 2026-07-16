simplicialBaseline = function(dataFrame, target){
  #A function to calculate the distance between two points
  distance = function(point1, point2){
    return(sqrt(sum((point2-point1)^2)))
  }
  
  #formula to calculate the area of a triangle (heron's)
  triangleArea = function(a, b, c){
    s = (a + b + c)/2
    return(sqrt(s*(s-a)*(s-b)*(s-c)))
  }
  
  df = dataFrame
  targetPoint = target
  n = nrow(df)
  totalTriangles = choose(n,3)
  counter = 0
  
  #algorithm (for every possible triangle, check if the triangle contains the point)
  for(i in 1:(n-2)){
    for(j in (i+1):(n-1)){
      for(k in (j+1):n){
        #get the lengths of the sides of the triangle
        a = distance(df[i,],df[j,])
        b = distance(df[i,],df[k,])
        c = distance(df[j,],df[k,])
        #get the lengths of the sides formed by the targe point 
          #and the triangle vertices
        targetA = distance(df[i,],targetPoint)
        targetB = distance(df[j,],targetPoint)
        targetC = distance(df[k,],targetPoint)
        
        #compute all areas and check if the sum of the subtriangles 
          #equals the total
        totalArea = triangleArea(a, b, c)
        subArea1 = triangleArea(a,targetB,targetA)
        subArea2 = triangleArea(b,targetC,targetA)
        subArea3 = triangleArea(c,targetB,targetC)
        totalObserved = subArea1 + subArea2 + subArea3
        
        #using a tolerance instead of exact comparison
        if(abs(totalArea - totalObserved) < 1e-6){ 
          counter = counter + 1
        }
      }
    }
  }
  
  simplicialBase = counter/totalTriangles
  return(simplicialBase)
}

simplicialOpt = function(dataFrame, target){
  df = dataFrame
  targetPoint = target
  n = nrow(df)
  totalTriangles = choose(n,3)
  
  #shift to origin
  shiftedMPG = df[,1] - targetPoint[1]
  shiftedWT = df[,2] - targetPoint[2]
  
  #we also need to calculate the true angle that the point forms with the x-axis
  thetaAngle = atan2(shiftedWT, shiftedMPG)
  #fix negative angles for points below the x axis (i.e. negative y-coordinate)
  thetaAngle = ifelse(thetaAngle < 0, ((2*pi) + thetaAngle), thetaAngle)
  
  #As per Rousseeuw and Ruts Paper algorithm, now we compute the antipodal angles
  switchAngle = ifelse(between(thetaAngle, 0, pi), thetaAngle + pi, thetaAngle - pi)
  
  #implementing updating mechanism
  decrease = rep.int(-1, n)
  increase = rep.int(1, n)
  
  #combine both into one list, sort the list ascending order by angle
  ListOfAngleSwitches = matrix(c(thetaAngle, switchAngle, increase, decrease), nrow = 2*n, ncol = 2)
  ListOfAngleSwitches = ListOfAngleSwitches[order(ListOfAngleSwitches[,1]),]
  
  #we start at the smallest angle (first point in sweep)
  firstTheta = min(thetaAngle)
  
  #finding the corresponding antipodal angle
  antipodal1 = 0
  if(firstTheta >= 0 && firstTheta <= pi){
    antipodal1 = firstTheta + pi
  } else{
    antipodal1 = firstTheta - pi
  }
  
  #implementing the algorithm
  fArray = numeric(n)
  fArray[1] = 0
  
  #making sure we don't compare decimals using ==
  b1Index = which.min(abs(ListOfAngleSwitches[,1] - antipodal1))
  
  fArray[1] = sum((ListOfAngleSwitches[,2] == 1) & 
                    (ListOfAngleSwitches[,1] < antipodal1))
  
  globalTotal = fArray[1]
  k = 2 #keeping track of current index for F array
  
  #finding the range we traverse over
  if(b1Index == 1){
    temp = c((b1Index+1):nrow(ListOfAngleSwitches))
  }else if(b1Index == nrow(ListOfAngleSwitches)){
    temp = c(1:(b1Index-1))
  }else{
    temp = c((b1Index+1):nrow(ListOfAngleSwitches), 1:(b1Index-1))
  }
  
  #following the updating mechanism
  for(i in temp){ 
    if(ListOfAngleSwitches[i,2] == 1){
      globalTotal = globalTotal + 1
    }else{
      fArray[k] = globalTotal
      k = k+1
    }
  }
  
  #as per the algorithm, we remove the points that don't belong in the halfplane
  #to get the final halfspace counts
  right_count = fArray - (1:n)
  
  #calculating the final depth using combinatorial definition
  triangleCounter = sum(choose(right_count, 2))
  simplicialDepth = (totalTriangles - triangleCounter)/totalTriangles
  return(simplicialDepth)
}