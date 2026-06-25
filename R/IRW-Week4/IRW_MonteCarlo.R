#setup the 2d dataset
rows = seq(1, nrow(mtcars), by = 2) #limiting the number of observations
df = mtcars[rows,c(1,6)]
targetPoint = c(20.5, 3.1)

#number of directions
m = 5000
thetas = runif(m, 0, 2*pi) #generate uniformly randomly
univarDepth = numeric()

#function to calculate scalar projection
calcProjection = function(point, theta){
  return((point[1]*cos(theta))+(point[2]*sin(theta)))
}

#unit direction vector for each theta is given by <cos(theta), sin(theta)>

#for each theta, for each point, compute scalar projection of original first, 
#then the rest, for each increment the appropriate variable, 
#record the univariate depth
for(i in 1:length(thetas)){
  lessThan = 0
  greaterThan = 0 
  targetProjection = calcProjection(targetPoint, thetas[i])
  for (j in 1:nrow(df)) {
    projection = calcProjection(df[j,], thetas[i])
    if(projection < targetProjection){ 
      lessThan = lessThan + 1
    } else{
      greaterThan = greaterThan + 1
    }
  }
  univarDepth[i] = min(lessThan, greaterThan)
}

monteCarloIrwDepth = mean(univarDepth)
monteCarloIrwDepth