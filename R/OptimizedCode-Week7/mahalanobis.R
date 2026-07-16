#function to compute mahalanobis depth
mhlDepth = function(dataFrame, target){
  #Function that calculates mahalanobis distance
  mhlnbs_Dist = function(obs, mean, covMatrix){
    result = t(obs - mean) %*% 
      solve(covMatrix) %*% 
      (obs - mean)
    return(result)
  }

  df = dataFrame
  targetPoint = target
  S = cov(df)
  mu = t(t(colMeans(df)))
  
  distance = mhlnbs_Dist(targetPoint, mu, S)
  mhlnbsDepth_2d = 1 / (1 + distance)
  return(mhlnbsDepth_2d)
}