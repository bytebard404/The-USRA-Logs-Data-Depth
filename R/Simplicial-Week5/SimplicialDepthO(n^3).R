#A function to calculate the distance between two points
distance = function(point1, point2){
  return(sqrt(((point2[1]-point1[1])^2)+((point2[2]-point1[2])^2)))
}

#formula to calculate the area of a triangle (heron's)
triangleArea = function(a, b, c){
  s = (a + b + c)/2
  return(sqrt(s*(s-a)*(s-b)*(s-c)))
}

#setup the 2d dataset
df = mtcars[,c(1,6)]
df = matrix(c(df[,1], df[,2]), byrow = FALSE, ncol = 2)
targetPoint = c(20.5, 3.1)
n = nrow(df)
totalTriangles = choose(n,3)
counter = 0

# algorithm (for every possible triangle, check if the triangle contains the point)
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

simplicialDepth = counter/totalTriangles
print("my algorithm calculated: ")
simplicialDepth

#verifying results
ddalphaDepth = depth.simplicial(targetPoint, df)
print("ddalpha calculated: ")
ddalphaDepth