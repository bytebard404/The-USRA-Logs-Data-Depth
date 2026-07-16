#setup the 2d dataset
x = rnorm(32,0,1)
y = rnorm(32,0,1)
df = cbind(x,y)

targetPoint = c(-0.4220978,0.820559)
n = nrow(df)
totalTriangles = choose(n,3)

#Since the line is rotated around the origin, we would need to shift every point
#the following functions help with the shifting

#function to shift x
shiftX = function(targetX, currX){
  answerX = currX - targetX
  return(answerX)
}

#function to shift y
shiftY = function(targetY, currY){
  answerY = currY - targetY
  return(answerY)
}

#function to calculate theta
theta = function(x, y){
  answer = atan2(y,x) 
  if((x < 0 && y < 0) || 
     (x > 0 && y < 0) ||
     (x == 0 && y < 0)){ #4th quad and 3rd quad and negative y-axis
    answer = (2*pi) + answer
  } 
  return(answer)
}

shiftedMPG = numeric()
shiftedWT = numeric()

#we also need to calculate the true angle that the point forms with the x-axis
thetaAngle = numeric()

for(i in 1:n){
  shiftedMPG[i] = shiftX(targetPoint[1],df[i,1])
  shiftedWT[i] = shiftY(targetPoint[2],df[i,2])
  thetaAngle[i] = theta(shiftedMPG[[i]],shiftedWT[[i]])
}

#As per Rousseeuw and Ruts Paper algorithm, now we compute the 
#antipodal angles
switchAngle = numeric()

for(i in 1:n){
  if(between(thetaAngle[i], 0, pi)){
    switchAngle[i] = thetaAngle[i] + pi
  } else{
    switchAngle[i] = thetaAngle[i] - pi
  }
}

#implementing updating mechanism
decrease = numeric()
increase = numeric()

for(i in 1:n){
  increase[i] = 1
  decrease[i] = -1
}

#attaching the correct reference with each angle
originalAngles = cbind(thetaAngle, increase)
antipodalAngles = cbind(switchAngle, decrease)

#combine both into one list, sort the list ascending order by angle
ListOfAngleSwitches = rbind(originalAngles, antipodalAngles)
ListOfAngleSwitches = ListOfAngleSwitches[order(ListOfAngleSwitches[,1]),]

first = 0 #finding first theta
j = 1
while(first != 1){
  if(ListOfAngleSwitches[j,2] == 1){
    first = 1
  } else{
    j = j + 1
  }
}

antipodal1 = 0 #getting first antipodal angle
if(between(ListOfAngleSwitches[j,1], 0, pi)){
  antipodal1 = ListOfAngleSwitches[j,1] + pi
} else{
  antipodal1 = ListOfAngleSwitches[j,1] - pi
}

#implementing the algorithm 
  #here fArray, b_i are as defined in the paper
fArray = numeric()
fArray[1] = 0
b1Index = 0

#getting F[1] - counting how many 
#input data points we see before antipodal1
for(i in 1:nrow(ListOfAngleSwitches)){
  if((ListOfAngleSwitches[i,2] == 1) && 
     (ListOfAngleSwitches[i,1] < antipodal1)){
    fArray[1] = fArray[1] + 1
  }
  
  if(ListOfAngleSwitches[i,1] == antipodal1){
    b1Index = i
  }
}

globalTotal = fArray[1]
k = 2 #keeping index track for F array

#defining the range where we continue our sweep
temp = c((b1Index+1):nrow(ListOfAngleSwitches), 1:(b1Index-1))
for(i in temp){ 
  if(ListOfAngleSwitches[i,2] == 1){
    globalTotal = globalTotal + 1
  }else{
    fArray[k] = globalTotal
    k = k+1
  }
}

right_count = numeric()
right_count[1] = 0
for(i in 1:length(fArray)){
  right_count[i] = fArray[i] - i
}

#calculating final depth as per combinatorial definition
triangleCounter = sum(choose(right_count, 2))
simplicialDepth = (totalTriangles - triangleCounter)/totalTriangles
print("my algorithm calculated: ")
simplicialDepth

#verifying results
library(ddalpha)
ddalphaDepth = depth.simplicial(targetPoint, df)
print("ddalpha calculated: ")
ddalphaDepth