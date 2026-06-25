#setup the 2d dataset
rows = seq(1, nrow(mtcars), by = 2) #limiting the number of observations
df = mtcars[rows,c(1,6)]

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

#for now we are only finding the depth of a single point
#we start by shifting all the points so the target point 
#is equivalent to finding the depth of the origin
targetPoint = c(20.5, 3.1)
shiftedMPG = numeric()
shiftedWT = numeric()

#we also need to calculate the true angle that the point forms with the x-axis
thetaAngle = numeric()

for(i in 1:nrow(df)){
  shiftedMPG[i] = shiftX(targetPoint[1],df[i,1])
  shiftedWT[i] = shiftY(targetPoint[2],df[i,2])
  thetaAngle[i] = theta(shiftedMPG[[i]],shiftedWT[[i]])
}

#we calculate two angles for each observation
# 1 - the switch from projecting right to left
# 2 - the switch from projecting left to right
inFrom_thetaMinusPiby2 = numeric()
outAt_thetaPlusPiby2 = numeric()

for(i in 1:nrow(df)){
  inFrom_thetaMinusPiby2[i] = thetaAngle[i] - (pi/2)
  outAt_thetaPlusPiby2[i] = thetaAngle[i] + (pi/2)
}

#note that we will have some angles less than 0 and some greater than 2pi
#the line we rotate is only from 0 to 2pi
#thus we must fix those angles

fixRange = function(angle){
  fixedAngle = angle
  if(angle < 0){
    fixedAngle = (2*pi) + angle
  } else if(angle > (2*pi)){
    fixedAngle = angle - (2*pi)
  }
  return(fixedAngle)
}

for(i in 1:nrow(df)){
  inFrom_thetaMinusPiby2[i] = fixRange(inFrom_thetaMinusPiby2[i])
  outAt_thetaPlusPiby2[i] = fixRange(outAt_thetaPlusPiby2[i])
}

#For every angle in the 'in' list, we add (indicated by 1) one to our depth recorded so far
#and for every angle in the 'out' list, we subtract (indicated by 0) one from our depth recorded

actionForIn = numeric()
actionForOut = numeric()

for(i in 1:nrow(df)){
  actionForIn[i] = 1
  actionForOut[i] = 0
}

inFrom = cbind(inFrom_thetaMinusPiby2, actionForIn)
outAT = cbind(outAt_thetaPlusPiby2, actionForOut)

#combine both into one list, sort the list ascending order by angle
ListOfAngleSwitches = rbind(inFrom, outAT)
ListOfAngleSwitches = ListOfAngleSwitches[order(ListOfAngleSwitches[,1]),]

#function to check if we are in the right-half
checkQuad = function(x){
  if(x > 0){
    answer = 1
  } else{
    answer = 0
  }
  return(answer)
}

#calculating starting depth
minDepth = 0
right_count = 0
n = nrow(df)
i = 1
while(i <= nrow(df)){
  if(checkQuad(shiftedMPG[i]) == 1){
    right_count = right_count + 1
  }
  i = i + 1
}
left_count = n - right_count
minDepth = min(left_count,right_count)

univarDepth = numeric()
univarDepth[1] = minDepth

j = 1
while(j <= nrow(ListOfAngleSwitches)){
  if(ListOfAngleSwitches[j,2] == 0){ #one of the points can't project anymore
    right_count = right_count - 1
  } else{ #a new point can project now
    right_count = right_count + 1 
  }
  left_count = n - right_count
  univarDepth[j+1] = min(left_count, right_count)
  j = j + 1
}

slice_size = numeric()
slice_size[1] = ListOfAngleSwitches[1,1]
for(i in 2:nrow(ListOfAngleSwitches)){ 
  slice_size[i] = ListOfAngleSwitches[i,1] - ListOfAngleSwitches[i-1,1]
}

#update last slice size
slice_size[33] = (2*pi) - ListOfAngleSwitches[32,1]

slice_weight = slice_size/(2*pi)

irwDepth = sum(slice_weight*univarDepth)
irwDepth