rm(list=ls())
library(readr)
library(splines)
library(gam)

mydata <- read_csv("C:/mehrdad/R workspace/turbine data/Turbin_2.csv", 
                   col_types = cols(Direction = col_skip(), 
                                    Status = col_skip()))

mydata <-mydata[complete.cases(mydata), ]
mydata <- mydata[mydata$Power >= 0, ]
mydata <- mydata[mydata$Speed > 0, ]
mydata <- mydata[mydata$Speed <= 25, ]

y<-which(mydata$Speed > 12)
u<-which(mydata$Power<1000)
inds <-intersect(y,u)
mydata<- mydata[-inds,]
y<-which(mydata$Speed < 5)
u<-which(mydata$Power>1000)
inds <-intersect(y,u)
mydata<- mydata[-inds,]

Speed.grid = sort(unique(mydata$Speed))
histData <- hist(mydata$Speed, breaks = Speed.grid, freq = TRUE)

meanValues <- c(0)
varValues <- c(0)
ind <- 1
for (sp in Speed.grid ){
  powerss <- mydata[mydata$Speed == sp, 2]$Power
  meanValues[ind] <- mean(powerss)
  varValues[ind] <-var(powerss)
  ind <- ind + 1
}
totalVar <- mean(varValues[complete.cases(varValues)])


ind <- 1
pointsToRemoveSpeed <- c()
pointsToRemovePower <- c()

for (sp in Speed.grid ){
  if(ind<length(Speed.grid)){
    
    for(pointInd in mydata[mydata$Speed == sp, 2]$Power){
      if((!is.na(varValues[ind]))&&((pointInd - meanValues[ind])^2)> 3 * totalVar ){#varValues[ind]
        #remove that point
        pointsToRemoveSpeed <- c(pointsToRemoveSpeed , sp)
        pointsToRemovePower <- c(pointsToRemovePower , pointInd)
      }
    }
  }
  ind<-ind + 1
}
removed = c(mydata[1,])
combined <- cbind(pointsToRemoveSpeed,pointsToRemovePower)
for( i in 1:dim(combined)[1]){#dim(combined)[1]
  index1 <- which(mydata$Speed== combined[i,1] )
  index2 <- which(mydata[index1,]$Power==combined[i,2])
  for(j in index2){
    removed <- rbind(removed , mydata[ index1[j],])
    mydata<-mydata[- index1[j],]
  }
}
plot ( mydata$Speed , mydata$Power)
plot(removed$Speed , removed$Power) 
points(Speed.grid , meanValues, col="red")



meanValues2 <- c(0)
ind2 <- 1
for (sp in Speed.grid ){
  powerss <- mydata[mydata$Speed == sp, 2]$Power
  meanValues2[ind2] <- mean(powerss)
  ind2 <- ind2 + 1
}
points(Speed.grid , meanValues2, col="red")

plot ( mydata$Speed , mydata$Power)
points(Speed.grid , meanValues2, col="red")
