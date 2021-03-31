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








# fitting some methods on the data
dimmatrix = dim(mydata)
data <- as.data.frame( mydata )
attach(data)

fit=lm(Power~poly(Speed ,6) ,data=data)
plot(Speed ,Power)
coef(summary (fit))
Speed.grid = seq (from=0, to=20 , length=200)
preds=predict (fit ,newdata =list(Speed=Speed.grid),se=TRUE)
se.bands=cbind(preds$fit +2* preds$se.fit ,preds$fit -2* preds$se.fit)
title (" Degree-6 Polynomial ",outer =T)
lines(Speed.grid ,preds$fit ,lwd =2, col =" blue")
matlines (Speed.grid ,se.bands ,lwd =1, col =" blue",lty =3)


fit2=smooth.spline (Speed ,Power ,df =16)
fit3=smooth.spline (Speed ,Power ,cv=TRUE)
fit3$df
plot(Speed ,Power ,col ="darkgrey")
title ("Smoothing Spline")
lines(fit2 ,col ="red ",lwd =2)
lines(fit3 ,col =" blue",lwd =2)
legend ("topright",legend =c("16 DF" ,"cv"),col=c("red","blue"),lty =1, lwd =2, cex =.8)

fit4 =gam(Power~bs(Speed,  df = NULL, knots = NULL, degree = 3, intercept = FALSE ,Boundary.knots = range(Speed)),data=data)
preds=predict (fit4 ,newdata =list(Speed=Speed.grid),se=TRUE)
plot(Speed ,Power)
lines(Speed.grid ,preds$fit ,lwd =2, col =" blue")

fit5 =gam(Power~ns(Speed,  df = 20, knots = NULL,intercept = FALSE ,Boundary.knots = range(Speed)),data=data)
preds=predict (fit5 ,newdata =list(Speed=Speed.grid),se=TRUE)
plot(Speed ,Power)
lines(Speed.grid ,preds$fit ,lwd =2, col =" blue")


