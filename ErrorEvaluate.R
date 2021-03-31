rm(list=ls())
library(quadprog)
library(splines)


# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}

# Function that returns Mean Absolute Error
nmpae <- function(error , maxFitted)
{
  mean(abs(error)/maxFitted)*100
}

naToMean <- function (y){
  for(i in 1:length(y)){
    if(is.na(y[i])){
      if(i == 1){
        y[1] <- 0
      }else{
        y[i]<-y[i-1]
      }
    }
  }
  return(y)
}
quadProgSolve <- function (nonmono, n, y,maxPowerValue){
  i=1
  repeat{
    
    ll <- (rep(1,length(nonmono))*nonmono)+i
    ll2 <- (rep(1,length(nonmono))*nonmono)-i
    nonmono <- c(nonmono , ll ,ll2 )
    nonmono <- unique(sort(nonmono, decreasing = FALSE))
    nonmono <- nonmono[nonmono>=1]
    nonmono <- nonmono[nonmono<=n]
    #print(nonmono)
    m = length(nonmono)
    #Dmat
    Dmat =  diag(1,n , n)
    Dmat <- Dmat * (n^2) *2
    
    #dvec
    dvec <- rep(0 , n)
    
    #Amat 
    Amat <- matrix(c(0),nrow=n , ncol = 1) #col1
    Amat[nonmono , 1] = 1
    Amat <- cbind(Amat, diag(1,n,n))#col 1->n+1 (n)
    Amat[nonmono,nonmono+1]=0
    tt11 <-diag(1,n,n)
    mono <- seq(1,n)
    tt11[which (mono %in% nonmono == FALSE),which (mono %in% nonmono == FALSE)] <- 0
    Amat <- cbind(Amat,tt11)#diag(1,n,n) 
    testMat <- diag(-1,n,n-1) * y
    testMat2 <- diag(1,n,n) * y
    testMat2 <- testMat2[,2:n]
    
    matyy <- testMat + testMat2
    
    Amat <- cbind(Amat,matyy)
    ttMat <- diag(-1,n,n) * (y*n)
    Amat <- cbind(Amat,ttMat)
    #bvec
    bvec <- sum(p[nonmono])
    tempP = rep(1/n , n)
    tempP[nonmono] = 0
    bvec<-c(bvec , tempP)
    tempP = rep(0 , n)
    bvec<-c(bvec , tempP)
    tempP = rep(eta,n-1)
    bvec<-c(bvec , tempP)
    bvec <- c(bvec ,rep((-maxPowerValue) , n))
    quad <- 0
    save(Amat, bvec, matyy,  file = "tempVars.RData")
    result = tryCatch({
      quad<-solve.QP(Dmat, dvec, Amat, bvec, meq=n+1, factorized = FALSE)
      break
    }, error = function(e) {
      
    })
    i<-i+1
    
    if(i>n-1){
      stop("not found")
    }
  }
  return(quad)
  
}

#debug(quadProgSolve)
xNumber = 200
load("C:/mehrdad/R workspace/cleanData22.RData")


## 75% of the sample size

cvNumbers<-10
## set the seed to make your partition reproductible
set.seed(77)

splits <-split(seq_len(nrow(newdata)), sample(1:cvNumbers, nrow(newdata), replace=T))

xlims =range(newdata$Speed)
x.grid=seq (from=xlims [1], to=xlims [2] , by=0.1)#length.out = ((xlims [2]-xlims [1])/0.1))

maxPowerValue <- max(newdata$Power)
par(mfcol=c(2,2))
lens <-c(0.2 , 0.5, 1, 1.5 , 2 , 2.5 , 3 , 5)
imroveBinsMae <- c()
imroveBinsRmse <- c()
imroveBinsNmpae <- c()

for(bandsW in lens){
  maeBefore <- c()
  maeAfter <- c()
  rmseBefore <- c()
  rmseAfter <- c()
  nmpaeBefore <- c()
  nmpaeAfter <- c()
  for(iterat in 1:cvNumbers){
  
    train <- newdata[- unlist(splits[iterat]),]
    test <- newdata[unlist(splits[iterat]),]
    
    x<-train$Speed
    y<-train$Power
    
    
    new <- data.frame(x = x.grid)
    n <- length(x.grid)#xNumber
    h <- 0.25
    data<-train
    

    #fit2=smooth.spline (x ,y ,df =12)
    plot(x,y)
    title(paste("Local Regression - Nadaraya Watson \n bandwidth = ",bandsW , " CV-" , iterat ))
    ndwat <- ksmooth(x, y, kernel = c("box", "normal"), bandwidth =bandsW,
                     range.x = range(x),
                     x.points=x.grid)
    ndwat$y <- naToMean(ndwat$y)
    model <- ndwat
      
      
      
    #pred=predict (fit2 ,new,se=T)
    pred<-ndwat
    maxFit <- max(pred$y)
    
    lines(model ,lwd =2, col =" blue")
    
    diff <- pred$y[2:xNumber] - pred$y[1:xNumber-1]#,1
    p<-rep(1/n , n)
    A<-pred$y * n
    errors <- c()
    for( t in 1:dim(test)[1]){
      ind <- which(x.grid == test[t,]$Speed)
      predictedVal <- pred$y[ind]
      err <- test[t,]$Power - predictedVal
      errors <- c(errors , err)
    }
     
    mae1 <- mae(errors)
    rmse1 <- rmse(errors)
    
    maeBefore <- c( maeBefore , mae1)
    rmseBefore <- c( rmseBefore , rmse1)
    nmpaeBefore <- c(nmpaeBefore , nmpae(errors , maxFit))
    eta = 0
    nonmono = which(diff<eta)
    i=1
    qd <- quadProgSolve(nonmono, n, pred$y, maxPowerValue)
    p <- qd$solution
    newres <- A * p
    maxFitted <- max(newres)
    errors <- c()
    
    for( t in 1:dim(test)[1]){
      ind <- which(x.grid == test[t,]$Speed)
      predictedVal <- newres[ind]
      err <- test[t,]$Power - predictedVal
      errors <- c(errors , err)
    }
    mae2 <- mae(errors)
    rmse2 <- rmse(errors)
    maeAfter <- c( maeAfter , mae2)
    rmseAfter <- c( rmseAfter , rmse2)
    nmpaeAfter<- c(nmpaeAfter , nmpae(errors , maxFitted))
    lines(pred$x, newres ,lwd =1, col =" red")
  }
  ImproveMae <- mean(maeAfter-maeBefore) #mean(maeAfter) - mean(maeBefore)
  ImproveRmse <- mean(rmseAfter - rmseBefore)
  ImproveNmpae <-mean(nmpaeAfter - nmpaeBefore)
  
  imroveBinsMae <- c(imroveBinsMae , ImproveMae)
  imroveBinsRmse <- c(imroveBinsRmse , ImproveRmse)
  imroveBinsNmpae <- c(imroveBinsNmpae , mean(nmpaeAfter - nmpaeBefore))
}






imroveBinsMae2 <- c()
imroveBinsRmse2 <- c()
imroveBinsNmpae2 <- c()

#ns
par(mfcol=c(2,2))
knotNumbers <- c(3,10,15,20)
for(knotNumbers in c(3,10,15,20) ){


#title(paste("Natural Spline ", knotNumbers ," knots"))
  xNumber = length(x.grid)
  maeBefore2 <- c()
  maeAfter2 <- c()
  rmseBefore2 <- c()
  rmseAfter2 <- c()
  nmpaeBefore2 <- c()
  nmpaeAfter2 <- c()
  for(iterat in 1:cvNumbers){
    train <- newdata[- unlist(splits[iterat]),]
    test <- newdata[unlist(splits[iterat]),]
    
    x<-train$Speed
    y<-train$Power
    
    
    new <- data.frame(x = x.grid)
    n <- length(x.grid)#xNumber
    h <- 0.25
    data<-train
    
    fit <- lm(y~ns(x ,df =knotNumbers) ,data=data)
    pred=predict (fit ,new,se=T)
    maxFit <- max(pred$fit)
    plot(x,y)
    title(paste("Natural Spline \n #bins = ",knotNumbers , " CV-" , iterat ))
    lines(x.grid , pred$fit ,col ="red",lwd =2)
    diff <- pred$fit[2:xNumber] - pred$fit[1:xNumber-1]#,1
    p<-rep(1/n , n)
    A<-pred$fit * n
    errors2 <- c()
    for( t in 1:dim(test)[1]){
      ind <- which(x.grid == test[t,]$Speed)
      predictedVal <- pred$fit[ind]
      err <- test[t,]$Power - predictedVal
      errors2 <- c(errors2 , err)
    }
    
    mae12 <- mae(errors2)
    rmse12 <- rmse(errors2)
    
    maeBefore2 <- c( maeBefore2 , mae12)
    rmseBefore2 <- c( rmseBefore2 , rmse12)
    nmpaeBefore2 <- c(nmpaeBefore2 , nmpae(errors2 , maxFit))
    
    eta = 0
    nonmono = which(diff<eta)
    i=1
    qd <- quadProgSolve(nonmono, n, pred$fit , maxPowerValue)
    p <- qd$solution
    newres <- A * p
    lines(x.grid, newres ,lwd =1, col =" green")
    maxFitted <- max(newres)
    errors <- c()
    
    for( t in 1:dim(test)[1]){
      ind <- which(x.grid == test[t,]$Speed)
      predictedVal <- newres[ind]
      err <- test[t,]$Power - predictedVal
      errors <- c(errors , err)
    }
    mae22 <- mae(errors)
    rmse22 <- rmse(errors)
    maeAfter2 <- c( maeAfter2 , mae22)
    rmseAfter2 <- c( rmseAfter2 , rmse22)
    nmpaeAfter2<- c(nmpaeAfter2 , nmpae(errors , maxFitted))
  
  }
  ImproveMae <- mean(maeAfter2-maeBefore2) #mean(maeAfter) - mean(maeBefore)
  ImproveRmse <- mean(rmseAfter2 - rmseBefore2)
  ImproveNmpae <-mean(nmpaeAfter2 - nmpaeBefore2)
  
  imroveBinsMae2 <- c(imroveBinsMae2 , ImproveMae)
  imroveBinsRmse2 <- c(imroveBinsRmse2 , ImproveRmse)
  imroveBinsNmpae2 <- c(imroveBinsNmpae2 , ImproveNmpae)

}

