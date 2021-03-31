rm(list=ls())
library(quadprog)

xNumber = 30
x <- runif(xNumber, -1, 1)
epsilon <- rnorm(xNumber, 0, (0.65^2))
y <- -3*(x^6) + (3*x) + epsilon
plot(x,y)

n <- 100
h <- 0.25
mat <- matrix(c(x,y),nrow=length(x))
data<- as.data.frame( mat )

ndwat <- ksmooth(x, y, kernel = c("box", "normal"), bandwidth = 0.45,
        range.x = range(x),
        n.points = n)

lines(ndwat$x, ndwat$y ,lwd =2, col =" green")


diff <- ndwat$y[2:length(ndwat$x)] - ndwat$y[1:length(ndwat$x)-1]

p<-rep(1/n , n)
A<-ndwat$y * n


eta = 0
nonmono = which(diff<eta)
i=1;
repeat{
  
  ll <- (rep(1,length(nonmono))*nonmono)+i
  ll2 <- (rep(1,length(nonmono))*nonmono)-i
  nonmono <- c(nonmono , ll ,ll2 )
  nonmono <- unique(sort(nonmono, decreasing = FALSE))
  nonmono <- nonmono[nonmono>1]
  nonmono <- nonmono[nonmono<n]
  #print(nonmono)
  m = length(nonmono)
  #Dmat
  Dmat =  diag(1,n , n)
  Dmat <- Dmat * (n^2) *2
  
  #dvec
  dvec <- rep(0 , n)
  
  #Amat
  Amat <- matrix(c(0),nrow=n , ncol = 1)
  Amat[nonmono , 1] = 1
  Amat <- cbind(Amat, diag(1,n,n))
  Amat[nonmono,nonmono+1]=0
  Amat <- cbind(Amat, matrix(c(0),nrow=n , ncol = n))#diag(1,n,n)
  Amat[nonmono,nonmono+1+n]=1
  testMat <- diag(-1,n,n-1) * ndwat$y
  testMat2 <- diag(1,n,n) * ndwat$y
  testMat2 <- testMat2[,2:n]
  matyy <- testMat + testMat2
  Amat <- cbind(Amat,matyy)
  
  #bvec
  bvec <- sum(p[nonmono])
  tempP = rep(1/n , n)
  tempP[nonmono] = 0
  bvec<-c(bvec , tempP)
  tempP = rep(0 , n)
  bvec<-c(bvec , tempP)
  tempP = rep(eta,n-1)
  bvec<-c(bvec , tempP)
  quad <- 0
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
i

p <- quad$solution
newres <- A * p
lines(ndwat$x, newres ,lwd =1, col =" black")
newDiff = newres[2:length(ndwat$x)] - newres[1:length(ndwat$x)-1]
plot(ndwat$x,ndwat$y)
lines(ndwat$x, ndwat$y ,lwd =1, col =" green")
lines(ndwat$x,newres ,lwd =1, col =" red")

