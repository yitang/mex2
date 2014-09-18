BigNumber <- 10^40
WeeNumber <- 10^(-10)
aLow <- -1 + 10^(-10)
obj <- function(para, y, y0, ...){
  BigNumber <- 10^40
  WeeNumber <- 10^(-10)
  aLow <- -1 + 10^(-10)
  a <- para[1]
  b <- para[2]
  z <- (y-y0*a)/y0^b
  m <- mean(z)
  s <- sd(z)
  mu <- y0*a + m * y0^b
  sig <- s * y0^b
  if(a < aLow[1] | s < WeeNumber | a > 1-WeeNumber  | b > 1-WeeNumber)
    return(1e40)
  res <- sum(0.5 * log(2*pi) + log(sig) + 0.5 * ((y - mu)/sig)^2)
  if (is.infinite(res))
    return(sign(res))
  return(res)
} 

HTMLE_obj <- function(para, yi, yj){
  BigNumber <- 10^40
  WeeNumber <- 10^(-10)
  aLow <- -1 + 10^(-10)
  a <- para[1]
  b <- para[2]

  z = (yj - yi * a ) / yi ^ b
  m <- mean(z)
  s <- sd(z)
       if(a < aLow[1] | s < WeeNumber | a > 1-WeeNumber  | b > 1-WeeNumber)
          return(1e40)
  
  nllh <- log(s * yi^b) + 1/2 * ((yj - a * yi - m * yi ^ b) / (s * yi ^ b))^2
  nllh <- sum(nllh)
  return(nllh)
}

HTMLE_obj2 <- function(para, yi, yj){
  BigNumber <- 10^40
  WeeNumber <- 10^(-10)
  aLow <- -1 + 10^(-10)
  a <- para[1]
  b <- para[2]
  
##  if (a < -1 | a > 1 | b >= 1)
##    return(1e4)
  z = (yj - yi * a ) / yi ^ b
  m <- mean(z)
  s <- sd(z)
      if(a < aLow[1] | s < WeeNumber | a > 1-WeeNumber  | b > 1-WeeNumber)
          return(1e40)
  
  llh <- 0.5 * log(2 * pi) + log(s * yi^b) + 1/2 * ((yj - a * yi - m * yi ^ b) / (s * yi ^ b))^2
llh <-  log(s * yi^b) + 1/2 * ((yj - a * yi - m * yi ^ b) / (s * yi ^ b))^2
  nllh <- -sum(llh)
  return(-nllh)
}
