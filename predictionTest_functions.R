#functions for testing predictions


#CRPS 
CPRS <- function(pred, trueobs){
  z <- as.numeric((trueobs - pred$mean) / pred$sd)
  scores <- pred$sd * (z * (2 * pnorm(z,0,1) - 1) + 
                         2 * dnrom(z,0,1) - 1/sqrt(pi))
  return(scores)
}


#IS95 (95% interval scores)
intscore <- function(x ,y, alpha = 0.05){
  hw <- -qnorm(alpha/2) * x$sd
  scores <- 2 * hw + (2/alpha) * (((x$mean - hw) - y) * (y < x$mean - hw) +
                                    (y - (x$mean +hw)) * (y > x$mean + hw))
  return(scores)
}


#coverage (coverage scores)
cvg <- function(x, y, alpha = 0.05){
  hw <- -qnorm(alpha/2) * x$sd
  scores <- y >= (x$mean - hw) & y <= (x$mean + hw)
  return(scores)
}