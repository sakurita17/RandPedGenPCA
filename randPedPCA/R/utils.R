
# Utility functions

# The oracle function that returns the value of A * G (but taking L^-1 as input)
oraculumLi <- function(Li, G, center=FALSE){
  if(center) G <- apply(G, 2, function(col) col - mean(col))
  Y <- spam::backsolve(t(Li), G)
  if(center) {
    return(apply(spam::forwardsolve(Li, Y), 2, function(col) col - mean(col)))
  } else {
    return(spam::forwardsolve(Li, Y))
  }
}
