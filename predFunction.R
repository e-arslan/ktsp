predictktsp <- function(testdta,classifier){
  n = dim(testdta)[2]
  k = dim(classifier$Indices)[1]
  estgrp = c()
  for(i in 1:n){
    c0 = 0
    c1 = 0
    for(j in 1:k){
      indi = classifier$Indices[j,]
      if(testdta[indi[1],i] > testdta[indi[2],i]){
        c0 = c0+1
      }else{
        c1 = c1+1
      }
      
    }
    
    if(c0 > c1){
      estgrp[i] = 0
    }else{
      estgrp[i] = 1
    }
    
  }
  return(estgrp)
}
