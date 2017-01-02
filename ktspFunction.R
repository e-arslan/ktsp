ktsp<-function(training_data,training_grp,k=1){
  
  p<-dim(training_data)[1]
  c_0<-training_data[,training_grp==0]
  c_1<-training_data[,training_grp==1]
  if(!is.null(dim(c_0)[1])){
  Rank0 = apply(c_0,2,rank)
  }else{
  Rank0 = rank(c_0)
  }
  
  if(!is.null(dim(c_1)[1])){
    Rank1 = apply(c_1,2,rank)
  }else{
    Rank1 = rank(c_1)
  }
  
  kktsp<-c()
  kktsp$Indices = matrix(0,k,2)
  kktsp$FistTPScore = matrix(0,k,1)
  fscore<-matrix(0,p,p)
  sscore<-matrix(0,p,p)
  
  dat_w<-function(dat){
    p1<-dim(dat)[1]
    w<-matrix(0,p1,p1)
    newd <- t(dat)
    for(i in 1:p1) { 
      w[,i] <- colSums((newd - dat[i,]) > 0)
    } 
    return(w)
  }
  
  secscore  <-function(Rank0,Rank1,ind){
    
    if (dim(ind)[1]>1){
      sscore = matrix(0,dim(ind)[1])
      
      for(i in 1:dim(ind)[1]){
        muij_0 = sum(Rank0[ind[i,1],] - Rank0[ind[i,2],])/dim(Rank0)[2]
        muij_1 = sum(Rank1[ind[i,1],] - Rank1[ind[i,2],])/dim(Rank1)[2]
        sscore[i] = abs(muij_0 - muij_1)
      } 
      
    }else{
      muij_0 = sum(Rank0[ind[1],] - Rank0[ind[2],])/dim(Rank0)[2]
      muij_1 = sum(Rank1[ind[1],] - Rank1[ind[2],])/dim(Rank1)[2]
      sscore = abs(muij_0 - muij_1)
    }
    return(sscore)
  }
  
  w_0 = dat_w(c_0)
  w_1 = dat_w(c_1)
  pij_0 = w_0 / dim(c_0)[2] 
  pij_1 = w_1 / dim(c_1)[2]
  
  
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      fscore [i,j] = abs(pij_0[i,j]-pij_1[i,j]) 
    }
  }
  
  
  
  for (ii in 1:k){
    mscore = which(fscore==max(fscore),arr.ind = T)

    if (dim(mscore)[1] > 1) {
      sscore <- secscore(Rank0,Rank1,mscore)
      mscore<-matrix(mscore[which(sscore==max(sscore)),],nrow=1)
      #if we have multiple secondscore
      if(dim(mscore)[2]>2){
        mscore = matrix(mscore[,1:2],ncol=2)
      }
      
    }
    
    kktsp$Indices[ii,]<-mscore
    kktsp$FistTPScore[ii,1]<-fscore[mscore]
    fscore[mscore[,1],]<-fscore[mscore[,2],]<-fscore[,mscore[,1]]<-fscore[,mscore[,2]]<-0
  }
  
  
  # we want to be sure first column of indexes is bigger in class 0
  for (ii in 1:k){
    
    if(pij_0[kktsp$Indices[ii,1],kktsp$Indices[ii,2]]<pij_1[kktsp$Indices[ii,1],kktsp$Indices[ii,2]]){
      temp<-kktsp$Indices[ii,1]
      kktsp$Indices[ii,1]<-kktsp$Indices[ii,2]
      kktsp$Indices[ii,2]<-temp
    }
  }
  
  return(kktsp)
}


