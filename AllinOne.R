rm(list = ls())

#packages
{
  library(C50)
  require("ktspair")
  library(e1071)
  library(pamr)
  library(switchBox)
  library(tspair)
  
}

#functions
{
  ktsp<-function(training_data,training_grp,krange=seq(3,9,by=2)){
    
    p<-dim(training_data)[1]
    c_0<-training_data[,training_grp==0]
    c_1<-training_data[,training_grp==1]
    
    if(length(krange) > 1){
      kmax = ifelse(max(krange) %% 2 == 1, max(krange),max(krange)-1)
      kmin = ifelse(min(krange) %% 2 == 1, min(krange),min(krange)+1)
      krange = seq(kmin,kmax,by=2)
    } else{
      kmax = krange
    }
    
    
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
    kktsp$Indices = matrix(0,kmax,2)
    kktsp$FistTPScore = matrix(0,kmax,1)
    fscore<-matrix(0,p,p)
    sscore<-matrix(0,p,p)
    
    dat_w<-function(dat){
      p1<-dim(dat)[1]
      w<-matrix(0,p1,p1)
      newd <- t(dat)
      for(i in 1:p1) { 
        w[,i] <- colSums((newd - dat[i,]) >= 0)
      } 
      return(w)
    }
    tscore <- function(dta,grp,classifier,k){
      p<-dim(training_data)[1]
      c_0<-dta[,grp==0]
      c_1<-dta[,grp==1]
      
      s1 = c_0[classifier$Indices[1:k, 1], ] > c_0[classifier$Indices[1:k,2],]
      s2 = c_0[classifier$Indices[1:k, 1], ] < c_0[classifier$Indices[1:k,2],]
      s1 = apply(as.matrix(s1), 2, sum)
      s2 = apply(as.matrix(s2), 2, sum)
      stat0 = s1-s2
      
      s1 = c_1[classifier$Indices[1:k, 1], ] > c_1[classifier$Indices[1:k,2],]
      s2 = c_1[classifier$Indices[1:k, 1], ] < c_1[classifier$Indices[1:k,2],]
      s1 = apply(as.matrix(s1), 2, sum)
      s2 = apply(as.matrix(s2), 2, sum)
      stat1 = s1-s2
      
      tt <- ( abs(mean(stat0) - mean(stat1)) /
                sqrt(var(stat1) + var(stat0) + 0.000000001) )
      
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
    
    for (ii in 1:kmax){
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
    for (ii in 1:kmax){
      
      if(pij_0[kktsp$Indices[ii,1],kktsp$Indices[ii,2]]<pij_1[kktsp$Indices[ii,1],kktsp$Indices[ii,2]]){
        temp<-kktsp$Indices[ii,1]
        kktsp$Indices[ii,1]<-kktsp$Indices[ii,2]
        kktsp$Indices[ii,2]<-temp
      }
    }
    
    if (length(krange) > 1){
      tt = 0
      for (u in 1:length(krange)){
        k = krange[u]
        t = tscore(training_data, training_grp,kktsp,k)
        if(t > tt){
          kopt = k
          tt = t
        }
        
        
      }
      kktsp$Indices = kktsp$Indices[1:kopt,]
      kktsp$FistTPScore = kktsp$FistTPScore[1:kopt,]
      
    }
    
    return(kktsp)
  }
  ktspSecondScore<-function(training_data,training_grp,k=1){
    
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
        fscore [i,j] = secscore(Rank0,Rank1,matrix(c(i,j),nrow=1)) 
      }
    }
    
    
    for (ii in 1:k){
      mscore = which(fscore==max(fscore),arr.ind = T)
      
      if (dim(mscore)[1] > 1) {
        sscore <- secscore(Rank0,Rank1,mscore)
        mscore<-matrix(mscore[which(sscore==max(sscore)),],nrow=1)
        #if we have multiple secondscore
        c("noluyolan!")
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
  RowVar <- function(x) {
    rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
  }
  make.consecutive.int <- function (y){
    intialwarn = getOption("warn")
    options(warn = -1)
    if (is.null(y)) {
      return(NULL)
    }
    if (!is.vector(y)) {
      y = as.vector(as.character(y))
    }
    out <- as.integer(as.factor(as.character(y))) - 1
    options(warn = intialwarn)
    return(out)
  }
  svmrfeFeatureRanking = function(x,y){
    n = ncol(x)
    
    survivingFeaturesIndexes = seq(1:n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    
    while(length(survivingFeaturesIndexes)>0){
      #train the support vector machine
      svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,
                     scale=F, type="C-classification", kernel="linear" )
      
      #compute the weight vector
      w = t(svmModel$coefs)%*%svmModel$SV
      
      #compute ranking criteria
      rankingCriteria = w * w
      
      #rank the features
      ranking = sort(rankingCriteria, index.return = TRUE)$ix
      
      #update feature ranked list
      featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
      rankedFeatureIndex = rankedFeatureIndex - 1
      
      #eliminate the feature with smallest ranking criterion
      (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
      
    }
    
    return (featureRankedList)
  } 
}

#prepare data
{
  rvar = RowVar(data)
  qt <- quantile(rvar,probs = 1-2000/dim(data)[1],na.rm = T)
  data <- data[apply( data , 1 , function(x) var(x) > qt[1] ),]
  grp = data.matrix(grp)
  grp = make.consecutive.int(grp)
  n = dim(data)[2]
  p = dim(data)[1]
  ind = 1:n
  rownames(data) = 1:p
  
}


  
#define parameters
{
testperc = 0.8
rep = 100
}
  
#classifiers
{
  testPredDT = matrix(0,rep)
  nfeatsDT = matrix(0,rep)
  kktspair39 = matrix(0,rep)
  corpredicktspair39 = matrix(0,rep) 
  testPredNB = matrix(0,rep) 
  corpredicpam = matrix(0,rep) 
  nfeatsPAM = matrix(0,rep) 
  testPredSVM10 = matrix(0,rep) 
  testPredSVM100 = matrix(0,rep) 
  testPredSw39 = matrix(0,rep)
  kSw39 = matrix(0,rep)
  testPred9TSP = matrix(0,rep)
  testPredS9TSP = matrix(0,rep)
  corpredictsp = matrix(0,rep)
  kbtsp9 = matrix(0,rep)
  
  
  u=1
  while(u!=(rep+1)){
        test_ind = sample(ind,round(n*testperc))
        training_ind = ind[-test_ind]
        training_data = data.matrix(data[,training_ind])
        training_grp = grp[training_ind]
        test_data = data.matrix(data[,test_ind])
        test_grp = grp[test_ind]
        
        while (is.na(unique(training_grp)[2]) | table(training_grp)[1] < 2 | table(training_grp)[2] < 2){
          test_ind = sample(ind,round(n*testperc))
          training_ind = ind[-test_ind]
          training_data = data.matrix(data[,training_ind])
          training_grp = grp[training_ind]
          test_data = data.matrix(data[,test_ind])
          test_grp = grp[test_ind]
        }
        
        
        #kbtsp
        #feats <- k_sel_feats(training_data,training_grp,N_Gibbs=N_Gibbs,N_burn=N_burn,a=a,k=9)
        #kbtsp9[u] = mean(pre_kbtsp(test_data,feats) == test_grp)
        
        
        #decision tree
        #model = C5.0(t(training_data),as.factor(training_grp))
        #testPredDT[u] <- mean (predict(model,t(test_data),type="class") == as.factor(test_grp))
        #nfeatsDT[u] = sum(C5imp(model)!=0)
        
        #ktspair w k = 3:9
        # ccv <- cv(training_data, training_grp,cross=3)
        # ktspa <- ktspcalc(training_data, training_grp, ccv$k)
        # kktspair39[u] = ktspa$k
        # corpredicktspair39[u]<-mean(predict(ktspa,data.matrix(test_data))==test_grp)
        
        #naive bayes
        #classifier = naiveBayes(t(training_data),as.factor(training_grp))
        #testPredNB[u] <- mean (predict(classifier,t(test_data)) == as.factor(test_grp))
        
        #pam
        # tr.data = list(x = training_data, y = training_grp,
        #                geneid=as.character(1:nrow(training_data)),
        #                genenames=paste("g",as.character(1:nrow(training_data)),sep=""))
        # te.data = list(x=test_data, y = test_grp)
        # tr.train = pamr.train(tr.data)
        # new.scales <- pamr.adaptthresh(tr.train)
        # tr.train2 <- pamr.train(tr.data, threshold.scale=new.scales,prior = c(0.5,0.5))
        # nfeatsPAM[u] = dim(pamr.listgenes(tr.train2, tr.data, threshold=new.scales))[1]
        # corpredicpam[u] = mean(pamr.predict(tr.train2, te.data$x, threshold=new.scales) == te.data$y)
        
        #SVM-RFE
         featureRankedList = svmrfeFeatureRanking(t(training_data),as.factor(training_grp))
         svmModel10 = svm(t(training_data[featureRankedList[1:10],]), as.factor(training_grp), cost = 0.1, kernel="linear" ) 
         svmModel100 = svm(t(training_data[featureRankedList[1:100],]), as.factor(training_grp), cost = 0.1, kernel="linear" ) 
         testPredSVM10[u] <- mean (predict(svmModel10,t(test_data[featureRankedList[1:10],])) == as.factor(test_grp))
         testPredSVM100[u] <- mean (predict(svmModel100,t(test_data[featureRankedList[1:100],])) == as.factor(test_grp))
        
        #switchbox39
        #classifier = SWAP.KTSP.Train(training_data, as.factor(training_grp),FilterFunc = NULL,krange = seq(3,9,by=2))
        #testPredSw39[u] <- mean(SWAP.KTSP.Classify(test_data, classifier) == test_grp)
        #kSw39[u] = dim(classifier$TSPs)[1]
        
        #ktsp defined func
        # classifier9 = ktsp(training_data,training_grp,k=9)
        # testPred9TSP[u] = mean(predictktsp(test_data, classifier9) == test_grp)
        # classifierS9 = ktspSecondScore(training_data,training_grp,k=9)
        # testPredS9TSP[u] = mean(predictktsp(test_data, classifierS9) == test_grp)
        
        #tspair
        #tsp <- tspcalc(training_data, training_grp)
        #corpredictsp[u]<-mean(predict(tsp,data.matrix(test_data))==test_grp)
        
        
        
        cat(u,"\n")
        u=u+1
  }
  
}
  


#results
{
  mean(testPredDT)
  mean(nfeatsDT)
  
  mean(corpredicktspair39)
  mean(kktspair39)
  
  mean(testPredNB)
  
  mean(corpredicpam)
  mean(nfeatsPAM)
  
  mean(testPredSVM10)
  mean(testPredSVM100) 
  
  mean(testPredSw39)
  mean(kSw39)
  
  mean(testPred9TSP)
  mean(testPredS9TSP)
  
  mean(corpredictsp)
}

#write results
{
  write.csv(testPredDT,file= "DT.csv" )
  
  write.csv(corpredicktspair39,file= "ktspair.csv" )
  
  write.csv(testPredNB,file= "NB.csv" )
  
  write.csv(corpredicpam,file= "PAM.csv" )
  
  write.csv(testPredSVM10,file= "SVM10.csv" )
  
  write.csv(testPredSVM100,file= "SVM100.csv" )
  
  write.csv(testPredSw39,file= "switk39.csv" )
  
  write.csv(testPred9TSP,file= "9tsp.csv" )
  
  write.csv(testPredS9TSP,file= "S9tsp.csv" )
  
  write.csv(corpredictsp,file= "tsp.csv" )
}





































