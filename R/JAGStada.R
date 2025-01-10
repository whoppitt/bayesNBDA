
JAGStadaData<-function(nbdadata){

  i<-1

  tempNBDAdata<-nbdadata[[i]]

  noSParams<-dim(tempNBDAdata@stMetric)[2]
  if(tempNBDAdata@asoc_ilv[1]=="ILVabsent"){noAsocParams<-0}else{noAsocParams<-dim(tempNBDAdata@asocILVdata)[2]}
  if(tempNBDAdata@int_ilv[1]=="ILVabsent"){noIntParams<-0}else{noIntParams<-dim(tempNBDAdata@intILVdata)[2]}
  if(tempNBDAdata@multi_ilv[1]=="ILVabsent"){noMultiParams<-0}else{noMultiParams<-dim(tempNBDAdata@multiILVdata)[2]}
  if(tempNBDAdata@random_effects[1]=="REabsent"){noRandomEffects<-0}else{noRandomEffects<-dim(tempNBDAdata@randomEffectdata)[2]}


  stMetric<-NULL
  asocILVdata<-NULL
  intILVdata<-NULL
  multiILVdata<-NULL
  randomEffectdata<-NULL
  offsetMatrix<-NULL
  status<-NULL
  TADAtime1<-NULL
  TADAtime2<-NULL

  for(i in 1:length(nbdadata)){

    tempNBDAdata<-nbdadata[[i]]


    stMetric<-rbind(stMetric,tempNBDAdata@stMetric)
    asocILVdata<-rbind(asocILVdata,tempNBDAdata@asocILVdata)
    intILVdata<-rbind(intILVdata,tempNBDAdata@intILVdata)
    multiILVdata<-rbind(multiILVdata,tempNBDAdata@multiILVdata)
    randomEffectdata<-rbind(randomEffectdata,tempNBDAdata@randomEffectdata)
    status<-c(status,tempNBDAdata@status)
    offsetMatrix<-rbind(offsetMatrix,tempNBDAdata@offsetCorrection)
    TADAtime1<-c(TADAtime1,tempNBDAdata@TADAtime1)
    TADAtime2<-c(TADAtime2,tempNBDAdata@TADAtime2)

  }

  randomEffectsLevels<-max(randomEffectdata[,1])
  noEvents<-length(status)
  linesOfData<-dim(stMetric)[1]
  zeros<-rep(0,linesOfData)

  imputedAsocILVs<-which(apply(is.na(asocILVdata),2,sum)>0)
  if(length(imputedAsocILVs)>0){
    imputationAsocILVs<-matrix(0,nrow=randomEffectsLevels+1,ncol=length(imputedAsocILVs))
    for(i in 2:(randomEffectsLevels+1)){
      for(j in 1:length(imputedAsocILVs)){
        imputationAsocILVs[i,j]<-asocILVdata[randomEffectdata[,1]==i,imputedIntILVs[j]][1]
      }
    }
    dimnames(imputationAsocILVs)[[2]]<-list(nbdadata[[1]]@int_ilv[imputedAsocILVs])
  }else{imputationAsocILVs<-NULL}


  imputedIntILVs<-which(apply(is.na(intILVdata),2,sum)>0)
  if(length(imputedIntILVs)>0){
    imputationIntILVs<-matrix(0,nrow=randomEffectsLevels+1,ncol=length(imputedIntILVs))
    for(i in 2:(randomEffectsLevels+1)){
      for(j in 1:length(imputedIntILVs)){
        imputationIntILVs[i,j]<-intILVdata[randomEffectdata[,1]==i,imputedIntILVs[j]][1]
      }
    }
    dimnames(imputationIntILVs)[[2]]<-list(nbdadata[[1]]@int_ilv[imputedIntILVs])
  }else{imputationIntILVs<-NULL}

  imputedMultiILVs<-which(apply(is.na(multiILVdata),2,sum)>0)
  if(length(imputedMultiILVs)>0){
    imputationMultiILVs<-matrix(0,nrow=randomEffectsLevels+1,ncol=length(imputedMultiILVs))
    for(i in 2:(randomEffectsLevels+1)){
      for(j in 1:length(imputedMultiILVs)){
        imputationMultiILVs[i,j]<-multiILVdata[randomEffectdata[,1]==i,imputedIntILVs[j]][1]
      }
    }
    dimnames(imputationMultiILVs)[[2]]<-list(nbdadata[[1]]@int_ilv[imputedMultiILVs])
  }else{imputationMultiILVs<-NULL}

  imputationILVs<-cbind(imputationAsocILVs,imputationIntILVs,imputationMultiILVs)

  tada_jagsData<-list(
    noEvents=noEvents,
    linesOfData=linesOfData,
    noSParams=noSParams,
    noAsocParams=noAsocParams,
    noIntParams=noIntParams,
    noMultiParams=noMultiParams,
    noRandomEffects=noRandomEffects,
    randomEffectsLevels=randomEffectsLevels,

    stMetric=stMetric,
    asocILVdata=asocILVdata,
    intILVdata=intILVdata,
    multiILVdata=multiILVdata,
    randomEffectdata=randomEffectdata,
    status=status,
    offsetMatrix=offsetMatrix,
    time1=TADAtime1,
    time2=TADAtime2,
    zeros=zeros,
    imputationILVs=imputationILVs
  )

  return(tada_jagsData)

}

JAGStadaModel<-function(JAGStadaDataIn,modelFileName,invRatePriorUpper,baseline="constant",sPriorFromRatePrior=F,
                        shapePrior=c(0,10),upperS=1000, asocPriorVar=1000, intPriorVar=10000, multiPriorVar=10000, REhyperPriorUpper=10){

  noSParams<-JAGStadaDataIn$noSParams
  noAsocParams<-JAGStadaDataIn$noAsocParams
  noIntParams<-JAGStadaDataIn$noIntParams
  noMultiParams<-JAGStadaDataIn$noMultiParams
  noRandomEffects<-JAGStadaDataIn$noRandomEffects
  randomEffectsLevels<-JAGStadaDataIn$randomEffectsLevels
  noEvents<-JAGStadaDataIn$noEvents
  maxNoInd<-JAGStadaDataIn$maxNoInd
  linesOfData<-JAGStadaDataIn$linesOfData
  time1<-JAGStadaDataIn$TADAtime1
  time2<-JAGStadaDataIn$TADAtime1

  if(baseline=="constant"){
    baselinePrior<-paste("\n\tRate~dunif(0,",1/invRatePriorUpper,")\n\tShape<-1",sep="",collapse="")
  }else{
    baselinePrior<-paste("\n\tRate~dunif(0,",1/invRatePriorUpper,")\n\tShape~dunif(",shapePrior[1],",",shapePrior[2],")",sep="",collapse="")
  }
  if(sPriorFromRatePrior){
    maxSTimesRatePrior<-paste("\n\tmaxSTimesRate[",1:noSParams,"]~dunif(0,",1/invRatePriorUpper,")",sep="",collapse="")
    sItself<-paste("\n\ts[",1:noSParams,"]<-maxSTimesRate[",1:noSParams,"]/(Rate*",apply(JAGStadaDataIn$stMetric,2,max),")",sep="",collapse="")
    sPriors<-paste(maxSTimesRatePrior,sItself)
  }else{
    sPriors<-paste("\n\ts[",1:noSParams,"]~dunif(0,",upperS,")",sep="")
  }
  unscaledST<-paste("s[",1:noSParams,"]*stMetric[j,",1:noSParams,"]",sep="",collapse = "+")
  if(noAsocParams==0){
    asocPriors<-NULL
    asocialLP<-"0"
  }else{
    asocPriors<-paste("\n\tbetaAsoc[",1:noAsocParams,"]~dnorm(0,",1/asocPriorVar,")",sep="")
    asocialLP<-paste("betaAsoc[",1:noAsocParams,"]*asocILVdata[j,",1:noAsocParams,"]",sep="",collapse = "+")
  }
  if(noIntParams==0){
    intPriors<-NULL
    intLP<-"0"
  }else{
    intPriors<-paste("\n\tbetaInt[",1:noIntParams,"]~dnorm(0,",1/intPriorVar,")",sep="")
    intLP<-paste0("betaInt[",1:noIntParams,"]*intILVdata[j,",1:noIntParams,"]",sep="",collapse = "+")
  }
  if(noMultiParams==0){
    multiPriors<-NULL
    multiLP<-NULL
  }else{
    multiPriors<-paste("\n\tbetaMulti[",1:noMultiParams,"]~dnorm(0,",1/multiPriorVar,")",sep="")
    multiLP<-paste("betaAsoc[",1:noMultiParams,"]*asocILVdata[j,",1:noMultiParams,"]",sep="",collapse = "+")
  }
  if(noRandomEffects==0){
    REPriors<-NULL
    sampleRandomEffects<-paste("\n\tre[1,j]<-0",sep="")
    multiLP_RE<-NULL
  }else{
    REPriors<-paste(paste("\n\tsigma[",1:noRandomEffects,"]~dunif(0,",REhyperPriorUpper,")",sep=""),
                    paste("\n\ttau[",1:noRandomEffects,"]<-1/(sigma[",1:noRandomEffects,"]*sigma[",1:noRandomEffects,"])",sep=""))
    sampleRandomEffects<-paste("\n\tre[",1:noRandomEffects,",j]~dnorm(0,tau[",1:noRandomEffects,"])",sep="")
    multiLP_RE<-paste("re[",1:noRandomEffects,",(randomEffectdata[j,",1:noRandomEffects,"]+1)]",sep="",collapse = "+")
  }

  if(noMultiParams==0&noRandomEffects==0){
    multiLP<-"0"
  }else{
    multiLP<-paste(multiLP,multiLP_RE,sep="",collapse = "+")
  }

  if(baseline=="constant"){
    solveHazards<-"solveHazards[j]<-Rate"
    cumulativeHazards<-"cumulativeHazards[j]<-(time1[j]-time2[j])*Rate"
  }
  if(baseline=="weibull"){
    solveHazards<-"solveHazards[j]<-dweib(time2[j],Shape,Rate)/(1-pweib(time2[j],Shape,Rate))"
    cumulativeHazards<-"cumulativeHazards[j]<-- log(1 - pweib(time1[j],Shape,Rate))+ log(1 - pweib(time2[j],Shape,Rate))"
  }
  if(baseline=="gamma"){
    solveHazards<-"solveHazards[j]<-dgamma(time2[j],Shape,Rate)/(1-pgamma(time2[j],Shape,Rate))"
    cumulativeHazards<-"cumulativeHazards[j]<-- log(1 - pgamma(time1[j],Shape,Rate))+ log(1 - pgamma(time2[j],Shape,Rate))"
  }

  #Specify the model in JAGS format (saves as a text file)
  sink(modelFileName)
  cat("

model{
    #1. Priors

    #Priors for baseline parameters",
      baselinePrior,"

    #Priors for S parameters",
      sPriors,
      "

    #Normal priors for ILV parameters
    #Effect of ILVs on asocial learning",asocPriors,"
    #Effect of ILVs on social learning",intPriors,"
    #Multiplicative ILVs (asocial effect = social effect)",multiPriors,"

    #Random effects", REPriors,"

    re[1,1]<-0
    for(j in 2:(randomEffectsLevels+1)){",sampleRandomEffects,"
    }


    for(j in 1:linesOfData){
      #Get the linear predictor for ILV effects on asocial and social learning
      asocialLP[j]<-offsetMatrix[j,2]+",asocialLP,"
      intLP[j]<-offsetMatrix[j,3]+",intLP,"
      multiLP[j]<-offsetMatrix[j,4]+",multiLP,"

      #Get the unscaled social transmission rate across all networks
      unscaledST[j]<-offsetMatrix[j,1]+",unscaledST,"

      #Get the relative rate of learning
      relativeRate[j]<-(exp(asocialLP[j]+multiLP[j])+exp(intLP[j]+multiLP[j])*unscaledST[j])

      #Get the component of log-likelihood for the solvers from relative rate
      lComp1[j] <- log(relativeRate[j])*status[j]

      #Account for baseline rate and calculate likelihoods- to go here
      ",solveHazards,"
      lComp2[j]<-log(solveHazards[j])*status[j]

      ",cumulativeHazards,"
      lComp3[j]<-relativeRate[j]*cumulativeHazards[j]

      logLik[j]<-lComp1[j]+lComp2[j]+lComp3[j]

      phi[j] <- -logLik[j] + 1
      zeros[j] ~ dpois(phi[j])
    }

    totalLogLik<-sum(logLik[1:linesOfData])

    #Add in a node to record pYgivenTheta to enable calculation of WAIC
    for(j in 1:linesOfData){
      #A node used to get WAIC (lppd and effective parameters)
      pYgivenTheta[j]<-exp(logLik[j])
    }

    #Calculate propST

    for(j in 1:linesOfData){
      for(l in 1:noSParams){
        socialRateTemp[j,l]<-s[l]*stMetric[j,l]*exp(intLP[j]+multiLP[j])
        probST[j,l]<-status[j]*socialRateTemp[j,l]/relativeRate[j]
      }
    }
    for(l in 1:noSParams){
      propST[l]<-sum(probST[1:linesOfData,l])/sum(status[1:linesOfData])
    }
  }
",fill=FALSE)
  sink()
}

JAGStadaInits<-function(JAGSoadaDataIn){
  noSParams<-JAGSoadaDataIn$noSParams
  noAsocParams<-JAGSoadaDataIn$noAsocParams
  noIntParams<-JAGSoadaDataIn$noIntParams
  noMultiParams<-JAGSoadaDataIn$noMultiParams
  noRandomEffects<-JAGSoadaDataIn$noRandomEffects

  outList<-NULL

  if(noSParams>0) outList<-c(outList,list(s=runif(noSParams,0,100)))
  if(noAsocParams>0) outList<-c(outList,list(betaAsoc=rnorm(noAsocParams,0.01)))
  if(noIntParams>0) outList<-c(outList,list(betaInt=rnorm(noIntParams,0.01)))
  if(noMultiParams>0) outList<-c(outList,list(betaMulti=rnorm(noMultiParams,0.01)))
  if(noRandomEffects>0) outList<-c(outList,list(sigma=runif(noRandomEffects,0,10)))

  return(outList)
}

