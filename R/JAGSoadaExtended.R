
JAGSoadaDataExtended<-function(nbdadata){

  i<-1

  tempNBDAdata<-nbdadata[[i]]

  noSParams<-dim(tempNBDAdata@stMetric)[2]
  if(tempNBDAdata@asoc_ilv[1]=="ILVabsent"){noAsocParams<-0}else{noAsocParams<-dim(tempNBDAdata@asocILVdata)[2]}
  if(tempNBDAdata@int_ilv[1]=="ILVabsent"){noIntParams<-0}else{noIntParams<-dim(tempNBDAdata@intILVdata)[2]}
  if(tempNBDAdata@multi_ilv[1]=="ILVabsent"){noMultiParams<-0}else{noMultiParams<-dim(tempNBDAdata@multiILVdata)[2]}
  if(tempNBDAdata@random_effects[1]=="REabsent"){noRandomEffects<-0}else{noRandomEffects<-dim(tempNBDAdata@randomEffectdata)[2]}

  dataTemplate<-matrix(0,ncol=length(unique(tempNBDAdata@id)),nrow=length(unique(tempNBDAdata@event.id)))
  dimnames(dataTemplate)[[1]]<-unique(tempNBDAdata@event.id)
  dimnames(dataTemplate)[[2]]<-unique(tempNBDAdata@id)
  availabilityToLearn<-status<-dataTemplate

  stMetric<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@stMetric)[2]))
  asocILVdata<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@asocILVdata)[2]))
  intILVdata<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@intILVdata)[2]))
  multiILVdata<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@multiILVdata)[2]))
  randomEffectdata<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@randomEffectdata)[2]))
  offsetMatrix<-array(0,dim=c(dim(dataTemplate),4))

  assMatrix<-tempNBDAdata@assMatrix
  assMatrixIndex<-tempNBDAdata@assMatrixIndex

  for(j in 1: length(tempNBDAdata@id)){
    index1<-which(tempNBDAdata@event.id[j]==unique(tempNBDAdata@event.id))
    index2<-which(tempNBDAdata@id[j]==unique(tempNBDAdata@id))

    stMetric[index1,index2,]<-tempNBDAdata@stMetric[j,]
    asocILVdata[index1,index2,]<-tempNBDAdata@asocILVdata[j,]
    intILVdata[index1,index2,]<-tempNBDAdata@intILVdata[j,]
    multiILVdata[index1,index2,]<-tempNBDAdata@multiILVdata[j,]
    randomEffectdata[index1,index2,]<-tempNBDAdata@randomEffectdata[j,]
    availabilityToLearn[index1,index2]<-1
    status[index1,index2]<-tempNBDAdata@status[j]
    offsetMatrix[index1,index2,]<-tempNBDAdata@offsetCorrection[j,]
  }

  #Get size of data across diffusions

  dataLength<-maxNoInd<-maxAssMatrixIndex<-0
  for(i in 1:length(nbdadata)){
    dataLength<-dataLength+length(unique(nbdadata[[i]]@event.id))
    maxNoInd<-max(maxNoInd,length(unique(nbdadata[[i]]@id)))
    maxAssMatrixIndex<-max(max(nbdadata[[i]]@assMatrixIndex),maxAssMatrixIndex)
  }

  stMetric_allDiffusions<-array(0,dim=c(dataLength,maxNoInd,dim(stMetric)[3]))
  asocILVdata_allDiffusions<-array(0,dim=c(dataLength,maxNoInd,dim(asocILVdata)[3]))
  intILVdata_allDiffusions<-array(0,dim=c(dataLength,maxNoInd,dim(intILVdata)[3]))
  multiILVdata_allDiffusions<-array(0,dim=c(dataLength,maxNoInd,dim(multiILVdata)[3]))
  randomEffectdata_allDiffusions<-array(0,dim=c(dataLength,maxNoInd,dim(randomEffectdata)[3]))
  availabilityToLearn_allDiffusions<-array(0,dim=c(dataLength,maxNoInd))
  status_allDiffusions<-array(0,dim=c(dataLength,maxNoInd))
  offsetMatrix_allDiffusions<-array(0,dim=c(dataLength,maxNoInd,4))

  assMatrixIndex_allDiffusions<-matrix(0,nrow=dataLength,ncol=2)
  assMatrix_allDiffusions<-array(0,dim=c(maxNoInd,maxNoInd,noSParams,maxAssMatrixIndex,length(nbdadata)))
  presenceMatrix_allDiffusions<-array(0,dim=c(dataLength,maxNoInd))
  statusMatrix_allDiffusions<-array(0,dim=c(dataLength,maxNoInd))

  i<-1
  index3<-1
  index4<-dim(status)[1]

  stMetric_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-stMetric
  asocILVdata_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-asocILVdata
  intILVdata_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-intILVdata
  multiILVdata_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-multiILVdata
  randomEffectdata_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-randomEffectdata
  availabilityToLearn_allDiffusions[index3:index4,1:dim(stMetric)[2]]<-availabilityToLearn
  status_allDiffusions[index3:index4,1:dim(stMetric)[2]]<-status
  offsetMatrix_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-offsetMatrix

  assMatrixIndex_allDiffusions[index3:index4,1]<-assMatrixIndex
  assMatrixIndex_allDiffusions[index3:index4,2]<-1
  assMatrix_allDiffusions[,,,,1]<-assMatrix[1:dim(stMetric)[2],1:dim(stMetric)[2],,]


  #The stMetric array has been corrected to allow for "ties" individuals who learn too close in time for us to assume they might have learned from one another
  #To allow for this when the assMatrix is used directly, here we set all network connections between such individuals to zero
  tieCount<-1
  if(tempNBDAdata@ties[tieCount]==1){
    tieTracker<-T
    tieCount2<-0
    while(tieTracker&(tieCount+tieCount2)<length(tempNBDAdata@ties)){
      tieCount2<-tieCount2+1
      if(tempNBDAdata@ties[tieCount+tieCount2]==1){
        assMatrix_allDiffusions[tempNBDAdata@orderAcq[tieCount:(tieCount+tieCount2)],tempNBDAdata@orderAcq[tieCount:(tieCount+tieCount2)],,i]<-0
      }else{
        tieTracker<-F
      }
    }
  }

  for(tieCount in 2:length(tempNBDAdata@ties)){
    if(tempNBDAdata@ties[tieCount]==1&tempNBDAdata@ties[tieCount-1]==0){
      tieTracker<-T
      tieCount2<-0
      while(tieTracker&(tieCount+tieCount2)<length(tempNBDAdata@ties)){
        tieCount2<-tieCount2+1
        if(tempNBDAdata@ties[tieCount+tieCount2]==1){
          assMatrix_allDiffusions[tempNBDAdata@orderAcq[tieCount:(tieCount+tieCount2)],tempNBDAdata@orderAcq[tieCount:(tieCount+tieCount2)],,i]<-0
        }else{
          tieTracker<-F
        }
      }
    }
  }

  #We also need the presence matrix and statusMatrix if the assMatrix is being used directly
  presenceMatrix_allDiffusions[index3:index4,1:dim(stMetric)[2]]<-tempNBDAdata@presenceMatrix
  statusMatrix_allDiffusions[index3:index4,1:dim(stMetric)[2]]<-tempNBDAdata@statusMatrix[,-dim(tempNBDAdata@statusMatrix)[2]]



  if(length(nbdadata)>1){

    for(i in 2:length(nbdadata)){

      tempNBDAdata<-nbdadata[[i]]

      dataTemplate<-matrix(0,ncol=length(unique(tempNBDAdata@id)),nrow=length(unique(tempNBDAdata@event.id)))
      dimnames(dataTemplate)[[1]]<-unique(tempNBDAdata@event.id)
      dimnames(dataTemplate)[[2]]<-unique(tempNBDAdata@id)
      availabilityToLearn<-status<-dataTemplate

      stMetric<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@stMetric)[2]))
      asocILVdata<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@asocILVdata)[2]))
      intILVdata<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@intILVdata)[2]))
      multiILVdata<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@multiILVdata)[2]))
      randomEffectdata<-array(0,dim=c(dim(dataTemplate),dim(tempNBDAdata@randomEffectdata)[2]))
      offsetMatrix<-array(0,dim=c(dim(dataTemplate),4))


      for(j in 1: length(tempNBDAdata@id)){
        index1<-which(tempNBDAdata@event.id[j]==unique(tempNBDAdata@event.id))
        index2<-which(tempNBDAdata@id[j]==unique(tempNBDAdata@id))

        stMetric[index1,index2,]<-tempNBDAdata@stMetric[j,]
        asocILVdata[index1,index2,]<-tempNBDAdata@asocILVdata[j,]
        intILVdata[index1,index2,]<-tempNBDAdata@intILVdata[j,]
        multiILVdata[index1,index2,]<-tempNBDAdata@multiILVdata[j,]
        randomEffectdata[index1,index2,]<-tempNBDAdata@randomEffectdata[j,]
        availabilityToLearn[index1,index2]<-1
        status[index1,index2]<-tempNBDAdata@status[j]
        offsetMatrix[index1,index2,]<-tempNBDAdata@offsetCorrection[j,]

      }

      index3<-index4+1
      index4<-index4+dim(status)[1]

      stMetric_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-stMetric
      asocILVdata_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-asocILVdata
      intILVdata_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-intILVdata
      multiILVdata_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-multiILVdata
      randomEffectdata_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-randomEffectdata
      availabilityToLearn_allDiffusions[index3:index4,1:dim(stMetric)[2]]<-availabilityToLearn
      status_allDiffusions[index3:index4,1:dim(stMetric)[2]]<-status
      offsetMatrix_allDiffusions[index3:index4,1:dim(stMetric)[2],]<-offsetMatrix

      assMatrixIndex_allDiffusions[index3:index4,1]<-tempNBDAdata@assMatrixIndex
      assMatrixIndex_allDiffusions[index3:index4,2]<-i
      assMatrix_allDiffusions[,,,,i]<-tempNBDAdata@assMatrix

      presenceMatrix_allDiffusions[index3:index4,1:dim(stMetric)[2]]<-tempNBDAdata@presenceMatrix
      statusMatrix_allDiffusions[index3:index4,1:dim(stMetric)[2]]<-tempNBDAdata@statusMatrix[,-dim(tempNBDAdata@statusMatrix)[2]]

      #The stMetric array has been corrected to allow for "ties": individuals who learn too close in time for us to assume they might have learned from one another
      #To allow for this when the assMatrix is used directly, here we set all network connections between such individuals to zero
      tieCount<-1
      if(tempNBDAdata@ties[tieCount]==1){
        tieTracker<-T
        tieCount2<-0
        while(tieTracker&(tieCount+tieCount2)<length(tempNBDAdata@ties)){
          tieCount2<-tieCount2+1
          if(tempNBDAdata@ties[tieCount+tieCount2]==1){
            assMatrix_allDiffusions[tempNBDAdata@orderAcq[tieCount:(tieCount+tieCount2)],tempNBDAdata@orderAcq[tieCount:(tieCount+tieCount2)],,i]<-0
          }else{
            tieTracker<-F
          }
        }
      }

      for(tieCount in 2:length(tempNBDAdata@ties)){
        if(tempNBDAdata@ties[tieCount]==1&tempNBDAdata@ties[tieCount-1]==0){
          tieTracker<-T
          tieCount2<-0
          while(tieTracker&(tieCount+tieCount2)<length(tempNBDAdata@ties)){
            tieCount2<-tieCount2+1
            if(tempNBDAdata@ties[tieCount+tieCount2]==1){
              assMatrix_allDiffusions[tempNBDAdata@orderAcq[tieCount:(tieCount+tieCount2)],tempNBDAdata@orderAcq[tieCount:(tieCount+tieCount2)],,i]<-0
            }else{
              tieTracker<-F
            }
          }
        }
      }
    }
  }


  randomEffectsLevels<-max(randomEffectdata_allDiffusions)
  noEvents<-dim(status_allDiffusions)[1]


  imputedAsocILVs<-which(apply(is.na(asocILVdata_allDiffusions),3,sum)>0)
  if(length(imputedAsocILVs)>0){
    imputationAsocILVs<-matrix(0,nrow=randomEffectsLevels+1,ncol=length(imputedAsocILVs))
    for(i in 2:(randomEffectsLevels+1)){
      for(j in 1:length(imputedAsocILVs)){
        imputationAsocILVs[i,j]<-asocILVdata_allDiffusions[1,i,imputedAsocILVs[j]]
      }
    }
    dimnames(imputationAsocILVs)[[2]]<-nbdadata[[1]]@int_ilv[imputedAsocILVs]
  }else{imputationAsocILVs<-NULL}


  imputedIntILVs<-which(apply(is.na(intILVdata_allDiffusions),3,sum)>0)
  if(length(imputedIntILVs)>0){
    imputationIntILVs<-matrix(0,nrow=randomEffectsLevels+1,ncol=length(imputedIntILVs))
    for(i in 2:(randomEffectsLevels+1)){
      for(j in 1:length(imputedIntILVs)){
        imputationIntILVs[i,j]<-intILVdata_allDiffusions[1,i,imputedIntILVs[j]]
      }
    }
    dimnames(imputationIntILVs)[[2]]<-nbdadata[[1]]@int_ilv[imputedIntILVs]
  }else{imputationIntILVs<-NULL}

  imputedMultiILVs<-which(apply(is.na(multiILVdata_allDiffusions),3,sum)>0)
  if(length(imputedMultiILVs)>0){
    imputationMultiILVs<-matrix(0,nrow=randomEffectsLevels+1,ncol=length(imputedMultiILVs))
    for(i in 2:(randomEffectsLevels+1)){
      for(j in 1:length(imputedMultiILVs)){
        imputationMultiILVs[i,j]<-multiILVdata_allDiffusions[1,i,imputedMultiILVs[j]]
      }
    }
    dimnames(imputationMultiILVs)[[2]]<-nbdadata[[1]]@int_ilv[imputedMultiILVs]
  }else{imputationMultiILVs<-NULL}

  imputationILVs<-cbind(imputationAsocILVs,imputationIntILVs,imputationMultiILVs)



  if(is.null(imputationILVs)){
    oada_jagsData<-list(
      noSParams=noSParams,
      noAsocParams=noAsocParams,
      noIntParams=noIntParams,
      noMultiParams=noMultiParams,
      noRandomEffects=noRandomEffects,
      randomEffectsLevels=randomEffectsLevels,
      noEvents=noEvents,
      maxNoInd=maxNoInd,

      stMetric=stMetric_allDiffusions,
      asocILVdata=asocILVdata_allDiffusions,
      intILVdata=intILVdata_allDiffusions,
      multiILVdata=multiILVdata_allDiffusions,
      randomEffectdata=randomEffectdata_allDiffusions,
      availabilityToLearn=availabilityToLearn_allDiffusions,
      status=status_allDiffusions,
      offsetMatrix=offsetMatrix_allDiffusions,

      assMatrixIndex=assMatrixIndex_allDiffusions,
      assMatrix=assMatrix_allDiffusions,
      presenceMatrix=presenceMatrix_allDiffusions,
      statusMatrix=statusMatrix_allDiffusions

    )
  }else{

    oada_jagsData<-list(
    noSParams=noSParams,
    noAsocParams=noAsocParams,
    noIntParams=noIntParams,
    noMultiParams=noMultiParams,
    noRandomEffects=noRandomEffects,
    randomEffectsLevels=randomEffectsLevels,
    noEvents=noEvents,
    maxNoInd=maxNoInd,

    stMetric=stMetric_allDiffusions,
    asocILVdata=asocILVdata_allDiffusions,
    intILVdata=intILVdata_allDiffusions,
    multiILVdata=multiILVdata_allDiffusions,
    randomEffectdata=randomEffectdata_allDiffusions,
    availabilityToLearn=availabilityToLearn_allDiffusions,
    status=status_allDiffusions,
    offsetMatrix=offsetMatrix_allDiffusions,

    assMatrixIndex=assMatrixIndex_allDiffusions,
    assMatrix=assMatrix_allDiffusions,
    presenceMatrix=presenceMatrix_allDiffusions,
    statusMatrix=statusMatrix_allDiffusions,

    imputationILVs=imputationILVs
  )
  }

  return(oada_jagsData)

}

JAGSoadaModelExtended<-function(JAGSoadaDataIn,modelFileName,upperS=1000, asocPriorVar=1000, intPriorVar=10000, multiPriorVar=10000, REhyperPriorUpper=10){

  noSParams<-JAGSoadaDataIn$noSParams
  noAsocParams<-JAGSoadaDataIn$noAsocParams
  noIntParams<-JAGSoadaDataIn$noIntParams
  noMultiParams<-JAGSoadaDataIn$noMultiParams
  noRandomEffects<-JAGSoadaDataIn$noRandomEffects
  randomEffectsLevels<-JAGSoadaDataIn$randomEffectsLevels
  noEvents<-JAGSoadaDataIn$noEvents
  maxNoInd<-JAGSoadaDataIn$maxNoInd

  sPriors<-paste("\n\ts[",1:noSParams,"]~dunif(0,",upperS,")",sep="")
  unscaledST<-paste("+s[",1:noSParams,"]*sum(assMatrix[k,1:maxNoInd,",1:noSParams,",assMatrixIndex[j,1],assMatrixIndex[j,2]]*statusMatrix[j,1:maxNoInd]*presenceMatrix[j,1:maxNoInd])",sep="",collapse = "")
  if(noAsocParams==0){
    asocPriors<-NULL
    asocialLP<-"+0"
  }else{
    asocPriors<-paste("\n\tbetaAsoc[",1:noAsocParams,"]~dnorm(0,",1/asocPriorVar,")",sep="")
    asocialLP<-paste("+betaAsoc[",1:noAsocParams,"]*asocILVdata[j,k,",1:noAsocParams,"]",sep="",collapse = "")
  }
  if(noIntParams==0){
    intPriors<-NULL
    intLP<-"+0"
  }else{
    intPriors<-paste("\n\tbetaInt[",1:noIntParams,"]~dnorm(0,",1/intPriorVar,")",sep="")
    intLP<-paste0("+betaInt[",1:noIntParams,"]*intILVdata[j,k,",1:noIntParams,"]",sep="",collapse = "")
  }
  if(noMultiParams==0){
    multiPriors<-NULL
    multiLP<-NULL
  }else{
    multiPriors<-paste("\n\tbetaMulti[",1:noMultiParams,"]~dnorm(0,",1/multiPriorVar,")",sep="")
    multiLP<-paste("+betaMulti[",1:noMultiParams,"]*asocILVdata[j,k,",1:noMultiParams,"]",sep="",collapse = "")
  }
  if(noRandomEffects==0){
    REPriors<-NULL
    sampleRandomEffects<-paste("\n\tre[1,j]<-0",sep="")
    multiLP_RE<-NULL
  }else{
    REPriors<-paste(paste("\n\tsigma[",1:noRandomEffects,"]~dunif(0,",REhyperPriorUpper,")",sep=""),
                    paste("\n\ttau[",1:noRandomEffects,"]<-1/(sigma[",1:noRandomEffects,"]*sigma[",1:noRandomEffects,"])",sep=""))
    sampleRandomEffects<-paste("\n\tre[",1:noRandomEffects,",j]~dnorm(0,tau[",1:noRandomEffects,"])",sep="")
    multiLP_RE<-paste("+re[",1:noRandomEffects,",(randomEffectdata[j,k,",1:noRandomEffects,"]+1)]",sep="",collapse = "")
  }

  if(noMultiParams==0&noRandomEffects==0){
    multiLP<-"0"
  }else{
    multiLP<-paste(multiLP,multiLP_RE,sep="",collapse = "")
  }




  #Specify the model in JAGS format (saves as a text file)
  sink(modelFileName)
  cat("
model{
    #1. Priors

    #Uniform prior for S parameters",
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


    for(j in 1:noEvents){
      for(k in 1:maxNoInd){
        #Get the linear predictor for ILV effects on asocial and social learning
        asocialLP[j,k]<-offsetMatrix[j,k,2]",asocialLP,"
        intLP[j,k]<-offsetMatrix[j,k,3]",intLP,"
        multiLP[j,k]<-offsetMatrix[j,k,4]",multiLP,"

        #Get the unscaled social transmission rate across all networks
        unscaledST[j,k]<-offsetMatrix[j,k,1]",unscaledST,"


        #Get the relative rate of learning

      relativeRate[j,k]<-(exp(asocialLP[j,k]+multiLP[j,k])+exp(intLP[j,k]+multiLP[j,k])*unscaledST[j,k])*availabilityToLearn[j,k]
      }

    }

    for(j in 1:noEvents){
      for(k in 1:maxNoInd){
        probs[j,k]<-relativeRate[j,k]/sum(relativeRate[j,1:maxNoInd])
      }
      status[j,1:maxNoInd]~ dmulti(probs[j,1:maxNoInd],1)
    }

    #Calculate propST

    for(j in 1:noEvents){
      for(k in 1:maxNoInd){
        learnersRateTemp[j,k]<-relativeRate[j,k]*status[j,k]
      }
      learnersRate[j]<-sum(learnersRateTemp[j,1:maxNoInd])
      for(l in 1:noSParams){
          for(k in 1:maxNoInd){
             learnerSocialRateTemp[j,k,l]<-s[l]*sum(assMatrix[k,1:maxNoInd,l,assMatrixIndex[j,1],assMatrixIndex[j,2]]*statusMatrix[j,1:maxNoInd]*presenceMatrix[j,1:maxNoInd])*exp(intLP[j,k]+multiLP[j,k])*status[j,k]
          }
          learnerSocialRate[j,l]<-sum(learnerSocialRateTemp[j,1:maxNoInd,l])
          probST[j,l]<-learnerSocialRate[j,l]/learnersRate[j]
      }
    }
    for(l in 1:noSParams){
      propST[l]<-sum(probST[1:noEvents,l])/noEvents
    }
  }
",fill=FALSE)
  sink()
}

