
JAGSoadaData<-function(nbdadata){

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

  dataLength<-maxNoInd<-0
  for(i in 1:length(nbdadata)){
    dataLength<-dataLength+length(unique(nbdadata[[i]]@event.id))
    maxNoInd<-max(maxNoInd,length(unique(nbdadata[[i]]@id)))
  }

  stMetric_allDiffusions<-array(NA,dim=c(dataLength,maxNoInd,dim(stMetric)[3]))
  asocILVdata_allDiffusions<-array(NA,dim=c(dataLength,maxNoInd,dim(asocILVdata)[3]))
  intILVdata_allDiffusions<-array(NA,dim=c(dataLength,maxNoInd,dim(intILVdata)[3]))
  multiILVdata_allDiffusions<-array(NA,dim=c(dataLength,maxNoInd,dim(multiILVdata)[3]))
  randomEffectdata_allDiffusions<-array(NA,dim=c(dataLength,maxNoInd,dim(randomEffectdata)[3]))
  availabilityToLearn_allDiffusions<-array(NA,dim=c(dataLength,maxNoInd))
  status_allDiffusions<-array(NA,dim=c(dataLength,maxNoInd))
  offsetMatrix_allDiffusions<-array(NA,dim=c(dataLength,maxNoInd,4))

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

    }
  }

  randomEffectsLevels<-max(randomEffectdata_allDiffusions,na.rm=T)
  noEvents<-dim(status_allDiffusions)[1]

  #If there are unequal numbers of individuals in each diffusion then there will be NAs appearing in
  #the phantom slots for non-existent individuals. We need to replace all these with zeros to avoid errors
  for(i in 1:dim(stMetric_allDiffusions)[3]) stMetric_allDiffusions[,,i][is.na(availabilityToLearn_allDiffusions)]<-0
  for(i in 1:dim(asocILVdata_allDiffusions)[3]) asocILVdata_allDiffusions[,,i][is.na(availabilityToLearn_allDiffusions)]<-0
  for(i in 1:dim(intILVdata_allDiffusions)[3]) intILVdata_allDiffusions[,,i][is.na(availabilityToLearn_allDiffusions)]<-0
  for(i in 1:dim(multiILVdata_allDiffusions)[3]) multiILVdata_allDiffusions[,,i][is.na(availabilityToLearn_allDiffusions)]<-0
  for(i in 1:dim(randomEffectdata_allDiffusions)[3]) randomEffectdata_allDiffusions[,,i][is.na(availabilityToLearn_allDiffusions)]<-0
  for(i in 1:dim(offsetMatrix_allDiffusions)[3]) offsetMatrix_allDiffusions[,,i][is.na(availabilityToLearn_allDiffusions)]<-0
  status_allDiffusions[is.na(availabilityToLearn_allDiffusions)]<-0
  availabilityToLearn_allDiffusions[is.na(availabilityToLearn_allDiffusions)]<-0
  #Since the availabilityToLearn=0 for all these slots, the data is ignored when fitting the model


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
      offsetMatrix=offsetMatrix_allDiffusions

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

    imputationILVs=imputationILVs
  )
  }

  return(oada_jagsData)

}

JAGSoadaModel<-function(JAGSoadaDataIn,modelFileName,randomModel=T,upperS=1000, asocPriorVar=1000, intPriorVar=10000, multiPriorVar=10000, REhyperPriorUpper=10){

  noSParams<-JAGSoadaDataIn$noSParams
  noAsocParams<-JAGSoadaDataIn$noAsocParams
  noIntParams<-JAGSoadaDataIn$noIntParams
  noMultiParams<-JAGSoadaDataIn$noMultiParams
  noRandomEffects<-JAGSoadaDataIn$noRandomEffects
  randomEffectsLevels<-JAGSoadaDataIn$randomEffectsLevels
  noEvents<-JAGSoadaDataIn$noEvents
  maxNoInd<-JAGSoadaDataIn$maxNoInd

  sPriors<-paste("\n\ts[",1:noSParams,"]~dunif(0,",upperS,")",sep="")
  unscaledST<-paste("+s[",1:noSParams,"]*stMetric[j,k,",1:noSParams,"]",sep="",collapse = "")
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
  if(noRandomEffects==0|!randomModel){
    REPriors<-NULL
    sampleRandomEffectsFirst<-NULL
    sampleRandomEffects<-paste("\n\tre[1,j]<-0",sep="")
    multiLP_RE<-NULL
    #This just fills in a couple of entries for the random effects which is not used in the model anyway
    RElevelNumber="2"
  }else{
    REPriors<-paste(paste("\n\tsigma[",1:noRandomEffects,"]~dunif(0,",REhyperPriorUpper,")",sep=""),
                    paste("\n\ttau[",1:noRandomEffects,"]<-1/(sigma[",1:noRandomEffects,"]*sigma[",1:noRandomEffects,"])",sep=""))
    sampleRandomEffectsFirst<-paste("\n\tre[",1:noRandomEffects,",1]<-0",sep="")
    sampleRandomEffects<-paste("\n\tre[",1:noRandomEffects,",j]~dnorm(0,tau[",1:noRandomEffects,"])",sep="")
    multiLP_RE<-paste("+re[",1:noRandomEffects,",(randomEffectdata[j,k,",1:noRandomEffects,"]+1)]",sep="",collapse = "")
    RElevelNumber="(randomEffectsLevels+1)"
  }

  if(noMultiParams==0&noRandomEffects==0){
    multiLP<-"+0"
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
    ", sampleRandomEffectsFirst,
      "
    for(j in 2:",RElevelNumber,"){",sampleRandomEffects,"
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
        #A node used to get WAIC
        pYgivenThetaTemp[j,k]<-probs[j,k]*status[j,k]
      }
      status[j,1:maxNoInd]~ dmulti(probs[j,1:maxNoInd],1)
      #A node used to get WAIC (lppd and effective parameters)
      pYgivenTheta[j]<-sum(pYgivenThetaTemp[j,1:maxNoInd])
    }

    #Calculate propST
    #The log transformations are used here for numerical stability when rates are high, notably when gettting the prior for propST

    for(j in 1:noEvents){
      for(k in 1:maxNoInd){
        learnersRateTemp[j,k]<-relativeRate[j,k]*status[j,k]
      }
      learnersRate[j]<-sum(learnersRateTemp[j,1:maxNoInd])
      for(l in 1:noSParams){
          for(k in 1:maxNoInd){
             logLearnerSocialRateTemp[j,k,l]<-(log(s[l])+log(stMetric[j,k,l])+(intLP[j,k]+multiLP[j,k]))*status[j,k]
          }
          loglearnerSocialRate[j,l]<-sum(logLearnerSocialRateTemp[j,1:maxNoInd,l])
	  logProbST[j,l]<-loglearnerSocialRate[j,l]- log(learnersRate[j])
          probST[j,l]<-exp(logProbST[j,l])
      }
    }
    for(l in 1:noSParams){
      propST[l]<-sum(probST[1:noEvents,l])/noEvents
    }
  }
",fill=FALSE)
  sink()
}

#Needs updating to include changes above
JAGSoadaModel_priors<-function(JAGSoadaDataIn,modelFileName,upperS=1000, asocPriorVar=1000, intPriorVar=10000, multiPriorVar=10000, REhyperPriorUpper=10){

  noSParams<-JAGSoadaDataIn$noSParams
  noAsocParams<-JAGSoadaDataIn$noAsocParams
  noIntParams<-JAGSoadaDataIn$noIntParams
  noMultiParams<-JAGSoadaDataIn$noMultiParams
  noRandomEffects<-JAGSoadaDataIn$noRandomEffects
  randomEffectsLevels<-JAGSoadaDataIn$randomEffectsLevels
  noEvents<-JAGSoadaDataIn$noEvents
  maxNoInd<-JAGSoadaDataIn$maxNoInd

  sPriors<-paste("\n\ts[",1:noSParams,"]~dunif(0,",upperS,")",sep="")
  unscaledST<-paste("+s[",1:noSParams,"]*stMetric[j,k,",1:noSParams,"]",sep="",collapse = "")
  if(noAsocParams==0){
    asocPriors<-NULL
    asocialLP<-"0"
  }else{
    asocPriors<-paste("\n\tbetaAsoc[",1:noAsocParams,"]~dnorm(0,",1/asocPriorVar,")",sep="")
    asocialLP<-paste("+betaAsoc[",1:noAsocParams,"]*asocILVdata[j,k,",1:noAsocParams,"]",sep="",collapse = "")
  }
  if(noIntParams==0){
    intPriors<-NULL
    intLP<-"0"
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
    sampleRandomEffectsFirst<-paste("\n\tre[",1:noRandomEffects,",1]<-0)",sep="")
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

    ", sampleRandomEffectsFirst,
    "
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

# BLOCKED OUT IN ORDER TO SAMPLE THE PRIORS
#    for(j in 1:noEvents){
#      for(k in 1:maxNoInd){
#        probs[j,k]<-relativeRate[j,k]/sum(relativeRate[j,1:maxNoInd])
#      }
#      status[j,1:maxNoInd]~ dmulti(probs[j,1:maxNoInd],1)
#    }

    #Calculate propST
    #The log transformations are used here for numerical stability when rates are high, notably when gettting the prior for propST

    for(j in 1:noEvents){
      for(k in 1:maxNoInd){
        learnersRateTemp[j,k]<-relativeRate[j,k]*status[j,k]
      }
      learnersRate[j]<-sum(learnersRateTemp[j,1:maxNoInd])
      for(l in 1:noSParams){
          for(k in 1:maxNoInd){
             logLearnerSocialRateTemp[j,k,l]<-(log(s[l])+log(stMetric[j,k,l])+(intLP[j,k]+multiLP[j,k]))*status[j,k]
          }
          loglearnerSocialRate[j,l]<-sum(logLearnerSocialRateTemp[j,1:maxNoInd,l])
	  logProbST[j,l]<-loglearnerSocialRate[j,l]- log(learnersRate[j])
          probST[j,l]<-exp(logProbST[j,l])
      }
    }
    for(l in 1:noSParams){
      propST[l]<-sum(probST[1:noEvents,l])/noEvents
    }
  }
",fill=FALSE)
  sink()
}
