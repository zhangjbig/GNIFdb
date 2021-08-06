#!/usr/bin/Rscript

args <- commandArgs()

i <- which(args=="-i")
input <- read.csv(args[i+1])
WT.pep <-input[,1]
MT.pep <- input[,2]
j <- which(args=="-o")
path <- args[j+1]

print(args)

WT.pep <- c(WT.pep,"VALPINPTI")
MT.pep <- c(MT.pep,"VVLPINPTI")
print(WT.pep)
print(MT.pep)
print(path)
load(paste0(path,"featurenames.Rdata"))
Peptype<-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')

PF_calculation <- function(MT.pep,WT.pep){
  
  library(Peptides)
  library(doParallel)
  data(AAdata)
  library(entropy)
  
  IDHwt<-cbind(MT.pep,WT.pep)
  IDHwt<-data.frame(IDHwt)
  colnames(IDHwt)<-c("MT.pep","WT.pep")
  #   entropy
  # 1.2.1 sequence entropy
  entropyseq<-function(pep){
    Entropy<-NULL
    for (i in 1:length(pep)){
      pep1<-substring(pep[i],1:nchar(as.character(pep[i])),1:nchar(as.character(pep[i])))
      Density<-NULL
      for (j in 1:length(unique(pep1))){
        posibility<-which(pep1==unique(pep1)[j])
        Density[j]<-length(posibility)/nchar(as.character(pep[i]))
      }
      Entropy[i]<-entropy(Density,unit = "log2")
    }
    return(Entropy)
  }
  #  1.2.2 residue entropy
  entropyresidue<-function(pep){
    entropy_residue<-matrix(nrow=length(MT.pep),ncol=20)
    for (i in 1:length(pep)){
      pep1<-substring(pep[i],1:nchar(as.character(pep[i])),1:nchar(as.character(pep[i])))
      for (j in 1:length(Peptype)){
        posibility<-which(pep1==Peptype[j])
        if (length(posibility)==0){
          entropy_residue[i,j]<-0
        }else{
          entropy_residue[i,j]<-(-1)*(length(posibility)/nchar(as.character(pep[i])))*log2(length(posibility)/nchar(as.character(pep[i])))
          
        }
      }
    }
    return(entropy_residue)
  }
  
  # 1.2.3  mutated position determined
  MTposition<-NULL
  for (i in 1:length(WT.pep)){
    MTposition[i]<-which(substring(WT.pep[i],1:9,1:9)!=substring(MT.pep[i],1:9,1:9))
  }
  IDHwt<-data.frame(IDHwt,MTposition)
  
  ###### 1.2.4 R package
  ##(1)
  model_process = function(n){
    c(
      hydrophobicity(MT.pep[n]),
      instaIndex(MT.pep[n]),
      charge(MT.pep[n]),
      mw(MT.pep[n]),
      boman(MT.pep[n]),
      aIndex(MT.pep[n]),
      autoCorrelation(MT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,center = TRUE),
      autoCovariance(MT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,center = TRUE),
      hmoment(MT.pep[n]),
      crossCovariance(MT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,property2 = AAdata$Hydrophobicity$Eisenberg,center = TRUE),
      
      hydrophobicity(substring(WT.pep[n],MTposition[n],MTposition[n])),
      charge(substring(WT.pep[n],MTposition[n],MTposition[n])),
      mw(substring(WT.pep[n],MTposition[n],MTposition[n])),
      boman(substring(WT.pep[n],MTposition[n],MTposition[n])),
      aIndex(substring(WT.pep[n],MTposition[n],MTposition[n])),
      hmoment(substr(MT.pep[n],MTposition[n],MTposition[n])),
      
      ifelse(MTposition[n]==1,hydrophobicity(substr(MT.pep[n],1,2)),hydrophobicity(substr(MT.pep[n],MTposition[n]-1,MTposition[n]))),
      ifelse(MTposition[n]==1,charge(substr(MT.pep[n],1,2)),charge(substr(MT.pep[n],MTposition[n]-1,MTposition[n]))),
      ifelse(MTposition[n]==1,mw(substr(MT.pep[n],1,2)),mw(substr(MT.pep[n],MTposition[n]-1,MTposition[n]))),
      ifelse(MTposition[n]==1,boman(substr(MT.pep[n],1,2)),boman(substr(MT.pep[n],MTposition[n]-1,MTposition[n]))),
      ifelse(MTposition[n]==1,aIndex(substr(MT.pep[n],1,2)),aIndex(substr(MT.pep[n],MTposition[n]-1,MTposition[n]))),
      ifelse(MTposition[n]==1,hmoment(substr(MT.pep[n],1,2)),hmoment(substr(MT.pep[n],MTposition[n]-1,MTposition[n]))),
      ifelse(MTposition[n]==9,hydrophobicity(substr(MT.pep[n],8,9)),hydrophobicity(substr(MT.pep[n],MTposition[n],MTposition[n]+1))),
      ifelse(MTposition[n]==9,charge(substr(MT.pep[n],8,9)),charge(substr(MT.pep[n],MTposition[n],MTposition[n]+1))),
      ifelse(MTposition[n]==9,mw(substr(MT.pep[n],8,9)),mw(substr(MT.pep[n],MTposition[n],MTposition[n]+1))),
      ifelse(MTposition[n]==9,boman(substr(MT.pep[n],8,9)),boman(substr(MT.pep[n],MTposition[n],MTposition[n]+1))),
      ifelse(MTposition[n]==9,aIndex(substr(MT.pep[n],8,9)),aIndex(substr(MT.pep[n],MTposition[n],MTposition[n]+1))),
      ifelse(MTposition[n]==9,hmoment(substr(MT.pep[n],8,9)),hmoment(substr(MT.pep[n],MTposition[n],MTposition[n]+1))),
      ifelse(MTposition[n]==1|MTposition[n]==2,hydrophobicity(substr(MT.pep[n],1,3)),hydrophobicity(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==1|MTposition[n]==2,instaIndex(substr(MT.pep[n],1,3)),instaIndex(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==1|MTposition[n]==2,charge(substr(MT.pep[n],1,3)),charge(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==1|MTposition[n]==2,mw(substr(MT.pep[n],1,3)),mw(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==1|MTposition[n]==2,boman(substr(MT.pep[n],1,3)),boman(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==1|MTposition[n]==2,aIndex(substr(MT.pep[n],1,3)),aIndex(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==1|MTposition[n]==2,hmoment(substr(MT.pep[n],1,3)),hmoment(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==9|MTposition[n]==8,hydrophobicity(substr(MT.pep[n],7,9)),hydrophobicity(substr(MT.pep[n],MTposition[n],MTposition[n]+2))),
      ifelse(MTposition[n]==9|MTposition[n]==8,instaIndex(substr(MT.pep[n],7,9)),instaIndex(substr(MT.pep[n],MTposition[n],MTposition[n]+2))),
      ifelse(MTposition[n]==9|MTposition[n]==8,charge(substr(MT.pep[n],7,9)),charge(substr(MT.pep[n],MTposition[n],MTposition[n]+2))),
      ifelse(MTposition[n]==9|MTposition[n]==8,mw(substr(MT.pep[n],7,9)),mw(substr(MT.pep[n],MTposition[n],MTposition[n]+2))),
      ifelse(MTposition[n]==9|MTposition[n]==8,boman(substr(MT.pep[n],7,9)),boman(substr(MT.pep[n],MTposition[n],MTposition[n]+2))),
      ifelse(MTposition[n]==9|MTposition[n]==8,aIndex(substr(MT.pep[n],7,9)),aIndex(substr(MT.pep[n],MTposition[n],MTposition[n]+2))),
      ifelse(MTposition[n]==9|MTposition[n]==8,hmoment(substr(MT.pep[n],7,9)),hmoment(substr(MT.pep[n],MTposition[n],MTposition[n]+2))),
      
      hydrophobicity(substr(MT.pep[n],1,1)),
      hydrophobicity(substr(MT.pep[n],2,2)),
      hydrophobicity(substr(MT.pep[n],3,3)),
      hydrophobicity(substr(MT.pep[n],4,4)),
      hydrophobicity(substr(MT.pep[n],5,5)),
      hydrophobicity(substr(MT.pep[n],6,6)),
      hydrophobicity(substr(MT.pep[n],7,7)),
      hydrophobicity(substr(MT.pep[n],8,8)),
      hydrophobicity(substr(MT.pep[n],9,9)),
      hydrophobicity(substr(MT.pep[n],1,2)),
      hydrophobicity(substr(MT.pep[n],2,3)),
      hydrophobicity(substr(MT.pep[n],3,4)),
      hydrophobicity(substr(MT.pep[n],4,5)),
      hydrophobicity(substr(MT.pep[n],5,6)),
      hydrophobicity(substr(MT.pep[n],6,7)),
      hydrophobicity(substr(MT.pep[n],7,8)),
      hydrophobicity(substr(MT.pep[n],8,9)),
      hydrophobicity(substr(MT.pep[n],1,3)),
      hydrophobicity(substr(MT.pep[n],2,4)),
      hydrophobicity(substr(MT.pep[n],3,5)),
      hydrophobicity(substr(MT.pep[n],4,6)),
      hydrophobicity(substr(MT.pep[n],5,7)),
      hydrophobicity(substr(MT.pep[n],6,8)),
      hydrophobicity(substr(MT.pep[n],7,9)),
      charge(substr(MT.pep[n],1,1)),
      charge(substr(MT.pep[n],2,2)),
      charge(substr(MT.pep[n],3,3)),
      charge(substr(MT.pep[n],4,4)),
      charge(substr(MT.pep[n],5,5)),
      charge(substr(MT.pep[n],6,6)),
      charge(substr(MT.pep[n],7,7)),
      charge(substr(MT.pep[n],8,8)),
      charge(substr(MT.pep[n],9,9)),
      charge(substr(MT.pep[n],1,2)),
      charge(substr(MT.pep[n],2,3)),
      charge(substr(MT.pep[n],3,4)),
      charge(substr(MT.pep[n],4,5)),
      charge(substr(MT.pep[n],5,6)),
      charge(substr(MT.pep[n],6,7)),
      charge(substr(MT.pep[n],7,8)),
      charge(substr(MT.pep[n],8,9)),
      charge(substr(MT.pep[n],1,3)),
      charge(substr(MT.pep[n],2,4)),
      charge(substr(MT.pep[n],3,5)),
      charge(substr(MT.pep[n],4,6)),
      charge(substr(MT.pep[n],5,7)),
      charge(substr(MT.pep[n],6,8)),
      charge(substr(MT.pep[n],7,9)),
      mw(substr(MT.pep[n],1,1)),
      mw(substr(MT.pep[n],2,2)),
      mw(substr(MT.pep[n],3,3)),
      mw(substr(MT.pep[n],4,4)),
      mw(substr(MT.pep[n],5,5)),
      mw(substr(MT.pep[n],6,6)),
      mw(substr(MT.pep[n],7,7)),
      mw(substr(MT.pep[n],8,8)),
      mw(substr(MT.pep[n],9,9)),
      mw(substr(MT.pep[n],1,2)),
      mw(substr(MT.pep[n],2,3)),
      mw(substr(MT.pep[n],3,4)),
      mw(substr(MT.pep[n],4,5)),
      mw(substr(MT.pep[n],5,6)),
      mw(substr(MT.pep[n],6,7)),
      mw(substr(MT.pep[n],7,8)),
      mw(substr(MT.pep[n],8,9)),
      mw(substr(MT.pep[n],1,3)),
      mw(substr(MT.pep[n],2,4)),
      mw(substr(MT.pep[n],3,5)),
      mw(substr(MT.pep[n],4,6)),
      mw(substr(MT.pep[n],5,7)),
      mw(substr(MT.pep[n],6,8)),
      mw(substr(MT.pep[n],7,9)),
      boman(substr(MT.pep[n],1,1)),
      boman(substr(MT.pep[n],2,2)),
      boman(substr(MT.pep[n],3,3)),
      boman(substr(MT.pep[n],4,4)),
      boman(substr(MT.pep[n],5,5)),
      boman(substr(MT.pep[n],6,6)),
      boman(substr(MT.pep[n],7,7)),
      boman(substr(MT.pep[n],8,8)),
      boman(substr(MT.pep[n],9,9)),
      boman(substr(MT.pep[n],1,2)),
      boman(substr(MT.pep[n],2,3)),
      boman(substr(MT.pep[n],3,4)),
      boman(substr(MT.pep[n],4,5)),
      boman(substr(MT.pep[n],5,6)),
      boman(substr(MT.pep[n],6,7)),
      boman(substr(MT.pep[n],7,8)),
      boman(substr(MT.pep[n],8,9)),
      boman(substr(MT.pep[n],1,3)),
      boman(substr(MT.pep[n],2,4)),
      boman(substr(MT.pep[n],3,5)),
      boman(substr(MT.pep[n],4,6)),
      boman(substr(MT.pep[n],5,7)),
      boman(substr(MT.pep[n],6,8)),
      boman(substr(MT.pep[n],7,9)),
      aIndex(substr(MT.pep[n],1,1)),
      aIndex(substr(MT.pep[n],2,2)),
      aIndex(substr(MT.pep[n],3,3)),
      aIndex(substr(MT.pep[n],4,4)),
      aIndex(substr(MT.pep[n],5,5)),
      aIndex(substr(MT.pep[n],6,6)),
      aIndex(substr(MT.pep[n],7,7)),
      aIndex(substr(MT.pep[n],8,8)),
      aIndex(substr(MT.pep[n],9,9)),
      aIndex(substr(MT.pep[n],1,2)),
      aIndex(substr(MT.pep[n],2,3)),
      aIndex(substr(MT.pep[n],3,4)),
      aIndex(substr(MT.pep[n],4,5)),
      aIndex(substr(MT.pep[n],5,6)),
      aIndex(substr(MT.pep[n],6,7)),
      aIndex(substr(MT.pep[n],7,8)),
      aIndex(substr(MT.pep[n],8,9)),
      aIndex(substr(MT.pep[n],1,3)),
      aIndex(substr(MT.pep[n],2,4)),
      aIndex(substr(MT.pep[n],3,5)),
      aIndex(substr(MT.pep[n],4,6)),
      aIndex(substr(MT.pep[n],5,7)),
      aIndex(substr(MT.pep[n],6,8)),
      aIndex(substr(MT.pep[n],7,9)),
      hmoment(substr(MT.pep[n],1,1)),
      hmoment(substr(MT.pep[n],2,2)),
      hmoment(substr(MT.pep[n],3,3)),
      hmoment(substr(MT.pep[n],4,4)),
      hmoment(substr(MT.pep[n],5,5)),
      hmoment(substr(MT.pep[n],6,6)),
      hmoment(substr(MT.pep[n],7,7)),
      hmoment(substr(MT.pep[n],8,8)),
      hmoment(substr(MT.pep[n],9,9)),
      hmoment(substr(MT.pep[n],1,2)),
      hmoment(substr(MT.pep[n],2,3)),
      hmoment(substr(MT.pep[n],3,4)),
      hmoment(substr(MT.pep[n],4,5)),
      hmoment(substr(MT.pep[n],5,6)),
      hmoment(substr(MT.pep[n],6,7)),
      hmoment(substr(MT.pep[n],7,8)),
      hmoment(substr(MT.pep[n],8,9)),
      hmoment(substr(MT.pep[n],1,3)),
      hmoment(substr(MT.pep[n],2,4)),
      hmoment(substr(MT.pep[n],3,5)),
      hmoment(substr(MT.pep[n],4,6)),
      hmoment(substr(MT.pep[n],5,7)),
      hmoment(substr(MT.pep[n],6,8)),
      hmoment(substr(MT.pep[n],7,9)),
      instaIndex(substr(MT.pep[n],1,3)),
      instaIndex(substr(MT.pep[n],2,4)),
      instaIndex(substr(MT.pep[n],3,5)),
      instaIndex(substr(MT.pep[n],4,6)),
      instaIndex(substr(MT.pep[n],5,7)),
      instaIndex(substr(MT.pep[n],6,8)),
      instaIndex(substr(MT.pep[n],7,9)),
      
      aaComp(MT.pep[n])[[1]][1:9],
      aaComp(substring(MT.pep[n],1,3))[[1]][1:9],
      aaComp(substring(MT.pep[n],2,4))[[1]][1:9],
      aaComp(substring(MT.pep[n],3,5))[[1]][1:9],
      aaComp(substring(MT.pep[n],4,6))[[1]][1:9],
      aaComp(substring(MT.pep[n],5,7))[[1]][1:9],
      aaComp(substring(MT.pep[n],6,8))[[1]][1:9],
      aaComp(substring(MT.pep[n],7,9))[[1]][1:9],
      ifelse(MTposition[n]==1|MTposition[n]==2,aaComp(substr(MT.pep[n],1,3)),aaComp(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:9],
      ifelse(MTposition[n]==9|MTposition[n]==8,aaComp(substr(MT.pep[n],7,9)),aaComp(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:9],
      
      crucianiProperties(MT.pep[n])[[1]][1:3],
      kideraFactors(MT.pep[n])[[1]][1:10],
      zScales(MT.pep[n])[[1]][1:5],
      fasgaiVectors(MT.pep[n])[[1]][1:6],
      tScales(MT.pep[n])[[1]][1:5],
      vhseScales(MT.pep[n])[[1]][1:8],
      protFP(MT.pep[n])[[1]][1:8],
      stScales(MT.pep[n])[[1]][1:8],
      blosumIndices(MT.pep[n])[[1]][1:10],
      mswhimScores(MT.pep[n])[[1]][1:3],
      
      aaDescriptors(substr(MT.pep[n],MTposition[n],MTposition[n])),
      
      ifelse(MTposition[n]==1,crucianiProperties(substr(MT.pep[n],1,2)),crucianiProperties(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:3],
      ifelse(MTposition[n]==1, kideraFactors(substr(MT.pep[n],1,2)), kideraFactors(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:10],
      ifelse(MTposition[n]==1,zScales(substr(MT.pep[n],1,2)),zScales(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:5],
      ifelse(MTposition[n]==1,fasgaiVectors(substr(MT.pep[n],1,2)),fasgaiVectors(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:6],
      ifelse(MTposition[n]==1,tScales(substr(MT.pep[n],1,2)),tScales(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:5],
      ifelse(MTposition[n]==1,vhseScales(substr(MT.pep[n],1,2)),vhseScales(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:8],
      ifelse(MTposition[n]==1, protFP(substr(MT.pep[n],1,2)), protFP(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:8],
      ifelse(MTposition[n]==1,stScales(substr(MT.pep[n],1,2)),stScales(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:8],
      ifelse(MTposition[n]==1,blosumIndices(substr(MT.pep[n],1,2)),blosumIndices(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:10],
      ifelse(MTposition[n]==1, mswhimScores(substr(MT.pep[n],1,2)), mswhimScores(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:3],
      #i;i+1
      ifelse(MTposition[n]==9,crucianiProperties(substr(MT.pep[n],8,9)),crucianiProperties(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:3],
      ifelse(MTposition[n]==9, kideraFactors(substr(MT.pep[n],8,9)), kideraFactors(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:10],
      ifelse(MTposition[n]==9,zScales(substr(MT.pep[n],8,9)),zScales(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:5],
      ifelse(MTposition[n]==9,fasgaiVectors(substr(MT.pep[n],8,9)),fasgaiVectors(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:6],
      ifelse(MTposition[n]==9,tScales(substr(MT.pep[n],8,9)),tScales(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:5],
      ifelse(MTposition[n]==9,vhseScales(substr(MT.pep[n],8,9)),vhseScales(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:8],
      ifelse(MTposition[n]==9, protFP(substr(MT.pep[n],8,9)), protFP(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:8],
      ifelse(MTposition[n]==9,stScales(substr(MT.pep[n],8,9)),stScales(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:8],
      ifelse(MTposition[n]==9,blosumIndices(substr(MT.pep[n],8,9)),blosumIndices(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:10],
      ifelse(MTposition[n]==9, mswhimScores(substr(MT.pep[n],8,9)), mswhimScores(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:3],
      
      ifelse(MTposition[n]==1|MTposition[n]==2,crucianiProperties(substr(MT.pep[n],1,3)),crucianiProperties(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:3],
      ifelse(MTposition[n]==1|MTposition[n]==2, kideraFactors(substr(MT.pep[n],1,3)), kideraFactors(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:10],
      ifelse(MTposition[n]==1|MTposition[n]==2,zScales(substr(MT.pep[n],1,3)),zScales(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:5],
      ifelse(MTposition[n]==1|MTposition[n]==2,fasgaiVectors(substr(MT.pep[n],1,3)),fasgaiVectors(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:6],
      ifelse(MTposition[n]==1|MTposition[n]==2,tScales(substr(MT.pep[n],1,3)),tScales(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:5],
      ifelse(MTposition[n]==1|MTposition[n]==2,vhseScales(substr(MT.pep[n],1,3)),vhseScales(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:8],
      ifelse(MTposition[n]==1|MTposition[n]==2, protFP(substr(MT.pep[n],1,3)), protFP(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:8],
      ifelse(MTposition[n]==1|MTposition[n]==2,stScales(substr(MT.pep[n],1,3)),stScales(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:8],
      ifelse(MTposition[n]==1|MTposition[n]==2,blosumIndices(substr(MT.pep[n],1,3)),blosumIndices(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:10],
      ifelse(MTposition[n]==1|MTposition[n]==2, mswhimScores(substr(MT.pep[n],1,3)),mswhimScores(substr(MT.pep[n],MTposition[n]-2,MTposition[n])))[[1]][1:3],
      
      ifelse(MTposition[n]==9|MTposition[n]==8,crucianiProperties(substr(MT.pep[n],7,9)),crucianiProperties(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:3],
      ifelse(MTposition[n]==9|MTposition[n]==8, kideraFactors(substr(MT.pep[n],7,9)), kideraFactors(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:10],
      ifelse(MTposition[n]==9|MTposition[n]==8,zScales(substr(MT.pep[n],7,9)),zScales(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:5],
      ifelse(MTposition[n]==9|MTposition[n]==8,fasgaiVectors(substr(MT.pep[n],7,9)),fasgaiVectors(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:6],
      ifelse(MTposition[n]==9|MTposition[n]==8,tScales(substr(MT.pep[n],7,9)),tScales(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:5],
      ifelse(MTposition[n]==9|MTposition[n]==8,vhseScales(substr(MT.pep[n],7,9)),vhseScales(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:8],
      ifelse(MTposition[n]==9|MTposition[n]==8, protFP(substr(MT.pep[n],7,9)), protFP(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:8],
      ifelse(MTposition[n]==9|MTposition[n]==8,stScales(substr(MT.pep[n],7,9)),stScales(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:8],
      ifelse(MTposition[n]==9|MTposition[n]==8,blosumIndices(substr(MT.pep[n],7,9)),blosumIndices(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:10],
      ifelse(MTposition[n]==9|MTposition[n]==8, mswhimScores(substr(MT.pep[n],7,9)), mswhimScores(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))[[1]][1:3],
      
      ##############################################
      
      hydrophobicity(WT.pep[n])-hydrophobicity(MT.pep[n]),
      instaIndex(WT.pep[n])-instaIndex(MT.pep[n]),
      charge(WT.pep[n])-charge(MT.pep[n]),
      mw(WT.pep[n])- mw(MT.pep[n]),
      boman(WT.pep[n])- boman(MT.pep[n]),
      aIndex(WT.pep[n])-aIndex(MT.pep[n]),
      autoCorrelation(WT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,center = TRUE)-autoCorrelation(MT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,center = TRUE),
      autoCovariance(WT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,center = TRUE)-autoCovariance(MT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,center = TRUE),
      hmoment(WT.pep[n])- hmoment(MT.pep[n]),
      crossCovariance(WT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,property2 = AAdata$Hydrophobicity$Eisenberg,center = TRUE)-crossCovariance(MT.pep[n],lag = 1,property = AAdata$Hydrophobicity$KyteDoolittle,property2 = AAdata$Hydrophobicity$Eisenberg,center = TRUE),
      hydrophobicity(substring(WT.pep[n],MTposition[n],MTposition[n]))-hydrophobicity(substring(MT.pep[n],MTposition[n],MTposition[n])),
      charge(substring(WT.pep[n],MTposition[n],MTposition[n]))-charge(substring(MT.pep[n],MTposition[n],MTposition[n])),
      mw(substring(WT.pep[n],MTposition[n],MTposition[n]))-mw(substring(MT.pep[n],MTposition[n],MTposition[n])),
      boman(substring(WT.pep[n],MTposition[n],MTposition[n]))-boman(substring(MT.pep[n],MTposition[n],MTposition[n])),
      aIndex(substring(WT.pep[n],MTposition[n],MTposition[n]))-aIndex(substring(MT.pep[n],MTposition[n],MTposition[n])),
      hmoment(substring(WT.pep[n],MTposition[n],MTposition[n]))-hmoment(substring(MT.pep[n],MTposition[n],MTposition[n])),
      
      aaComp(WT.pep[n])[[1]][1:9]-aaComp(MT.pep[n])[[1]][1:9],
      aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][1:9],
      aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][1:9]-aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][1:9],
      ifelse(MTposition[n]==1,aaComp(substr(MT.pep[n],1,2)),aaComp(substr(MT.pep[n],MTposition[n]-1,MTposition[n])))[[1]][1:9],
      ifelse(MTposition[n]==9,aaComp(substr(MT.pep[n],8,9)),aaComp(substr(MT.pep[n],MTposition[n],MTposition[n]+1)))[[1]][1:9],
      aaComp(substring(MT.pep[n],1,2))[[1]][1:9],
      aaComp(substring(MT.pep[n],2,3))[[1]][1:9],
      aaComp(substring(MT.pep[n],3,4))[[1]][1:9],
      aaComp(substring(MT.pep[n],4,5))[[1]][1:9],
      aaComp(substring(MT.pep[n],5,6))[[1]][1:9],
      aaComp(substring(MT.pep[n],6,7))[[1]][1:9],
      aaComp(substring(MT.pep[n],7,8))[[1]][1:9],
      aaComp(substring(MT.pep[n],8,9))[[1]][1:9],
      aaComp(substring(MT.pep[n],1,1))[[1]][1:9],
      aaComp(substring(MT.pep[n],2,2))[[1]][1:9],
      aaComp(substring(MT.pep[n],3,3))[[1]][1:9],
      aaComp(substring(MT.pep[n],4,4))[[1]][1:9],
      aaComp(substring(MT.pep[n],5,5))[[1]][1:9],
      aaComp(substring(MT.pep[n],6,6))[[1]][1:9],
      aaComp(substring(MT.pep[n],7,7))[[1]][1:9],
      aaComp(substring(MT.pep[n],8,8))[[1]][1:9],
      aaComp(substring(MT.pep[n],9,9))[[1]][1:9],
      ifelse(aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][1] & aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][2],1,0),
      ifelse(aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][2] & aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][1],1,0),
      ifelse(aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][3] & aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][4],1,0),
      ifelse(aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][4] & aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][3],1,0),
      ifelse(aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][5] & aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][6],1,0),
      ifelse(aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][6] & aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][5],1,0),
      ifelse(aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][8] & aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][9],1,0),
      ifelse(aaComp(substring(WT.pep[n],MTposition[n],MTposition[n]))[[1]][9] & aaComp(substring(MT.pep[n],MTposition[n],MTposition[n]))[[1]][8],1,0),
      
      ifelse(substr(MT.pep[n],MTposition[n],MTposition[n]) == Peptype,1,0),    #20??
      ifelse(Peptype %in% unlist(strsplit(substr(MT.pep[n],1,3),"|")),1,0),  #\u{37b}????\u{1f0}????��?\u{f4}???XX????
      ifelse(Peptype %in% unlist(strsplit(substr(MT.pep[n],4,6),"|")),1,0),  #\u{37b}?????��?????��?\u{f4}???XX????
      ifelse(Peptype %in% unlist(strsplit(substr(MT.pep[n],7,9),"|")),1,0),  #\u{37b}????????????��?\u{f4}???XX????
      ifelse(Peptype %in% unlist(strsplit(substr(MT.pep[n],1,2),"|")),1,0),  #\u{37b}????\u{1f0}2??��?\u{f4}???XX????
      ifelse(Peptype %in% unlist(strsplit(substr(MT.pep[n],8,9),"|")),1,0),  #\u{37b}???\u{13a}?2??��?\u{f4}???XX????
      ifelse(substr(MT.pep[n],1,1) == Peptype,1,0),    #20??
      ifelse(substr(MT.pep[n],9,9) == Peptype,1,0),    #20
      
      crucianiProperties(WT.pep[n])[[1]][1:3]-crucianiProperties(MT.pep[n])[[1]][1:3],
      kideraFactors(WT.pep[n])[[1]][1:10]-kideraFactors(MT.pep[n])[[1]][1:10],
      zScales(WT.pep[n])[[1]][1:5]-zScales(MT.pep[n])[[1]][1:5],
      fasgaiVectors(WT.pep[n])[[1]][1:6]-fasgaiVectors(MT.pep[n])[[1]][1:6],
      tScales(WT.pep[n])[[1]][1:5]-tScales(MT.pep[n])[[1]][1:5],
      vhseScales(WT.pep[n])[[1]][1:8]-vhseScales(MT.pep[n])[[1]][1:8],
      protFP(WT.pep[n])[[1]][1:8]-protFP(MT.pep[n])[[1]][1:8],
      stScales(WT.pep[n])[[1]][1:8]-stScales(MT.pep[n])[[1]][1:8],
      blosumIndices(WT.pep[n])[[1]][1:10]-blosumIndices(MT.pep[n])[[1]][1:10],
      mswhimScores(WT.pep[n])[[1]][1:3]-mswhimScores(MT.pep[n])[[1]][1:3],   #  above  all  66
      aaDescriptors(substr(WT.pep[n],MTposition[n],MTposition[n]))- aaDescriptors(substr(MT.pep[n],MTposition[n],MTposition[n])),   #66
      
      ifelse(MTposition[n]<=3,1,0),
      ifelse(MTposition[n]>=4 & MTposition[n]<=6,1,0),
      ifelse(MTposition[n]>=7 & MTposition[n]<=9,1,0)
    )
  }
  model_mat_all1 = foreach(n = 1:length(MT.pep),.combine = rbind) %dopar% model_process(n)
  model_mat_all1<-data.frame(model_mat_all1)
  #### (2)
  There_PF1<-NULL  #i-1;i;i+1
  for (n in 1:length(MT.pep)){
    Ai<-NULL
    if(MTposition[n]==1){
      Ai<-c(
        hydrophobicity(substr(MT.pep[n],1,3)),
        charge(substr(MT.pep[n],1,3)),
        mw(substr(MT.pep[n],1,3)),
        boman(substr(MT.pep[n],1,3)),
        aIndex(substr(MT.pep[n],1,3)),
        hmoment(substr(MT.pep[n],1,3)),
        aaComp(substring(MT.pep[n],1,3))[[1]][1:9],
        crucianiProperties(substr(MT.pep[n],1,3))[[1]][1:3],
        kideraFactors(substr(MT.pep[n],1,3))[[1]][1:10],
        zScales(substr(MT.pep[n],1,3))[[1]][1:5],
        fasgaiVectors(substr(MT.pep[n],1,3))[[1]][1:6],
        tScales(substr(MT.pep[n],1,3))[[1]][1:5],
        vhseScales(substr(MT.pep[n],1,3))[[1]][1:8],
        protFP(substr(MT.pep[n],1,3))[[1]][1:8],
        stScales(substr(MT.pep[n],1,3))[[1]][1:8],
        blosumIndices(substr(MT.pep[n],1,3))[[1]][1:10],
        mswhimScores(substr(MT.pep[n],1,3))[[1]][1:3]
      )
    }else if(MTposition[n]==9){
      Ai<-c(
        hydrophobicity(substr(MT.pep[n],7,9)),
        charge(substr(MT.pep[n],7,9)),
        mw(substr(MT.pep[n],7,9)),
        boman(substr(MT.pep[n],7,9)),
        aIndex(substr(MT.pep[n],7,9)),
        hmoment(substr(MT.pep[n],7,9)),
        aaComp(substring(MT.pep[n],7,9))[[1]][1:9],
        crucianiProperties(substr(MT.pep[n],7,9))[[1]][1:3],
        kideraFactors(substr(MT.pep[n],7,9))[[1]][1:10],
        zScales(substr(MT.pep[n],7,9))[[1]][1:5],
        fasgaiVectors(substr(MT.pep[n],7,9))[[1]][1:6],
        tScales(substr(MT.pep[n],7,9))[[1]][1:5],
        vhseScales(substr(MT.pep[n],7,9))[[1]][1:8],
        protFP(substr(MT.pep[n],7,9))[[1]][1:8],
        stScales(substr(MT.pep[n],7,9))[[1]][1:8],
        blosumIndices(substr(MT.pep[n],7,9))[[1]][1:10],
        mswhimScores(substr(MT.pep[n],7,9))[[1]][1:3]
      )
    } else{
      Ai<-c(
        hydrophobicity(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1)),
        charge(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1)),
        mw(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1)),
        boman(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1)),
        aIndex(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1)),
        hmoment(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1)),
        aaComp(substring(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:9],
        
        crucianiProperties(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:3],
        kideraFactors(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:10],
        zScales(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:5],
        fasgaiVectors(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:6],
        tScales(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:5],
        vhseScales(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:8],
        protFP(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:8],
        stScales(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:8],
        blosumIndices(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:10],
        mswhimScores(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))[[1]][1:3]
      )
    }
    There_PF1<-rbind(There_PF1,Ai)  #mean
  }
  There_PF1<-data.frame(There_PF1)
  #### (3)
  MT.Entropy<-entropyseq(MT.pep)
  WT.Entropy<-entropyseq(WT.pep)
  WMT.Entropy<-WT.Entropy-MT.Entropy
  WTR.Entropy<-entropyresidue(WT.pep)
  entropy_0<-cbind(WMT.Entropy,WTR.Entropy) #21
  model_process = function(n){  #4
    c(
      #mean
      ifelse(MTposition[n]==1,entropyseq(substr(MT.pep[n],1,2)),entropyseq(substr(MT.pep[n],MTposition[n]-1,MTposition[n]))),
      ifelse(MTposition[n]==9,entropyseq(substr(MT.pep[n],8,9)),entropyseq(substr(MT.pep[n],MTposition[n],MTposition[n]+1))),
      ifelse(MTposition[n]==1|MTposition[n]==2,entropyseq(substr(MT.pep[n],1,3)),entropyseq(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==9|MTposition[n]==8,entropyseq(substr(MT.pep[n],7,9)),entropyseq(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))
      
    )
  }
  model_entropy_mean<- foreach(n = 1:length(MT.pep),.combine = rbind) %dopar% model_process(n)
  model_entropy_mean<-data.frame(model_entropy_mean)
  
  There_entropy<-NULL
  for (n in 1:length(MT.pep)){
    Bi<-NULL
    if(MTposition[n]==1){
      Bi<-c(
        entropyseq(substr(MT.pep[n],1,3)),
        entropyseq(substr(WT.pep[n],1,3))-entropyseq(substr(MT.pep[n],1,3))  #  0
        
      )
    }else if(MTposition[n]==9){
      Bi<-c(
        entropyseq(substr(MT.pep[n],7,9)),
        entropyseq(substr(WT.pep[n],7,9))-entropyseq(substr(MT.pep[n],7,9))  #  0
        
      )
    } else{
      Bi<-c(
        entropyseq(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1)),
        entropyseq(substr(WT.pep[n],MTposition[n]-1,MTposition[n]+1))-entropyseq(substr(MT.pep[n],MTposition[n]-1,MTposition[n]+1))  #  0
        
      )
    }
    There_entropy<-rbind(There_entropy,Bi)
  }
  model_process = function(n){
    c(
      ifelse(MTposition[n]==1,entropyseq(substr(WT.pep[n],1,2))-entropyseq(substr(MT.pep[n],1,2)),entropyseq(substr(WT.pep[n],MTposition[n]-1,MTposition[n]))-entropyseq(substr(MT.pep[n],MTposition[n]-1,MTposition[n]))),
      ifelse(MTposition[n]==9,entropyseq(substr(WT.pep[n],8,9))-entropyseq(substr(MT.pep[n],8,9)),entropyseq(substr(WT.pep[n],MTposition[n],MTposition[n]+1))-entropyseq(substr(MT.pep[n],MTposition[n],MTposition[n]+1))),
      ifelse(MTposition[n]==1|MTposition[n]==2,entropyseq(substr(WT.pep[n],1,3))-entropyseq(substr(MT.pep[n],1,3)),entropyseq(substr(WT.pep[n],MTposition[n]-2,MTposition[n]))-entropyseq(substr(MT.pep[n],MTposition[n]-2,MTposition[n]))),
      ifelse(MTposition[n]==9|MTposition[n]==8,entropyseq(substr(WT.pep[n],7,9))-entropyseq(substr(MT.pep[n],7,9)),entropyseq(substr(WT.pep[n],MTposition[n],MTposition[n]+2))-entropyseq(substr(MT.pep[n],MTposition[n],MTposition[n]+2)))
      
    )
  }
  model_entropy_zero<- foreach(n = 1:length(MT.pep),.combine = rbind) %dopar% model_process(n)
  model_entropy_zero<-data.frame(model_entropy_zero)
  Entropy_all1<-cbind(MT.Entropy,model_entropy_mean,There_entropy,model_entropy_zero,entropy_0)
  
  #### (4)
  MTR.Entropy<-entropyresidue(MT.pep)
  Entropy2<-MTR.Entropy
  #### (5)
  model_process = function(n){
    c(
      aaDescriptors(substr(MT.pep[n],1,1)),
      aaDescriptors(substr(MT.pep[n],2,2)),
      aaDescriptors(substr(MT.pep[n],3,3)),
      aaDescriptors(substr(MT.pep[n],4,4)),
      aaDescriptors(substr(MT.pep[n],5,5)),
      aaDescriptors(substr(MT.pep[n],6,6)),
      aaDescriptors(substr(MT.pep[n],7,7)),
      aaDescriptors(substr(MT.pep[n],8,8)),
      aaDescriptors(substr(MT.pep[n],9,9)),
      #\u{1f0}\u{fffd}\u{fffd}��
      crucianiProperties(substr(MT.pep[n],1,2))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],1,2))[[1]][1:10],
      zScales(substr(MT.pep[n],1,2))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],1,2))[[1]][1:6],
      tScales(substr(MT.pep[n],1,2))[[1]][1:5],
      vhseScales(substr(MT.pep[n],1,2))[[1]][1:8],
      protFP(substr(MT.pep[n],1,2))[[1]][1:8],
      stScales(substr(MT.pep[n],1,2))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],1,2))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],1,2))[[1]][1:3],         #  above  all  66
      #??\u{fffd}\u{fffd}��
      crucianiProperties(substr(MT.pep[n],8,9))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],8,9))[[1]][1:10],
      zScales(substr(MT.pep[n],8,9))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],8,9))[[1]][1:6],
      tScales(substr(MT.pep[n],8,9))[[1]][1:5],
      vhseScales(substr(MT.pep[n],8,9))[[1]][1:8],
      protFP(substr(MT.pep[n],8,9))[[1]][1:8],
      stScales(substr(MT.pep[n],8,9))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],8,9))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],8,9))[[1]][1:3],         #  above  all  66
      #\u{1f0}??��
      crucianiProperties(substr(MT.pep[n],1,3))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],1,3))[[1]][1:10],
      zScales(substr(MT.pep[n],1,3))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],1,3))[[1]][1:6],
      tScales(substr(MT.pep[n],1,3))[[1]][1:5],
      vhseScales(substr(MT.pep[n],1,3))[[1]][1:8],
      protFP(substr(MT.pep[n],1,3))[[1]][1:8],
      stScales(substr(MT.pep[n],1,3))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],1,3))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],1,3))[[1]][1:3],
      
      #?��???��
      crucianiProperties(substr(MT.pep[n],4,6))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],4,6))[[1]][1:10],
      zScales(substr(MT.pep[n],4,6))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],4,6))[[1]][1:6],
      tScales(substr(MT.pep[n],4,6))[[1]][1:5],
      vhseScales(substr(MT.pep[n],4,6))[[1]][1:8],
      protFP(substr(MT.pep[n],4,6))[[1]][1:8],
      stScales(substr(MT.pep[n],4,6))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],4,6))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],4,6))[[1]][1:3],
      
      #????��
      crucianiProperties(substr(MT.pep[n],7,9))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],7,9))[[1]][1:10],
      zScales(substr(MT.pep[n],7,9))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],7,9))[[1]][1:6],
      tScales(substr(MT.pep[n],7,9))[[1]][1:5],
      vhseScales(substr(MT.pep[n],7,9))[[1]][1:8],
      protFP(substr(MT.pep[n],7,9))[[1]][1:8],
      stScales(substr(MT.pep[n],7,9))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],7,9))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],7,9))[[1]][1:3],
      
      #2
      crucianiProperties(substr(MT.pep[n],2,3))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],2,3))[[1]][1:10],
      zScales(substr(MT.pep[n],2,3))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],2,3))[[1]][1:6],
      tScales(substr(MT.pep[n],2,3))[[1]][1:5],
      vhseScales(substr(MT.pep[n],2,3))[[1]][1:8],
      protFP(substr(MT.pep[n],2,3))[[1]][1:8],
      stScales(substr(MT.pep[n],2,3))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],2,3))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],2,3))[[1]][1:3],
      
      crucianiProperties(substr(MT.pep[n],3,4))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],3,4))[[1]][1:10],
      zScales(substr(MT.pep[n],3,4))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],3,4))[[1]][1:6],
      tScales(substr(MT.pep[n],3,4))[[1]][1:5],
      vhseScales(substr(MT.pep[n],3,4))[[1]][1:8],
      protFP(substr(MT.pep[n],3,4))[[1]][1:8],
      stScales(substr(MT.pep[n],3,4))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],3,4))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],3,4))[[1]][1:3],
      
      crucianiProperties(substr(MT.pep[n],4,5))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],4,5))[[1]][1:10],
      zScales(substr(MT.pep[n],4,5))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],4,5))[[1]][1:6],
      tScales(substr(MT.pep[n],4,5))[[1]][1:5],
      vhseScales(substr(MT.pep[n],4,5))[[1]][1:8],
      protFP(substr(MT.pep[n],4,5))[[1]][1:8],
      stScales(substr(MT.pep[n],4,5))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],4,5))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],4,5))[[1]][1:3],
      
      crucianiProperties(substr(MT.pep[n],5,6))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],5,6))[[1]][1:10],
      zScales(substr(MT.pep[n],5,6))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],5,6))[[1]][1:6],
      tScales(substr(MT.pep[n],5,6))[[1]][1:5],
      vhseScales(substr(MT.pep[n],5,6))[[1]][1:8],
      protFP(substr(MT.pep[n],5,6))[[1]][1:8],
      stScales(substr(MT.pep[n],5,6))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],5,6))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],5,6))[[1]][1:3],
      
      crucianiProperties(substr(MT.pep[n],6,7))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],6,7))[[1]][1:10],
      zScales(substr(MT.pep[n],6,7))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],6,7))[[1]][1:6],
      tScales(substr(MT.pep[n],6,7))[[1]][1:5],
      vhseScales(substr(MT.pep[n],6,7))[[1]][1:8],
      protFP(substr(MT.pep[n],6,7))[[1]][1:8],
      stScales(substr(MT.pep[n],6,7))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],6,7))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],6,7))[[1]][1:3],
      
      crucianiProperties(substr(MT.pep[n],7,8))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],7,8))[[1]][1:10],
      zScales(substr(MT.pep[n],7,8))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],7,8))[[1]][1:6],
      tScales(substr(MT.pep[n],7,8))[[1]][1:5],
      vhseScales(substr(MT.pep[n],7,8))[[1]][1:8],
      protFP(substr(MT.pep[n],7,8))[[1]][1:8],
      stScales(substr(MT.pep[n],7,8))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],7,8))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],7,8))[[1]][1:3],
      
      #3
      crucianiProperties(substr(MT.pep[n],2,4))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],2,4))[[1]][1:10],
      zScales(substr(MT.pep[n],2,4))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],2,4))[[1]][1:6],
      tScales(substr(MT.pep[n],2,4))[[1]][1:5],
      vhseScales(substr(MT.pep[n],2,4))[[1]][1:8],
      protFP(substr(MT.pep[n],2,4))[[1]][1:8],
      stScales(substr(MT.pep[n],2,4))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],2,4))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],2,4))[[1]][1:3],
      
      crucianiProperties(substr(MT.pep[n],3,5))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],3,5))[[1]][1:10],
      zScales(substr(MT.pep[n],3,5))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],3,5))[[1]][1:6],
      tScales(substr(MT.pep[n],3,5))[[1]][1:5],
      vhseScales(substr(MT.pep[n],3,5))[[1]][1:8],
      protFP(substr(MT.pep[n],3,5))[[1]][1:8],
      stScales(substr(MT.pep[n],3,5))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],3,5))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],3,5))[[1]][1:3],
      
      
      crucianiProperties(substr(MT.pep[n],5,7))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],5,7))[[1]][1:10],
      zScales(substr(MT.pep[n],5,7))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],5,7))[[1]][1:6],
      tScales(substr(MT.pep[n],5,7))[[1]][1:5],
      vhseScales(substr(MT.pep[n],5,7))[[1]][1:8],
      protFP(substr(MT.pep[n],5,7))[[1]][1:8],
      stScales(substr(MT.pep[n],5,7))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],5,7))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],5,7))[[1]][1:3],
      
      crucianiProperties(substr(MT.pep[n],6,8))[[1]][1:3],
      kideraFactors(substr(MT.pep[n],6,8))[[1]][1:10],
      zScales(substr(MT.pep[n],6,8))[[1]][1:5],
      fasgaiVectors(substr(MT.pep[n],6,8))[[1]][1:6],
      tScales(substr(MT.pep[n],6,8))[[1]][1:5],
      vhseScales(substr(MT.pep[n],6,8))[[1]][1:8],
      protFP(substr(MT.pep[n],6,8))[[1]][1:8],
      stScales(substr(MT.pep[n],6,8))[[1]][1:8],
      blosumIndices(substr(MT.pep[n],6,8))[[1]][1:10],
      mswhimScores(substr(MT.pep[n],6,8))[[1]][1:3]
    )
    
  }
  model_mat_all2 = foreach(n = 1:length(MT.pep),.combine = rbind) %dopar% model_process(n)
  model_mat_all2<-data.frame(model_mat_all2)
  
  #### (6)
  model_process = function(n){  #4
    c(
      #mean
      entropyseq(substr(MT.pep[n],1,2)),
      entropyseq(substr(MT.pep[n],2,3)),
      entropyseq(substr(MT.pep[n],3,4)),
      entropyseq(substr(MT.pep[n],4,5)),
      entropyseq(substr(MT.pep[n],5,6)),
      entropyseq(substr(MT.pep[n],6,7)),
      entropyseq(substr(MT.pep[n],7,8)),
      entropyseq(substr(MT.pep[n],8,9)),
      entropyseq(substr(MT.pep[n],1,3)),
      entropyseq(substr(MT.pep[n],2,4)),
      entropyseq(substr(MT.pep[n],3,5)),
      entropyseq(substr(MT.pep[n],4,6)),
      entropyseq(substr(MT.pep[n],5,7)),
      entropyseq(substr(MT.pep[n],6,8)),
      entropyseq(substr(MT.pep[n],7,9))
      
    )
  }
  Entropy_all2 = foreach(n = 1:length(MT.pep),.combine = rbind) %dopar% model_process(n)
  Entropy_all2<-data.frame(Entropy_all2)
  
  features<-data.frame(model_mat_all1,There_PF1,Entropy_all1,Entropy2,model_mat_all2,Entropy_all2)
  neoantigens<-data.frame(IDHwt,features)
  return (neoantigens)
  
}




library(dplyr)
library(readxl)


result <- PF_calculation(MT.pep,WT.pep)
result <- result[-nrow(result),]



colnames(result)[4:2931] <- a
fea <- result
des <- sapply(colnames(fea), function(x) unlist(strsplit(x,split=" "))[1])
library(stringr)

zScale <- grep("^zScales",colnames(fea))

zScale_result <-  cbind(result[,1:3],result[,zScale])

write.csv(zScale_result,file =  "zScale_result.csv")



