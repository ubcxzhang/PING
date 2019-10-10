#PING<-function(segReadsList,detail=0, rescale=1, PE=F)
PING<-function(segReadsList, paraEM=NULL, paraPrior=NULL, dataType="MNase", detail=0, rescale=1, nCores=1)
{
# detail is an integer indicate how much detail to print
# detail=0 print no details, The larger detail, the more detail I print 
# rescale is an integer indicate how reads was rescaled by yF/rescale and yR/rescale, default rescale=1
# PE=F or PE=0 means single-end sequencing data, otherwise Paired-Ends sequencing data.

  #get PE/SE
  if(class(segReadsList)=="segReadsListPE")
    PE<-TRUE
  else
    PE<-FALSE

  ### Constant used in the calculations
  cst<-gamma(3.5)/gamma(3)/sqrt(pi)
  minReads<-list(perPeak=3,perRegion=4)
  detail=as.integer(detail)
  rescale=as.integer(rescale)
  PE=as.integer(PE)

  ### Default parameters for EM algorithm
  if(length(paraEM)!=7)
  {
    paraEM<-setParaEM(dataType=dataType)
  }
  if(length(paraPrior)!=6)
  {
    if(PE>0) #==TRUE
    {
      paraPrior<-setParaPrior(dataType=dataType, PExi=segReadsList@paraSW$xi) #average reads length is already known
    }
    else
    {
      paraPrior<-setParaPrior(dataType=dataType)
    }
  }
  calpha <- 1.5
  
  if(nCores>1 & "parallel" %in% names(getLoadedDLLs()) )
  {
	#Number of cores
	availCores<-parallel:::detectCores()
	if(nCores > availCores)
        {
	  warning("The number of cores required is higher than the available cores on this machine (",availCores,").\n", immediate.=TRUE)
	  nCores<-availCores
        }
	message("Using the parallel version of PING with ", nCores, " cpus or cores")
	#Split into nCores segReadsList
	cl <- parallel:::makeCluster(getOption("cl.cores", nCores))
	segSplit<-split(segReadsList,cut(1:length(segReadsList),nCores))
	#Use parallel version of lapply
	res<-unlist(parallel:::parLapply(cl,segSplit,.fitModelAllkSplit,paraEM,paraPrior,minReads,detail,rescale, calpha, PE),recursive=FALSE) #::: in case parallel is imported by another package (in that case, parLapply would not be exposed)
	parallel:::stopCluster(cl)
  }
  else
  {
	message("Using the serial version of PING")
	res<-.Call("fitPING", segReadsList, paraEM, paraPrior, minReads, detail,rescale, calpha, PE, PACKAGE="PING")
  }

  myPingList<-newPingList(res,paraEM,paraPrior,minReads,segReadsList@N,segReadsList@Nc, PE)
  return(myPingList)
}

.fitModelAllkSplit<-function(segReadsList,paraEM,paraPrior,minReads,detail,rescale, calpha, PE)
{
  res<-.Call("fitPING", segReadsList, paraEM, paraPrior, minReads, detail, rescale, calpha, PE, PACKAGE="PING")
}



# ================================================================================
# = Purpose: Screen out abnormal predictions whose parameters are too far from others
# = Input:
# 		ping.df: a data.frame converted from PING result using as(ping,"data.frame") 
# 		alpha: proportion of outliers used to decide filter boundary
# 		detail: indicator of whether or not print/plot detail information
#		deltaB,sigmaB1, sigmaB2,seB: conservative bound, use these boundaries, if calculated boundaries from quantile is narrower. If all estimated parameters are good, we do not want to filter out 5% good nucleosomes
#									 filter "if one sigma above the sigma B1" or "Both sigma F/R above sigmaB2"
#		score: lower bound for scoreF and scoreR, note the default bound is for the case without control, if we have control data, differnt bound should be used.
# = Ouput:
#		ping.df: the filtered dataframe
#		myFilter: the filter  
# ==========
FilterPING <- function(ping.df,detail=F,alpha=0.05,deltaB=c(80,250),sigmaB1=10000,sigmaB2=4900,seB=100,score=0.05,use.min=F)
{
	if (use.min)
	{
		delta.m=max(deltaB[1],quantile(ping.df$delta,alpha/2))
		delta.M=min(deltaB[2],quantile(ping.df$delta,1-alpha/2)) 
		se=c(0,min(seB,quantile(ping.df$se,1-alpha)))
		sigmaSqF=c(0,min(sigmaB1,quantile(ping.df$sigmaSqF,1-alpha)))
		sigmaSqR=c(0,min(sigmaB1,quantile(ping.df$sigmaSqR,1-alpha)))
	}else
	{
		delta.m=min(deltaB[1],quantile(ping.df$delta,alpha/2))
		delta.M=max(deltaB[2],quantile(ping.df$delta,1-alpha/2)) 
		se=c(0,max(seB,quantile(ping.df$se,1-alpha)))
		sigmaSqF=c(0,max(sigmaB1,quantile(ping.df$sigmaSqF,1-alpha)))
		sigmaSqR=c(0,max(sigmaB1,quantile(ping.df$sigmaSqR,1-alpha)))	
	}
	#sigmaSq2=c(0,max(sigmaB2,quantile(ping.df$sigmaSqF,1-alpha),quantile(ping.df$sigmaSqR,1-alpha)))
	sigmaSq2=c(0,sigmaB2)
	delta=c(delta.m,delta.M)
	myFilter=list(score=c(score,Inf),delta=delta,se=se,sigmaSqF=sigmaSqF, sigmaSqR=sigmaSqR,sigmaSq2=sigmaSq2,scoreF=score,scoreR=score)

	if(detail)
	{
		print(myFilter)
		par(mfrow=c(2,3),oma=c(2,2,5,2),mar=c(5,5,2,0))
		plot(density(ping.df$delta),xlim=c(50,300),main="delta")
		abline(v=delta,col=2)
		# plot(density(ping.df$occup),xlim=c(0,quantile(ping.df$occup,0.999)),main="occupancy")
		# abline(v=occup,col=2)
		ping.df$se[ping.df$se>100]=100
		plot(density(ping.df$se),xlim=c(0,min(100,quantile(ping.df$se,0.999))),main="SE")
		abline(v=se,col=2)
		plot(density(ping.df$sigmaSqF),xlim=c(0,max(3600,quantile(ping.df$sigmaSqF,0.999))),main="sigmaSqF")
		abline(v=sigmaSqF,col=2)
		abline(v=sigmaSq2,col=3,lty=2)
		plot(density(ping.df$sigmaSqR),xlim=c(0,max(3600,quantile(ping.df$sigmaSqR,0.999))),main="sigmaSqR")
		abline(v=sigmaSqR,col=2)
		abline(v=sigmaSq2,col=3,lty=2)
		plot(density(ping.df$scoreF),main="scoreF")
		abline(v=score,col=2)	
		plot(density(ping.df$scoreR),main="scoreR")
		abline(v=score,col=2)		
		title("Density and Filter Bound of Estimated Parameters",outer=T)
	}

	ind1=(ping.df$delta<=delta[2])&(ping.df$delta>=delta[1])		
	ind2=ping.df$se<=se[2]
	ind3=(ping.df$sigmaSqF<=sigmaSqF[2])&(ping.df$sigmaSqR<=sigmaSqR[2])
	ind4=(ping.df$sigmaSqF<=sigmaSq2[2])|(ping.df$sigmaSqR<=sigmaSq2[2])
	ind5=(ping.df$scoreF>score)&(ping.df$scoreR>score)
	ind=ind1&ind2&ind3&ind4&ind5
	ping.df=ping.df[ind,]
	return(list(ping.df=ping.df,myFilter=myFilter))
}#temp1=FilterPING(as(PING1,"data.frame"),detail=T)


# 
# #it filter the data.frame converted from ping object
# filterPING3 <- function(ss,filter=list(delta=c(50,250),sigmaSq=22500, se=50, mu=c(0,Inf), chr=NULL))
# {
# 	ind1	<- (ss$delta>=filter$delta[1])&(ss$delta<=filter$delta[2])
# 	ind2	<- (ss$sigmaSqF<filter$sigmaSq)&(ss$sigmaSqR<filter$sigmaSq)
# 	ind3	<- (ss$mu>=filter$mu[1])&(ss$mu<filter$mu[2])
# 	ind4	<- (ss$se<filter$se)
# 	ind4[is.na(ind4)]	<- T  # do not filter by SE if it is not calculatable
# 	ind5	<- rep(T,nrow(ss))
# 	if (length(filter$chr)>0) ind5 <- (ss$chr %in% filter$chr)
# 	ans		<- ss[ind1&ind2&ind3&ind4&ind5,]
# 	return(ans)
# }
# 
# filterPING <- function(ss,delta=c(50,250),sigmaSq=c(0,22500), se=50, mu=c(0,Inf), chr=NULL, score=c(0,Inf))
# {
# 	ind1 <- ss$delta>delta[1] & ss$delta<delta[2]
# 	ind2 <- ss$sigmaSqF>sigmaSq[1] & ss$sigmaSqF<sigmaSq[2]
# 	ind3 <- ss$sigmaSqR>sigmaSq[1] & ss$sigmaSqR<sigmaSq[2]
# 	ind5 <- is.finite(ss$score) & ss$score>score[1] & ss$score<score[2]
# 	ind  <- ind1 & ind2 & ind3 & ind5
# 	if(!is.null(se))	ind<-ind & ss$se>se[1] & ss$se<se[2]
#     if(!is.null(chr))	ind=ind & (ss$chr %in% chr)  
# 	if((mu[1]>0)|is.finite(mu[2]))	ind=ind & ss$mu>mu[1] & ss$mu<mu[2]
# 	ss=ss[ind,]
# 	return(ss)
# }
# 

# make.thickthin2 <- function(ping, Length=73,nSe=3, rescale=NULL, shift=NULL, filter=list(delta=c(0,Inf),se=c(0,Inf),sigmaSqF=c(0,Inf),sigmaSqR=c(0,Inf),score=c(0,Inf)))
# {
# ###input
# # ping: PING output
# # Length: half length of thick bar
# # nSe: used to define nucleosome region
# # rescale, shift: used to rescale the score to [0,1000], if they are missing, we calculate the quantiles for each score and multiply the quantile by 850 than plus 150, so all of the scores are non-linearly scaled to [150,1000]
# # filter: the filters 
# ###output
# # thickthin: a data frame to be visualized in UCSC
# ###note if rescale or shift is not given, user need to rescale the score by themselves
# 	mu<-mu(ping)
#     delta<-delta(ping)
#     se<-se(ping)
#     seF<-seF(ping)
#     seR<-seR(ping)
#     start<-(mu-delta/2-nSe*seF);    
#     end<-(mu+delta/2+nSe*seR);      
# 
# 	#construct temp data frame
# 	temp=data.frame(chr=chromosome(ping), start=round(start), end=round(end), score=score(ping),
# 					thickstart=round(pmax(mu-Length,start)), thickend=round(pmin(mu+Length,end)),
# 					strand=".", se=se(ping), sf=sigmaSqF(ping), sr=sigmaSqR(ping), delta=delta )
# 	
# 	#apply the filters
# 	if(!is.null(filter))
#     {
#       ind1 <- temp$delta>filter$delta[1] & temp$delta<filter$delta[2]
#       ind2 <- temp$sf>filter$sigmaSqF[1] & temp$sf<filter$sigmaSqF[2]
#       ind3 <- temp$sr>filter$sigmaSqR[1] & temp$sr<filter$sigmaSqR[2]
#       ind4 <- is.finite(temp$score) & temp$score>filter$score[1] & temp$score<filter$score[2]
#       ind5 <- is.finite(temp$se) & temp$se>filter$se[1] & temp$se<filter$se[2]
# 	  ind  <- ind1 & ind2 & ind3 & ind4 & ind5
#     }
#     else
#     {
#       ind<-is.finite(score)
#     }
# 	temp=temp[ind,]
# 	temp=temp[order(temp$score,decreasing=T),]
# 	NN=nrow(temp)
# 	temp$name=paste("ping",1:NN,sep="")
# 	temp$chr=as.character(temp$chr)
# 	temp$chr[temp$chr=="MT"]="M"
# 	temp$chr=paste("chr",temp$chr,sep="")
# 
# 
# 	#rescale the score to [0,1000]
# 	if (is.null(rescale)|is.null(shift))
# 	{
# 		temp$score=round(NN:1/NN*850+150)
# 	}else
# 	{
# 		tt=temp$score*rescale+shift; tt=pmin(tt,1000); tt=pmax(tt,100); temp$score=round(tt)
# 	}
#     thickthin = temp[,c("chr","start","end","name","score","strand","thickstart","thickend")]
# 	return(thickthin)
# }

make.thickthin <- function(ping.df, Length=73,nSe=3, rescale=NULL, shift=NULL, 
                           Add.Original.Score=F, roman=F)
{
###input
# ping.df: dataframe converted from PING output
# Length: half length of thick bar
# nSe: used to define nucleosome region
# rescale, shift: used to rescale the rank of score to [150,1000], if they are missing, 
#                 we calculate the quantiles for each score and multiply the quantile by 850
#	  					    then plus 150, so all of the scores are non-linearly scaled to [150,1000]
#	Add.Original.Score: indicator of whether or not add an extra column of PING score in the right
###output
# thickthin: a data frame to be visualized in UCSC
###note if rescale or shift is not given, user need to rescale the score by themselves

	#construct temp data frame 
 	temp=ping.df
  temp$pingScore=temp$score
	temp$start=round(ping.df$mu-ping.df$delta/2-nSe*ping.df$seF)
	temp$end=round(ping.df$mu+ping.df$delta/2+nSe*ping.df$seR)
	temp$thickstart=ping.df$mu-Length
	temp$thickend=ping.df$mu+Length
	temp$thickstart=round(pmax(temp$thickstart,temp$start))
	temp$thickend=round(pmin(temp$thickend,temp$end))
	temp$strand="."
	temp=temp[order(temp$score,decreasing=T),]
	NN=nrow(temp)
	
	#temp$name=paste("ping",1:NN,sep="")
  temp$name=round(temp$pingScore,3)
  if (roman) temp$chr=paste("chr",as.roman(substring(temp$chr,4)),sep="")
  
	#rescale the score to [0,1000]
	if (is.null(rescale)|is.null(shift))
	{
		temp$score=round(NN:1/NN*850+150)
	}else
	{
		tt=temp$score*rescale+shift; tt=pmin(tt,1000); tt=pmax(tt,100); temp$score=round(tt)
	}

  thickthin = temp[,c("chr","start","end","name","score","strand","thickstart","thickend")]
  scores=data.frame(chr=temp$chr,start=round(temp$mu-5),end=round(temp$mu+5),score=temp$pingScore)
  if (Add.Original.Score) thickthin=data.frame(thickthin,pingScore=temp$pingScore)
	return(list(thickthin=thickthin,scores=scores))
}



#filter nucleosomes predicted outside of segment range.
filterPING2 <- function(ss)
{
	ind1	<- (ss$mu<=ss$maxRange)&(ss$mu>=ss$minRange)
	ans		<- ss[ind1,]
	return(ans)
}



#Given a data.frame, and the chrInfo, we truncate certain columns specified by "positions"
# according the chr length, and return the data frame with columns modified
#truncate.result <- function(dat, chr.info.name="../chromInfo.mm9.txt", positions=c("start","end"))
truncateResult <- function(dat, chr.info.name="../chromInfo.mm9.txt", positions=c("start","end"))
{
	chr.info <- read.table(chr.info.name,colClasses=c("character", "integer"))
	names(chr.info)=c("chr","max.pos")
	for(idx.chr in 1:nrow(chr.info))
	{
	  idx <- (dat$chr==chr.info$chr[idx.chr])
		if(!any(idx)) next
		for(idx.pos in 1:length(positions))
		{
			dat[idx, positions[idx.pos]] <- pmin(chr.info$max.pos[idx.chr], dat[idx, positions[idx.pos]])
		}
	}
	return(dat)
}
