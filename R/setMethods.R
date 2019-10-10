
setAs("pingList", "RangedData",
      function(from) {
            makeRangedDataOutput(from, type="bed", filter=list(delta=c(50,300),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500),score=c(1,Inf)),length=100)
      })

setAs("pingList", "data.frame",
  function(from)
  {
    ans<-data.frame(ID=rep(1:length(from),K(from)),chr=chromosome(from),w=w(from), mu=mu(from),
                    delta=delta(from), sigmaSqF=sigmaSqF(from), sigmaSqR=sigmaSqR(from),se=se(from),
                    score=score(from), scoreF=scoreForward(from),scoreR=scoreReverse(from),
                    minRange=minRange(from), maxRange=maxRange(from), seF=seF(from), seR=seR(from) )
    ans$chr       <- as.character(ans$chr)
    ans           <- ans[is.finite(ans$mu),]
    return(ans)
  })

setMethod("as.data.frame", "pingList", function(x, ...){as(x, "data.frame")})

setMethod("show", "pingList",
          function(object)
          {
            cat("Object of class 'pingList'","\n")
            cat("This object has the following slots: \n")
            cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
            cat("List is a list of 'ping' or pingError ojects\n")
          })

# setGeneric("density", function(x, ...) standardGeneric("density"))
setMethod("density", "ping",
          function(x,strand="+",step=10,sum=FALSE,filter=NULL,scale=TRUE)
          {
            
            # Check that all filters are passed
            missingNames<-!c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")%in%names(filter)
            filter[c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")[missingNames]]<-list(c(0,Inf))
            if(strand=="+")
            {
              strand<-1
            }
            else if(strand=="-")
            {
              strand<--1
            }
            else if(strand=="*")
            {
              strand<-0
            }
            else
            {
              stop("Strand must be either '+', '-' or '*'")
            }            
            strand<-as.double(paste(strand,"1",sep=""))
            ans<-.Call("getDensity", x, strand, step, filter, sum, scale, PACKAGE="PICS")
            return(ans);
          }
)

setMethod("density", "pingList",
          function(x,strand="+",step=10,sum=FALSE,filter=NULL,scale=TRUE)
          {
            # Check that all filters are passed
            missingNames<-!c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")%in%names(filter)
            filter[c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")[missingNames]]<-list(c(0,Inf))

            if(strand=="+")
            {
              strand<-1
            }
            else if(strand=="-")
            {
              strand<--1
            }
            else if(strand=="*")
            {
              strand<-0
            }
            else
            {
              stop("Strand must be either '+', '-' or '*'")
            }
            ans<-.Call("getDensityList", x, strand, step, filter, sum, scale, PACKAGE="PICS")
            return(ans);
          }
)

setMethod("density", "pingError",
          function(x,strand=NULL,step=NULL,sum=NULL,filter=NULL)
          {
            return(NULL)
          }
)


setMethod("summary", "segReads",
      function(object)
      {
        m<-min(object@yF[1],object@yR[1])
        M<-max(tail(object@yF,1),tail(object@yR,1))        
        cat("** Region summary ** \n")                    
        cat("Summary on Forward reads:\n",summary(object@yF,digits=100),"\n")
        cat("Summary on Reverse reads:\n",summary(object@yR,digits=100),"\n")
        cat("Summary on control Forward reads:\n",summary(object@cF,digits=100),"\n")
        cat("Summary on control Reverse reads:\n",summary(object@cR,digits=100),"\n")        
        cat("Non mappable intervals cover ", sum(diff(t(object@map)))/(M-m),"% of the region \n")
      })



setMethod("[","pingList",
		function(x,i, j,..., drop=FALSE)
		{
			if(missing(i))
			{
				return(x)
			}
			if(!missing(j))
			{
			  stop("incorrect number of dimensions")
			}
      else
      {
        newPingList(x@List[i], x@paraEM, x@paraPrior, x@minReads, x@N, x@Nc)        
      }
		})


setMethod("summary", "pingList",
          function(object)
          {
            cat("** Experiment information ** \n")
            cat("Chromosomes interogated: ")
            cat(unique(chromosome(object)),"\n")
            cat("Number of reads:")
            cat("In IP: ",object@N," in control: ",object@Nc,"\n")
            cat("** Prior parameters ** \n")
            cat("The following settings were used:\n")          
            cat("  Hyper parameters for the fragment length distribution:\n")
            cat("  xi, rho, alpha, beta: ", object@paraPrior$xi,",", object@paraPrior$rho, ",", object@paraPrior$alpha, ",", object@paraPrior$beta,"\n")          
            cat("** Score summary ** \n")                    
            print(summary(score(object)))
            cat("** Fragment length distribution summary ** \n")                    
            print(summary(delta(object)))
            cat("** Summary on the number of binding events per candidate region** \n")
            summary(K(object))
      })


setMethod("summary", "ping",
          function(object)
          {
            cat("** Score ** \n")                    
            cat(score(object),"\n")
            cat("** Fragment length estimate ** \n")                    
            cat(delta(object),"\n")
            cat("** Number of binding events in the candidate region** \n")
            cat(K(object),"\n")
          })



setMethod("plot", signature("ping", "segReads"),
          function(x, y, addKernel=FALSE, addNucleosome=FALSE, addSe=TRUE, main=NULL, rescale=1, ...)
      {
         #Set outer and figure margins to reduce gap between plots
        if(addNucleosome)
        {
          nG<-3
          par(oma=c(2.5,5,5,5),mar=c(0,5,0,0),cex.lab=2)
          layout(matrix(1:nG,ncol=1), heights = c(.4,.5,.1))
        }
        else
        {
          nG<-2
          par(oma=c(2.5,5,5,5),mar=c(0,5,0,0),cex.lab=2)
          layout(matrix(1:nG,ncol=1), heights = c(.5,.5))
        }
        
        step<-5/rescale
        .densityMix<-function(x,para)
        {
          v<-4
          w<-para$w
          mu<-para$mu
          sigmaSq<-para$sigmaSq
          sigma<-sqrt(sigmaSq)
          xNorm<-outer(-mu,x,"+")/sigma 
          return(colSums(w*dt(xNorm,df=v)/sigma))
        }

        yF<-y@yF
        yR<-y@yR
        cF<-y@cF
        cR<-y@cR        
        map<-y@map
        m<-min(yF[1],yR[1])-100/rescale
        M<-max(tail(yF,1),tail(yR,1))+100/rescale

        paraR<-list(w=x@estimates$w, mu=x@estimates$mu+x@estimates$delta/2, sigmaSq=x@estimates$sigmaSqR)
        paraF<-list(w=x@estimates$w, mu=x@estimates$mu-x@estimates$delta/2, sigmaSq=x@estimates$sigmaSqF)

        dR<-.densityMix(seq(m,M,step),paraR)
        dF<-.densityMix(seq(m,M,step),paraF)
        maxRange<-max(c(dF,dR))
        plot(seq(m,M,step),dF,xlim=c(m,M),ylim=c(0,maxRange),lty=2,type="l",xlab="",ylab="density",xaxt='n',axes=FALSE)
        title(main=main,outer=TRUE,cex.main=2)
        axis(2)
        axis(1)
        
        
        lines(seq(m,M,step),dR,lty=2,col=2)

        # if(length(map)>0)
        # {
        #   nMap<-nrow(map)
        #   for(i in 1:nMap)
        #   {
        #     segments(map[i,1], 0, map[i,2], 0,lwd=3,col=3)
        #   }
        # }
        
        # Add kernel density estimate
        if((addKernel==TRUE) & (length(yF)>1 & length(yR)>1))
        {
          dkF<-density(yF)
          dkR<-density(yR)
          lines(dkF,lty=3)
          lines(dkR,col=2,lty=3)
        }

        #Add single components and se's
        K<-length(x@estimates$w)
        for(k in 1:K)
        {
          paraR<-list(w=x@estimates$w[k], mu=x@estimates$mu[k]+x@estimates$delta[k]/2, sigmaSq=x@estimates$sigmaSqR[k])
          paraF<-list(w=x@estimates$w[k], mu=x@estimates$mu[k]-x@estimates$delta[k]/2, sigmaSq=x@estimates$sigmaSqF[k])

          dsR<-.densityMix(seq(m,M,step),paraR)
          dsF<-.densityMix(seq(m,M,step),paraF)

          lines(seq(m,M,step),dsF,lty=1)
          lines(seq(m,M,step),dsR,col=2,lty=1)
        }

        stripchart(yF[1],pch=">",method="stack",cex=2,at=.55,add=FALSE,axes=FALSE,xlim=c(m,M),ylim=c(0,1))        
        if(length(map)>0)
        {
          nMap<-nrow(map)
          symbols((map[,1]+map[,2])/2,rep(.35,nMap),rectangle=cbind(map[,2]-map[,1],rep(.6,nMap)), inches=FALSE, bg=grey(.6), fg=0, add=TRUE,xlim=c(m,M),ylim=c(0,1))
        }

        stripchart(yF,pch=">",method="stack",cex=2,at=.55,axes=FALSE,xlim=c(m,M),ylim=c(0,1),add=TRUE)
        mtext("IP",cex=1.2,side=2,las=2,at=.45)
        stripchart(yR,pch="<",method="stack",cex=2,at=.35,col=2,add=TRUE,offset=-1/3)
        
        abline(h=.45,lty=1,col="grey")
        if(addSe)
        {
          points(x@estimates$mu,rep(.45,K),pch="+",cex=2)
          if (any(x@estimates$seMu!=0))
          {
            points(x@estimates$mu-2*x@estimates$seMu,rep(.45,K),pch="[",cex=1)
            points(x@estimates$mu+2*x@estimates$seMu,rep(.45,K),pch="]",cex=1)
            segments(x@estimates$mu-2*x@estimates$seMu,rep(.45,K),x@estimates$mu+2*x@estimates$seMu,rep(.45,K),lwd=1,lty=rep(1,K))
          }
        }

        if(length(cF)>0)
        {
          stripchart(cF,pch=">",method="stack",at=0.15,cex=2,add=TRUE,xlim=c(m,M),ylab="Cont.",axes=FALSE)
          abline(h=.1,lty=1,col="grey")
        }
        if(length(cR)>0)
        {
          stripchart(cR,pch="<",method="stack",at=0.05,cex=2,col=2,add=TRUE,offset=-1/3)
        }
        if((length(cR)==0) & (length(cF)==0))
        {
        }
        mtext("Cont.",cex=1.2,side=2,las=2,at=.1)
        
        if(addNucleosome)
        {
          plot(c(m,M),c(0,1),axes=FALSE,col=0,ylim=c(0,1),xlim=c(m,M),ylab="")
          seq<-c(seq(-3,0,.1),seq(3,0,-.1))
          sapply(seq,function(shift,x,m,M,K){symbols(x@estimates$mu+shift*se(x),rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, bg=0, fg=grey(abs(shift)*se(x)/(3*(se(x)))), add=TRUE,xlim=c(m,M),ylim=c(0,1),lwd=2)},x,m,M,K)
          symbols(x@estimates$mu,rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, bg="white", fg=grey(abs(0)), add=TRUE,xlim=c(m,M),ylim=c(0,1))
          
          # if(addSe)
          # {
          #   symbols(x@estimates$mu,rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, bg=grey(.5*pmin(se(x)/50,1)), fg=0, add=TRUE,xlim=c(m,M),ylim=c(0,1))
          # }
          # else
          # {
          #   
          #   symbols(x@estimates$mu,rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, fg=0, bg=1, add=TRUE,xlim=c(m,M),ylim=c(0,1))
          # }
          mtext("Nucl.",cex=1.2,side=2,las=2,at=.5)
        }
})

setMethod("plot", signature("pingError", "segReads"),
          function(x, y, addKernel=FALSE, main=NULL, rescale=1,  ...)
      {
        par(oma=c(2.5,5,5,5),mar=c(0,5,0,0))
        layout(matrix(1:2,ncol=1), heights = c(.2,.1))
        
        yF<-y@yF
        yR<-y@yR
        cF<-y@cF
        cR<-y@cR
        map<-y@map
        m<-min(yF[1],yR[1])-100/rescale
        M<-max(tail(yF,1),tail(yR,1))+100/rescale

        stripchart(yF,pch=">",method="stack",cex=2,at=.5,add=FALSE,axes=FALSE,xlim=c(m,M),ylab="Cont | Inp.",ylim=c(0,1))
        stripchart(yR,pch="<",method="stack",cex=2,at=.5,col=2,add=TRUE,offset=-1/3)
        title(main=main,outer=TRUE,cex.main=2)
        
        abline(h=.35,lty=3)

        # Add kernel density estimate
        if((addKernel==TRUE) & (length(yF)>1 & length(yR)>1))
        {
          dkF<-density(yF,bw=75/rescale)
          dkR<-density(yR,bw=75/rescale)
          plot(dkF,lty=3)
          lines(dkR,col=2,lty=3)
        }

        if(length(cF)>0)
        {
          stripchart(cF,pch=">",method="stack",at=0.2,cex=2,add=TRUE,xlim=c(m,M),ylab="Cont.",axes=FALSE)
        }
        if(length(cR)>0)
        {
          stripchart(cR,pch="<",method="stack",at=0.2,cex=2,col=2,add=TRUE)
        }
})

setMethod("plot", signature("pingList", "segReadsList"),
          function(x, y, regionIndex=NULL, addKernel=FALSE, addNucleosome=FALSE, addSe=TRUE,main=NULL, rescale=1, ...)
{
  setMain<-is.null(main)
  if(is.null(regionIndex))
  {
    regionIndex<-1:length(x@List)
  }
  for(i in regionIndex)
  {    
    if(setMain)
    {
      main<-paste(as.character(i)," (",y@List[[i]]@chr,")",sep="")
    }
    if(class(x@List[[i]])!="pingError")
    {
      plot(x@List[[i]],y@List[[i]],addKernel=addKernel, addNucleosome=addNucleosome, addSe=addSe,main=main, rescale=rescale, ...)
    }
    else
    {
      plot(x@List[[i]],y@List[[i]],addKernel=addKernel, main=paste(as.character(i)," (",y@List[[i]]@chr,")",sep=""), rescale=rescale, ...)
      warning("Object of class pingError, no PING density displayed")
    }
  }
})

setMethod("plot", signature("pingList", "pingList"),
          function(x, y, filter=NULL, h=.1, ...)
{
  FDR<-pingFDR(x,y,filter=filter)
  plot(FDR[,2],FDR[,1],xlab="score",ylab="FDR",panel.first=grid(nx=50),...)
  # points(FDR[,2],FDR[,3]/max(FDR[,3]),xaxt="n",yaxt="n",lty=3,col=3,pch=2)
  # axis(4,at=seq(0,1,.05),labels=max(FDR[,3])*seq(0,1,.05))
  FDRex<-FDR[FDR[,1]>0,]
  notDup<-rev(!duplicated(rev(FDRex[,1])))
  lines(FDRex[notDup,2],FDRex[notDup,1],col=2,lty=2,lwd=1.5)
  abline(h=h,lw=1.5,col="grey")
})

# plot function for two ping results in data.frame format
setMethod("plot", signature("data.frame", "data.frame"), 
          function(x, y, h=.1, logscale=F, ...)
{
  FDR<-pingFDR2(x,y)
  if(logscale) {FDR$score=log(FDR$score); xlab="log(score)"} else xlab="score"
  plot(FDR[,"score"],FDR[,"FDR"],xlab=xlab,ylab="FDR",panel.first=grid(nx=50), 
  		 ylim=range(tail(head(FDR$FDR,-1),-1)), ...)
  FDRex<-FDR[FDR[,"FDR"] > 0,]
  notDup<-rev(!duplicated(rev(FDRex[,"FDR"])))
  lines(FDRex[notDup,"score"],FDRex[notDup,"FDR"],col=2,lty=2,lwd=1.5)
  lines(FDR[,"score"],FDR[,"FDR"],lty=1,lwd=1.5)
  abline(h=h,lw=1.5,col="grey")
})



############################################################################################



## CoverageTrack
CoverageTrack<-function(ping, reads, chr, gen="gen", FragmentLength=200, name="XSET")#, PE=FALSE
{
	PE<-ping@PE
	if(class(reads)!="GRanges")
	{
		stop("The reads should be of class 'GRanges'")
	}
	else if(!isTRUE(PE))
	{#no need to resize for PE sequencing data (width known)
		EXT<-resize(reads[seqnames(reads)==chr], width=FragmentLength)
		EXT<-EXT[start(EXT)>0]
	}
	else
	{
		EXT<-reads[seqnames(reads)==chr]
	}
	XSET = coverage(EXT)
	covTrack<-DataTrack(data=XSET[[chr]]@values, start=start(XSET[[chr]]), width=width(XSET[[chr]]),
	col="black",
	chromosome=chr, genome=gen, type=c("g","s"), v=0, col.axis="black", cex.axis=1, col.title="black", col.grid="gray", name=name)
	return(covTrack)
}


##
# RawReadsTrack
#   INPUT: The reads used in the segmentation step
#   OUTPUT: An AnnotationTrack object showing the starting position of forward and reverse reads
##
RawReadsTrack<-function(ping, reads, chr, gen="gen", from=NULL, to=NULL,...)# PE=FALSE, 
{
	PE<-ping@PE
	reads<-reads[seqnames(reads)==chr]
	if(length(c(from,to))==2 | is.numeric(c(from, to)))
	{
		reads<-reads[start(reads)>from & end(reads)<to]
	}
	if(!isTRUE(PE))
	{
		idxF<-which(as.character(strand(reads))=="+")
		idxR<-which(as.character(strand(reads))=="-")
		coordsF<-start(reads[idxF]) 
		coordsR<-end(reads[idxR])
	}
	else
	{
		coordsF<-start(reads)
		coordsR<-end(reads)
	}
	pos<-unique(c(coordsF, coordsR))

	val<-vector("list", length(pos))
	s1<-system.time({	
	for(idx in 1:length(pos))
	{
		val[[idx]]<-c(length(which(coordsF==pos[[idx]])), -length(which(coordsR==pos[[idx]])))# get m$M
	}
	})#end s1
	m<-min(unlist((val)), -1)
	M<-max(unlist((val)), 1)
	#add the NA
	s2<-system.time({
	vec<-unlist(lapply(val, function(x){
						if(x[[1]]==0)
						{
							c(rep(NA, abs(x[[2]]-m)), seq(x[[2]], -1), rep(NA, M))
						}
						else if(x[[2]]==0)
						{
							c(rep(NA, abs(m)), seq(1,x[[1]]), rep(NA, M-x[[1]]))
						}
						else
						{
							c(rep(NA, abs(x[[2]]-m)), seq(x[[2]], -1), seq(1,x[[1]]), rep(NA, M-x[[1]]))
						}
					}))
	})#end s2
	s3<-system.time({
	mat<-matrix(vec, nrow=M+abs(m), ncol=length(pos), dimnames=list(c(seq(m,-1),seq(1,M)) ,pos))
	})#end s3
	s4<-system.time({
			RawReadsTrack<-DataTrack(data=mat, start=pos, end=pos, chromosome=chr, genome=gen,
			groups=rep(c("rev", "fwd"),c(abs(m),M)),
			col=c("black","black"), pch=c(">", "<"), font="sans",
			ylim=c(m, M), size=1,
			showAxis=FALSE, col.axis="transparent", col.title="black", ...)	
	})#end s4
return(RawReadsTrack)
}





##
# NucleosomeTrack
#   Shows nucleosome positioning prediction with standard error
##
NucleosomeTrack<-function(PS, chr, gen="gen", scoreTrack=TRUE, scoreThreshold=0.05, name="PING", ...) #from=NULL, to=NULL)
{
	if(!is.null(scoreThreshold))
		PS<-FilterPING(PS, score=scoreThreshold)$ping.df

	#If object is ping, coerce into df
	if(class(PS)=="pingList")
		PS<-as(PS, "data.frame")
	#Order the data.frame by start
	PS<-PS[with(PS, order(mu)),]	
	#wide nucleosomes are nucleosomes +se
	smallWidth<-rep(147, nrow(PS))
	wideWidth<-147+2*PS$se
	smallStart<-PS$mu-(smallWidth/2)
	wideStart<-PS$mu-(wideWidth/2)

	#Drawing the wide first is important
	starts<-c(wideStart, smallStart)
	widths<-c(wideWidth, smallWidth)

	idList<-c(rep("", length(wideStart)), PS$ID)

	tList<-list()
	if(isTRUE(scoreTrack))
	{
		NSTrack<-DataTrack(data=score(PS), start=mu(PS)-5, width=10, chromosome=chr, genome=gen,
			type=c("histogram","g"), size=1,
			col.axis="black", cex.axis=0.5, col.title="black", name="score", v=0, col.grid="gray", ...)
		tList<-c(tList, NSTrack)
	}	
	nucTrack<-AnnotationTrack(start=starts, 
			width=widths,
			chromosome=chr, genome=gen,
			lwd=1, col="darkgray", alpha=1, size=0.5,# featureID=idList, showFeatureId=TRUE,
			stacking="dense",feature=c(rep("wide", nrow(PS)), rep("small", nrow(PS))),
			wide="darkgray",small="white", 
			shape="ellipse", col.title="black", collapse=FALSE, name=name,
			...) #collapse=FALSE is important. Prevents Gviz from changing the element order during prepare mode.
	tList<-c(tList, nucTrack)
	return(tList)
}


##
# Plot a summary of PING estimates using Gviz
#   PS : The ouput of PING or postPING or a list.
#   reads : The reads used for the segmentation (GRanges)
#   chr : The chromosome to display
#   from,to : The range to display
#   title : A main title for the plot
#   FragmentLength : Length of the DNA fragments used. For CoverageTrack.
#   scoreThreshold : Nucleosomes with a score < scoreThreshold are not displayed on the plot.
#   GRT : Doesn't  work in the current version of Gviz.
#   PE : Set to TRUE for Paired-End Sequencing data.
##
plotSummary<-function(PS, ping, reads, chr, gen="gen", from=NULL, to=NULL, FragmentLength=200, title="", scoreThreshold=0.05)
{
	GRT<-FALSE
	PE<-ping@PE 

	if(class(PS)!="list")
		PS<-list(PS)
	gt<-GenomeAxisTrack(add53=TRUE, add35=TRUE, col="black")
	tList<-list(gt)
	
	if(isTRUE(GRT))
	{	
		bgrTrack<-tryCatch(expr={
					r<-BiomartGeneRegionTrack(chromosome=chr, genome=gen, start=from, end=to, col.title="black", name="Exons", size=0.5, showId=TRUE)
				}, error=function(err){
					print(paste("WARNING: Unable to create BiomartGeneRegionTrack:", err))
					print("The 'gen' and 'chr' argument should be named as in the UCSC website")
					r<-NULL
				})# end tryCatch. r is returned (last assigned value)
		tList<-c(tList, bgrTrack)
	}
			
	covTrack<-CoverageTrack(ping=ping, reads=reads, chr=chr, gen=gen, FragmentLength=FragmentLength)
	tList<-c(tList, covTrack)
	if(!isTRUE(PE))
	{
		rrTrack<-RawReadsTrack(ping=ping, reads=reads, chr=chr, gen=gen, from=from, to=to, name="Aligned reads")
		tList<-c(tList, rrTrack)
	}
	for(idxPS in 1:length(PS))
	{
		tList<-c(tList,
			NucleosomeTrack(PS=PS[[idxPS]], chr=chr, gen=gen, scoreThreshold=scoreThreshold)
			)
	}
	plotTitle<-paste(title,chr,":",from,"-",to,"(",to-from,"bps)", sep="")
	# Plotting
	plotTracks(trackList=tList, from=from, to=to, main=plotTitle)
	#Return the tracks so they can be added to other plots. invisible=no print()
	return(invisible(tList))
}
