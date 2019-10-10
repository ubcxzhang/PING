segmentPING<-function(data, dataC=NULL, map=NULL,
		minReads=2, minReadsInRegion=3, jitter=FALSE, maxLregion=1200, minLregion=80, step=NULL, width=NULL, #dataType="H",
		islandDepth=3, min_cut=50, max_cut=1000, maxReadsWidth=500, PE=FALSE )
{
  ##Determine PE/SE based on reads width
  #if(length(unique(width(data)))==1)
    #PE<-FALSE
  #else
    #PE<-TRUE

  if(!isTRUE(PE))
  {
    cat("Performing segmentation for single end reads\n")
    step<-2L
    width<-150L
	  #if(dataType=="TF"){
		  #if(is.null(step)){
			  #step<-20L
	          #}
		  #if(is.null(width)){
			  #width<-150L
		  #}
	  #}
	  #if(dataType=="H"){
		  #if(is.null(step)){
			  #step<-2L
		  #}
		  #if(is.null(width)){
			  #width<-150L
		  #}
	  #}
	  newSet<-segReadsGeneric(data, dataC=dataC, map=map, minReads=minReads, minReadsInRegion=minReadsInRegion, jitter=jitter, maxLregion=maxLregion,minLregion=minLregion, step=step, width=width, package="PING")
	  
	  
  }
  else
  {
	  cat("Performing segmentation for paired-end reads\n")
	  if(var(width(data))<5)
    		warning("This data seems to be Single-End reads")
	  
	  chr<-seqlevels(data)
	  #if(!is.character(chr))
		  #stop("Argument chr should be a character\n")
	  if(length(chr)>1)
		  stop("Paired-end sequencing data segmentation does not support multiple chromosomes\n") 
	  if(!is.numeric(islandDepth))
		  stop("Argument islandDepth should be an integer\n")
	  if(!is.numeric(min_cut) | !is.numeric(max_cut))
		  stop("Arguments min_cut and max_cut should be integers, provided: ", class(min_cut),", ", class(max_cut), "\n")
	  
	  #with GRanges as input
	  #data<-data[seqnames(data)==chr]
          data<-data[width(data)<maxReadsWidth,]
	  #Save xi for later
	  xi<-mean(width(data))
	  PE.RD<-IRanges(start=start(data), end=end(data))
	  PE.RD<-PE.RD[which(width(PE.RD)>min_cut | width(PE.RD)<max_cut)]		
	  PE.RD<-IRanges(start=start(PE.RD), end=end(PE.RD))
	  #Get the candidate regions (IRanges)
	  candidate_RD<-candidate.region(PE.RD, islandDepth, min_cut, max_cut)
	  

	  #Getting parameters for segmentation
	  #PEMF.RD <- IRanges(start=data$yFm$"pos.+", end=data$yFm$"pos.+"+1)
	  #PEMR.RD <- IRanges(start=data$yRm$"pos.-"-1, end=data$yRm$"pos.-")
	  PEMF.RD<-IRanges()
	  PEMR.RD<-IRanges()

	  #TODO: Add the treatment in case !is.null(dataC)
	  N   <- length(PE.RD)
	  NFm <- length(PEMF.RD)
	  NRm <- length(PEMR.RD)
	  Nc <- NcFm <- NcRm <- as.integer(0)#integer(0)
	  
	  #Perform the segmentation
	  if(is.null(map)){
		  segmentation <- segChrRead(candidate_RD, PE.RD, PEMF.RD, PEMR.RD, PEC.RD=NULL, PECMF.RD=NULL, 
				  PECMR.RD=NULL, map.Start=NULL, map.End=NULL, chr=chr)
	  }else{ 
		  segmentation <- segChrRead(candidate_RD, PE.RD, PEMF.RD, PEMR.RD, PEC.RD=NULL, PECMF.RD=NULL, PECMR.RD=NULL, 
				  map.Start=start(map[chr]), map.End=end(map[chr]), chr=chr)
	  }
	  paraSW<-list(islandDepth=islandDepth, min_cut=min_cut, max_cut=max_cut, xi=xi)
	  newSet<-segReadsListPE(segmentation, paraSW=paraSW, N=N, NFm, NRm, Nc, NcFm, NcRm)
  }
  

  return(newSet)
}


