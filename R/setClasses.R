## =========================================================================##
## =========================================================================##
##                    Class definitions and contructors                     ##
## =========================================================================##
## =========================================================================##

setClass("ping", 
	contains="pics"
	)

setClass("pingError",
	contains="picsError"
	)

setClass("pingList", 
	representation=representation(PE="logical"),
	contains="picsList"
	)

### Constructor
newPing<-function(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,scoreF,scoreR,Nmerged,converge,range,chr)
{
  if(!all(is.double(w)))
  {
    stop("Argument 'mu' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(mu)))
  {
    stop("Argument 'mu' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(delta)))
  {
    stop("Argument 'delta' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(sigmaSqF)) | !all(is.double(sigmaSqR)))
  {
    stop("Argument 'sigmaSqF/sigmaSqR' must be numeric ", call.=FALSE)
  }
#  if(!all(is.double(seMu)) | !all(is.double(seMuF)) | !all(is.double(seMuR)))
#  {
#    stop("Argument 'seMu/seMuF/seMuR' must be numeric ", call.=FALSE)
#  }
  if(!all(is.double(score)))
  {
    stop("Argument 'score' must be numeric ", call.=FALSE)
  }
  if(!is.numeric(Nmerged))
  {
    stop("Argument 'Nmerged' must be numeric ", call.=FALSE)
  }
  if(!is.logical(converge))
  {
    stop("Argument 'converge' must be logical ", call.=FALSE)
  }
  if(!is.character(chr))
  {
    stop("Argument 'chr' must be a character string", call.=FALSE)
  }
  # if(!all(is.numeric(range)))
  # {
  #   stop("Argument 'range' must be numeric ", call.=FALSE)
  # }  
  new("ping", estimates=list(w=w,mu=mu,delta=delta,sigmaSqF=sigmaSqF,sigmaSqR=sigmaSqR,seMu=seMu,seMuF=seMuF,seMuR=seMuR),converge=converge,score=score,scoreF=scoreF,scoreR=scoreR,Nmerged=Nmerged,range=range,chr=chr)
}

# In case the algorithm does not converge
newPingError<-function(string)
{
  if(!is.character(string))
  {
    stop("Argument 'errorCode' must be of class character", call.=FALSE)
  }
  new("pingError", errorCode=string)
}

newPingList<-function(List, paraEM, paraPrior, minReads, N, Nc, PE)
{
  if(!is.list(paraEM) & !all(sapply(paraEM,"is.numeric")))
  {
    stop("Argument 'paraEM' must be a list of numeric arguments", call.=FALSE)
  }
  if(!is.list(paraPrior) & !all(sapply(paraPrior,"is.numeric")))
  {
    stop("Argument 'paraPrior' must be a list of numeric arguments", call.=FALSE)
  }
  if(!is.list(minReads) & !all(sapply(minReads,"is.numeric")))
  {
    stop("Argument 'minReads' must be a list of numeric arguments", call.=FALSE)
  }
  if(!all((lapply(List,"class")=="ping" | lapply(List,"class")=="pingError")))
  {
    stop("Argument 'List' must be a list of 'ping' or 'pingError' arguments", call.=FALSE)
  }
  if(!is.integer(N) | !is.integer(Nc))
  {
    stop("Argument 'N' and 'Nc' must be integers", call.=FALSE)    
  }
  PE<-as.logical(PE)
  if(!is.logical(PE))
  {
    stop("Argument 'PE' must be logical", call.=FALSE)    
  }
  new("pingList", List=List, paraEM=paraEM, paraPrior=paraPrior, minReads=minReads, N=N, Nc=Nc, PE=PE)
}
