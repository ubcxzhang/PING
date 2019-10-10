\name{pingList-class}
\docType{class}
\alias{pingList}
\alias{pingList-class}
\alias{chromosome,pingList-method}
\alias{mu,pingList-method}
\alias{delta,pingList-method}
\alias{w,pingList-method}
\alias{score,pingList-method}
\alias{density,pingList-method}
\alias{se,pingList-method}
\alias{seF,pingList-method}
\alias{seR,pingList-method}
\alias{[,pingList,ANY,ANY-method}
\alias{[[,pingList,ANY,ANY-method}
\alias{newPingList,pingList-method}
\alias{as.data.frame,pingList-method}
\alias{length}
\alias{length,pingList-method}
\alias{newPingList}
\title{The ping class}

\description{
This object is used to gather all parameters from fitting PING to multiple candidate regions (as returned by the 'segmentReads' function). 
The objet contains the following slots: 'List', 'paraPrior', 'paraEM', 'minReads', 'N', 'Nc'. 'List' is a list of 'ping' or 'pingError' objects. 'paraPrior' is a list containing the hyperparameters used for the prior, 'paraEM' is a list of convergence parameters for the EM, 'minReads' is a list containing the minimum number of reads used to fit a region with 'PING', 'N' is the total number of reads in the ChIP samples while 'Nc' is the total number of reads in the control sample.
}

\section{Accessors}{
The PING package provide accessors to directly access to most of the parameters/standard errors and chromosomes. 
In the code snippets below, 'x' is a 'pingList' object. For all accessors, the 'pingError' objects are omitted, so that the accessors only return values for the 'ping' objects (i.e. all valid binding events).
\describe{
\item{'chromosome(x)'}{Gets the chromosome names of all candidate regions.}
\item{'mu(x)'}{Gets the position estimates of all binding sites identified in all candidate regions.}
\item{'delta(x)'}{Gets the average fragment lengths of all binding sites identified in all candidate regions.}
\item{'sigmaSqF(x)'}{Gets the F peak variances of all binding sites identified in all candidate regions.}
\item{'sigmaSqR(x)'}{Gets the R peak variances of all binding sites identified in all candidate regions.}
\item{'seF(x)'}{Gets the standard errors of all binding site position estimates identified in all candidate regions.}
\item{'seF(x)'}{Gets the standard errors of all F peak modes identified in all candidate regions.}
\item{'seR(x)'}{Gets the standard errors of all R peak modes identified in all candidate regions.}
\item{'score(x)'}{Gets the scores of all binding events identified in all candidate regions.}
}
}

\section{Constructor}{
newPingList(List, paraEM, paraPrior, minReads, N, Nc)
\describe{
\item{List}{The mixture weights (a vector)}
\item{paraEM}{The binding site positions (a vector)}
\item{paraPrior}{The DNA fragment lengths (a vector)}
\item{N}{The variance parameters for the forward distribution (vector)}
\item{Nc}{The variance parameters for the forward distribution (vector)}
}
}


\section{Methods}{
  \describe{
 \item{[}{\code{signature(x = "ping")}: subset PING object.}
}
}


\section{Methods}{
  \describe{
 \item{length}{\code{signature(x = "ping")}: subset PING object.}
}
  \describe{
    \item{as.data.frame}{\code{signature(x = "pingList")}: Coerce a \code{pingList} to a \code{data.frame}}.
  }
}

\arguments{
		\item{object}{An object of class \code{ping}.}	
}


\section{Constructor}{
newPingList<-function(List, paraEM, paraPrior, minReads, N, Nc) constructs a new 'pingList' object with the following arguments. 
\describe{
\item{newPingList}{}
\item{w}{The mixture weights (a vector)}
\item{mu}{The binding site positions (a vector)}
\item{delta}{The DNA fragment lengths (a vector)}
\item{sigmaSqF}{The variance parameters for the forward distribution (vector)}
\item{sigmaSqR}{The variance parameters for the reverse distribution (vector)}
\item{seMu}{The standard errors for mu (vector)}
\item{seMuF}{The standard errors for muF (vector)}
\item{seMuR}{The standard errors for muR (vector)}
\item{seMuR}{The standard errors for muR (vector)}
\item{score}{The scores for each binding event (vector)}
\item{Nmerged}{The number of peaks that were merged (integer)}
\item{converge}{A logical value, TRUE, if the EM as converged}
\item{infMat}{The information matrix}
\item{chr}{The chromosome for the region}
}
}

\author{
  Xuekui Zhang <\email{xzhang@stat.ubc.ca}>, Sangsoon Woo, \email{swoo@fhcrc.org}
 and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>
}

\references{
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" PlosONE, under review.
}

\seealso{
\code{\link{ping}}
}


\examples{
# Here is an example of how to construct such a region
# Typically, you would not do this manually, you would use the ping function to return a 'pingList'
# that contains a list of 'ping' or 'pingError' object.
w<-1
mu<-10000
delta<-150
sigmaSqF<-5000
sigmaSqR<-5000
seMu<-10
seMuF<-10
seMuR<-10
score<-5
Nmerged<-0
converge<-TRUE
infMat<-matrix(0)
chr<-"chr1"
range<-c(1000,2000)
# Contructor
#myPING1<-newPing(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,Nmerged,converge,infMat,as.integer(range),chr)
#myPING2<-newPing(w,mu+1000,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,Nmerged,converge,infMat,as.integer(range),chr)

#minReads<-list(perPeak=2,perRegion=5)
#paraPrior<-list(xi=200,rho=1,alpha=20,beta=40000)
#paraEM<-list(minK=1,maxK=15,tol=10e-6,B=100)
#N<-100
#Nc<-200

#mynewPingList<-newPingList(list(myPING1,myPING2), paraEM, paraPrior, minReads, as.integer(100), as.integer(200))
# Accessors
# Get the standard error of Mu
#se(mynewPingList)
# Get the standard error of MuF
#seF(mynewPingList)
# Get the scores
#score(mynewPingList)
}


\keyword{models}
