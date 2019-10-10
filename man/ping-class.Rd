\name{ping-class}
\docType{class}
\alias{newPing}
\alias{chromosome,ping-method}
\alias{mu,ping-method}
\alias{delta,ping-method}
\alias{w,ping-method}
\alias{score,ping-method}
\alias{se,ping-method}
\alias{seF,ping-method}
\alias{seR,ping-method}
\alias{newPing,ping-method}
\alias{seF}
\alias{seR}
\alias{se}
\alias{w}
\alias{delta}
\alias{chromosome}
\alias{mu}

\title{The ping class}

\description{
This object is used to gather all parameters from fitting PING to a single candidate region. 
The objet contains the following slots: 'estimates', 'infMat', 'Nmerged', 'converge', 'chr'. 'estimates' is a list containing all parameters estimates as well as standard errors. 'infMat' is the Cholesky decomposition of the information matrix, 'converge' is a logical value indicating whether the EM algorithm has converged, while 'chr' is a character string corresponding to a candidate region's chromosome. 'Nmerged' gives the number of binding events that were merged; binding events that overlap are merged (see the cited paper below for details).
}

\section{Accessors}{
The PING package provide accessors to directly access to most of the parameters/standard errors and chromosome. In the code snippets below, 'x' is a 'ping' object.
\describe{
\item{`chromosome(x)'}{Gets the chromosome name of the candidate region.}
\item{`mu(x)'}{Gets the position estimates of all binding sites identified in the region.}
\item{`delta(x)'}{Gets the average fragment lengths of all binding sites identified in the region.}
\item{`sigmaSqF(x)'}{Gets the F peak variances of all binding sites identified in the region.}
\item{`sigmaSqR(x)'}{Gets the R peak variances of all binding sites identified in the region.}
\item{`se(x)'}{Gets the standard errors of all binding site position estimates identified in the region.}
\item{`seF(x)'}{Gets the standard errors of all F peak modes identified in the region.}
\item{`seR(x)'}{Gets the standard errors of all R peak modes identified in the region.}
\item{score}{\code{signature(x = "ping")}: return the score for each binding event.} 
\item{scoreF}{\code{signature(x = "ping")}: return the score  of the forward (F)  for each binding event.}
\item{scoreR}{\code{signature(x = "ping")}: return the score  of the forward (R)  for each binding event.}
}
}

\section{Constructor}{
newPing(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,Nmerged,converge,infMat,chr) construct a new 'ping' object with the following arguments:

\describe{
\item{w}{The mixture weights (a vector)}
\item{mu}{The binding site positions (a vector)}
\item{delta}{The DNA fragment lengths (a vector)}
\item{sigmaSqF}{The variance parameters for the forward distribution (vector)}
\item{sigmaSqR}{The variance parameters for the forward distribution (vector)}
\item{seMu}{The standard errors for mu (vector)}
\item{seMuF}{The standard errors for muF (vector)}
\item{seMuR}{The standard errors for muR (vector)}
\item{score}{The scores for each binding event (vector)}
\item{Nmerged}{The number of peaks that got merged (integer)}
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
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" GenomeBiology, under review.
}

\seealso{
\code{\link{ping}}
\code{\link{pingError}}
}


\examples{
# Here is an example of how to construct such a region.
# Typically, you would not do this manually, you would use the ping function to return a 'pingList'
# that contains a list of 'ping' or a 'pingError' object.
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
chr<-"chr1"
range<-c(1000,2000)
# Contructor
#myPING<-newPing(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,Nmerged,as.integer(range),chr)
}


\keyword{models}
