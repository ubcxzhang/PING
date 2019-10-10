\name{ping}
\alias{ping}
\alias{PING}
\alias{ping-class}
\alias{code}
\alias{code,ping-method}
\alias{code,pingError-method}
\alias{code,pingList-method}
\alias{coerce}
\alias{coerce,pingList,RangedData-method}
\alias{coerce,pingList,data.frame-method}
\alias{coerce,data.frame,pingList-method}
\alias{coerce,data.frame,pingList-method}
\alias{density}
\alias{density,ping-method}
\alias{density,pingError-method}
\alias{score}
\alias{scoreForward}
\alias{scoreForward,ping-method}
\alias{scoreForward,pingList-method}
\alias{scoreForward,pingError-method}
\alias{scoreReverse}
\alias{scoreReverse,ping-method}
\alias{scoreReverse,pingList-method}
\alias{scoreReverse,pingError-method}
\alias{sigmaSqR}
\alias{sigmaSqR,ping-method}
\alias{sigmaSqR,pingList-method}
\alias{sigmaSqR,pingError-method}
\alias{sigmaSqF}
\alias{sigmaSqF,ping-method}
\alias{sigmaSqF,pingList-method}
\alias{sigmaSqF,pingError-method}
\alias{K}
\alias{K,ping-method}
\alias{K,pingList-method}
\alias{K,pingError-method}
\alias{plot}
\alias{plot,ping,segReads-method}
\alias{plot,pingError,segReads-method}
\alias{plot,pingList,pingList-method}
\alias{plot,data.frame,data.frame-method}
\alias{plot,pingList,segReadsList-method}
\alias{maxRange}
\alias{maxRange,ping-method}
\alias{maxRange,pingList-method}
\alias{maxRange,pingList-method}
\alias{maxRange,pingError-method}
\alias{minRange}
\alias{minRange,ping-method}
\alias{minRange,pingList-method}
\alias{minRange,pingList-method}
\alias{minRange,pingError-method}

% I think a few methods are missing, e.g. [[<- and [<-

\title{Estimation of binding site positions} 

\description{
This object contains Estimation of binding site positions and has the following  slots: 
segReadsList.
}

\usage{
PING(segReadsList, paraEM=NULL, paraPrior=NULL, dataType="MNase", detail=0, rescale=1, nCores=1)
}

\arguments{
\item{segReadsList}{This object contains segmentation of Genome}

\item{paraEM}{A list of parameters for the EM algorithm. The default parameters should be good enough for most usages.}
	\describe{
		\item{}{\code{minK}: an \code{integer}, default=0. The minimum number of binding events per region. If the value is 0, the minimum number is automatically calculated.}
		\item{}{\code{maxK}: an \code{integer}, default=0. The maximum number of binding events per region. If the value is 0, the maximum number is automatically calculated.}
		\item{}{\code{tol}: a \code{numeric}, default=1e-4. The tolerance for the EM algorithm.}
		\item{}{\code{B}: an \code{integer}, default=100. The maximum number of iterations to be used.}
		\item{}{\code{mSelect}: a character string specifying the information criteria to be used when selecting the number of binding events. Default="AIC3"}
		\item{}{\code{mergePeaks}: a \code{logical} stating whether overlapping binding events should be picked. Default=TRUE}
		\item{}{\code{mapCorrect}: a \code{logical} stating whether mappability profiles should be incorporated in the estimation, i.e: missing reads estimated. Default=TRUE}
	}
\item{paraPrior}{A list of parameters for the prior distribution. The default parameters should be good enough for most usages.}
	\describe{
		\item{}{\code{xi}: an \code{integer}, default=150. The average DNA fragment size.}
		\item{}{\code{rho}: an \code{integer}, default=1.2. A variance parameter for the average DNA fragment size distribution.}
		\item{}{\code{alpha}: an \code{integer}, default=10. First hyperparameter of the inverse Gamma distribution for sigma^2 in the PICS model}
		\item{}{\code{beta}: an \code{integer}, default=20000. Second hyperparameter of the inverse Gamma distribution for sigma^2 in the PING model}
		\item{}{\code{lambda}: an \code{integer}, default=0.000064. The precision of the prior for mu used for histone data.}
		\item{}{\code{dMu}: an \code{integer}, default=200. Our best guess for the distance between two neighboring nucleosomes.}
	}
%\item{PE}{A \code{logical}. This should be set to TRUE for paired-end sequencing data.}
\item{dataType}{A character string that can be set to use selected default parameters for the algorithm.}
\item{detail}{An \code{integer}. Additional information are printed if set to a value > 0.}
\item{rescale}{An \code{integer}.}
\item{nCores}{An \code{integer}. The number of cores that should be used in parallel by the function.}
}

\section{Methods}{
  \describe{
 \item{code}{\code{signature(x = "ping")}: return the error code for each list element (i.e. candidate region) of a PING object. If the string is empty, there were no errors.}
 \item{plot}{\code{signature(x = "ping")}: Plot all regions in the PING object. This might be long, and should only be used to plot a few regions, so subset the object before plotting.}
 \item{sigmaSqR}{\code{signature(x = "ping")}: return the variance parameter of the reverse (R) distribution for each binding event.}
 \item{sigmaSqF}{\code{signature(x = "ping")}: return the variance parameter of the forward (F) distribution for each binding event.}
 \item{score}{\code{signature(x = "ping")}: return the score for each binding event.} 
 \item{scoreF}{\code{signature(x = "ping")}: return the score  of the forward (F)  for each binding event.}
 \item{scoreR}{\code{signature(x = "ping")}: return the score  of the forward (R)  for each binding event.}
 \item{maxRange}{\code{signature(x = "ping")}: return the range maximum.}
 \item{minRange}{\code{signature(x = "ping")}: return the range minimal.}
 \item{K}{\code{signature(x = "ping")}: subset PING object.}
 \item{density}{\code{signature(x = "ping")}: return the density for each binding event.}
}
}

\author{
  Xuekui Zhang, \email{ubcxzhang@gmail.com} and
  Raphael Gottardo, \email{raphael.gottardo@ircm.qc.ca}
  Sangsoon Woo, \email{swoo@fhcrc.org}
}

\seealso{
\code{\link{ping}}
}
\keyword{models}
