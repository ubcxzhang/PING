\name{postPING} 
\alias{postPING} 


\title{Post process Estimation of binding site positions obtained from PING} 
\description{ 
   Post process Estimation of binding site positions obtained from PING. Refit mixture models with stronger prior in candidate regions contain potential problems, and then convert final result into dataframe.}
\usage{
postPING(ping, seg, rho2=NULL, sigmaB2=NULL, alpha2=NULL, beta2=NULL, min.dist= 100, paraEM=NULL, paraPrior=NULL, score=0.05, dataType="MNase", nCores=1, makePlot=FALSE, FragmentLength=100, mart=NULL, seg.boundary=NULL, DupBound=NULL, IP=NULL, datname="")
}


\arguments{
  \item{ping}{A 'pingList' object containing estimation of nuclesome positions, result of 'PING' function.}
  \item{seg}{An object of class 'segmentReadsList' containing the results for all regions pre-processed, 'segmentReads' function.}
  \item{paraEM}{A list of parameters for the EM algorithm. The default parameters should be good enough for most usages.}
	\describe{
		\item{\code{minK}:}{ An \code{integer}, default=0. The minimum number of binding events per region. If the value is 0, the minimum number is automatically calculated.}
		\item{\code{maxK}:}{ An \code{integer}, default=0. The maximum number of binding events per region. If the value is 0, the maximum number is automatically calculated.}
		\item{\code{tol}:}{ A \code{numeric}, default=1e-4. The tolerance for the EM algorithm.}
		\item{\code{B}:}{ An \code{integer}, default=100. The maximum number of iterations to be used.}
		\item{\code{mSelect}:}{ A character string specifying the information criteria to be used when selecting the number of binding events. Default="AIC3"}
		\item{\code{mergePeaks}:}{ A \code{logical} stating whether overlapping binding events should be picked. Default=TRUE}
		\item{\code{mapCorrect}:}{ A \code{logical} stating whether mappability profiles should be incorporated in the estimation, i.e: missing reads estimated. Default=TRUE}
	}
  \item{paraPrior}{A list of parameters for the prior distribution. The default parameters should be good enough for most usages.}
	\describe{
		\item{\code{xi}:}{ An \code{integer}. The average DNA fragment size.}
		\item{\code{rho}:}{ An \code{integer}. A variance parameter for the average DNA fragment size distribution.}
		\item{\code{alpha}:}{ An \code{integer}. First hyperparameter of the inverse Gamma distribution for sigma^2 in the PICS model}
		\item{\code{beta}:}{ An \code{integer}. Second hyperparameter of the inverse Gamma distribution for sigma^2 in the PING model}
  		\item{\code{lambda}:}{ An \code{integer}. The lambda control Gaussian Markov Random Field prior on the distance of adjacent nucleosomes, we do not recommend user change the default value.}
		\item{\code{dMu}:}{ An \code{integer}. Our best guess for the distance between two neighboring nucleosomes.}
	}
  \item{rho2, sigmaB2, alpha2, beta2}{Integer values, the parameters in the prior of mixture models to be re-fitted.}
  \item{min.dist}{The minimum distance of two adjacent nucleosomes predicted from different candidate regions, smaller than that will be treated as duplicated predictions for the same nucleosomes.}
  \item{score}{A \code{numeric}. The score threshold used when calling \code{FilterPING}.}
  %\item{PE}{A \code{logical}. This should be set to TRUE for paired-end sequencing data.}
  \item{dataType}{A character string that can be set to use selected default parameters for the algorithm.}
  \item{nCores}{An \code{integer}. The number of cores that should be used in parallel by the function.}
  \item{makePlot}{A \code{logical}. Plot a summary of the output.}
  \item{FragmentLength, mart, seg.boundary, DupBound, IP, datname}{Plotting parameters and options.}
    \describe{
  		\item{\code{IP}:}{ A \code{GRanges} object. The reads used in segmentation process.}
  		\item{\code{FragmentLength}:}{ An \code{integer}. The length of XSET profile extension}
	}
}

\value{
A \code{data.frame} containing the estimation of binding site positions. 
}
\references{
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" PlosONE, under review.
}
\author{
  Xuekui Zhang <\email{xzhang@stat.ubc.ca}>, Sangsoon Woo, \email{swoo@fhcrc.org}
 and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>
}

\note{ Based on our experiemt on a few real data sets, we suggestion to use following values of parameters. For sonication data we use {rho1=1.2; sigmaB2=6400;rho=15;alpha1=10; alpha2=98; beta2=200000}.
For MNase data we use {rho1=3; sigmaB2=4900; rho=8; alpha1=20; alpha2=100; beta2=100000}. The value of xi depends on specy of sample, since that affect the length of linker-DNA. For example, we use xi=160 for yeast and xi=200 for mouse.
}

\seealso{
\code{\link{PING}}
\code{\link{plotSummary}}
}

\keyword{data}
\keyword{models}
