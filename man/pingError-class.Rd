\name{pingError-class}
\docType{class}
\alias{pingError}
\alias{pingError-class}
\alias{chromosome,pingError-method}
\alias{mu,pingError-method}
\alias{delta,pingError-method}
\alias{w,pingError-method}
\alias{score,pingError-method}
\alias{se,pingError-method}
\alias{seF,pingError-method}
\alias{seR,pingError-method}
\alias{newPingError}
\alias{newPingError,pingError-method}

\title{The ping class}

\description{
This object is used to return an error code when the PING function failed to return a valid set of estimates for a candidate regions. This could be due to non-convergence of the EM algorithm, a singular information matrix, or a number of reads below the limit specified by the user. All of these are typically due to too few reads in the region and do not affect the rest of the analysis, as such regions would most likely be labelled as false positives. 
}

\section{Accessors}{
All of the accessors defined for a 'ping' object still work for a 'pingError' object but will simply return a NULL pointer.
}

\section{Constructor}{
\describe{
newPingError(string) where 'string' is the error code.
}
}

\section{Constructor}{
newPingError<-function(string)
\describe{
\item{string}{The mixture weights (a vector)}
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
}


\examples{
# Here is an example on how to construct such a pingError object
# Typically, you would not do this manually, you would use the ping function to return a 'pingList'
# that contains a list of 'ping' or 'pingError' object.
# Contructor
myPingError<-newPingError("Singular information matrix")
# Accessors
# Get the standard error of Mu
se(myPingError)
# Get the standard error of MuF
seF(myPingError)
# Get the scores
score(myPingError)
}


\keyword{models}
