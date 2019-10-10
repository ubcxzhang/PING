\name{show}
\Rdversion{1.1}
\docType{methods}
\alias{show-methods}
\alias{show}
\alias{show,ping-method}
\alias{show,pingError-method}
\alias{show,pingList-method}
\alias{show,segReads-method}
\alias{show,segReadsList-method}
\title{show PING } 

\description{
This methods show the objects of PING
}

\usage{
\S4method{show}{ping}(object)
\S4method{show}{pingError}(object)
\S4method{show}{pingList}(object)
\S4method{show}{segReads}(object)
\S4method{show}{segReadsList}(object)
}

\arguments{
\item{object}{Object returned from \code{\link{ping}} .}
}

\details{
List of the slots include in the object
}

\author{
Xuekui Zhang <\email{xzhang@stat.ubc.ca}
Sangsoon Woo, \email{swoo@fhcrc.org}
Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}
}


\seealso{
\code{\link{summary}}
}

\keyword{print}
