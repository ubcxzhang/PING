\name{RawReadsTrack}
\alias{RawReadsTrack}

\title{Raw reads density track}

\description{
This track displays the start position of the forward and reverse reads used as input for the segmentation step}

\usage{
RawReadsTrack(ping, reads, chr, gen="gen", from=NULL, to=NULL, ...)
}

\arguments{
\item{ping}{A \code{\linkS4class{pingList}}. The output of \code{PING}.}
\item{reads}{An object of class \code{\linkS4class{GRanges}}. The reads used in the segmentation process.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}
\item{from, to}{Two \code{numeric}. A subset of the reads can be selected in order to speed up the computation.}
%\item{PE}{A \code{logical}. If set to TRUE, the start of the reads is used as forward reads, while the end is used as reverse reads.}
\item{...}{Optional arguments (\code{DisplayParameters}) that can be passed on to the \code{DataTrack}. Refer to Gviz help and vignette for a list of arguments.}

}

\value{
This function returns an object of class \code{DataTrack}.
}

\author{Renan Sauteraud}

\seealso{
\code{\link{segmentPING}}
\code{\link{plotSummary}}
\code{\link[Gviz]{DataTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}
