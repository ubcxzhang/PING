\name{CoverageTrack}
\alias{CoverageTrack}

\title{Reads coverage track}

\description{
This track displays the coverage of the genomic region by the reads used as input for the segmentation step.}

\usage{
CoverageTrack(ping, reads, chr, gen="gen", FragmentLength=200, name="XSET")
}

\arguments{
\item{ping}{An object of class \code{\linkS4class{pingList}}. The output of \code{PING}.}
\item{reads}{An object of class \code{\linkS4class{GRanges}}. The reads used in the segmentation process.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}
\item{FragmentLength}{An \code{integer}. The fixed length of the reads. If PE is set to TRUE, this argument can be ignored.}
\item{name}{A \code{character} string. The title of the track.}
%\item{PE}{A \code{logical}. Set to true for paired-end sequencing data, where the read length is already known.}
}

\value{
This function returns an object of class \code{AnnotationTrack}.
}

\author{Renan Sauteraud}

\seealso{
\code{\link{postPING}}
\code{\link{plotSummary}}
\code{\link[Gviz]{AnnotationTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}
