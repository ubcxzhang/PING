\name{NucleosomeTrack}
\alias{NucleosomeTrack}

\title{Track for the Nucleosome position}

\description{
This track actually creates two tracks: it displays the position of the predicted nucleosomes as well as the standard error associated with each prediction and an histogram of the score.
}

\usage{
NucleosomeTrack(PS, chr, gen="gen", scoreTrack=TRUE, scoreThreshold=0.05, name="PING", ...)
}

\arguments{
\item{PS}{A \code{data.frame} or \code{\linkS4class{pingList}}. The output of \code{postPING} or \code{PING}.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}
\item{scoreTrack}{A \code{logical}. If set to FALSE, the track showing the histogram of the score is removed.}
\item{scoreThreshold}{A \code{numeric}. Nucleosomes predicted with a score inferior to this threshold are removed.}
\item{name}{A character string for the title of the track.}
\item{...}{A list of supplemental argument that can be provided to the \code{\link[Gviz]{AnnotationTrack}} track constructor.}
}

\value{
This function returns a list containing an object of class \code{AnnotationTrack} and a \code{DataTrack} if scoreTrack is set to TRUE.
}


\author{Renan Sauteraud}

\seealso{
\code{\link{postPING}}
\code{\link{plotSummary}}
\code{\link[Gviz]{DataTrack}}
\code{\link[Gviz]{AnnotationTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}
