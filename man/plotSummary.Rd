\name{plotSummary}
\alias{plotSummary}

\title{
Plot a summary of the prediction for given ranges.
}

\description{
This function use Gviz package to summarize the nucleosome position prediction from postPING.
}

\usage{
plotSummary(PS, ping, reads, chr, gen="gen", from=NULL, to=NULL, FragmentLength=200, title="", scoreThreshold=0.05)
}

\arguments{
\item{PS}{A \code{data.frame}. The output \code{postPING}. If a list of \code{data.frame} is passed, multiple \code{NucleosomeTrack} will be plotted.}
\item{ping}{A \code{\linkS4class{pingList}}. The output of \code{PING}.}
\item{reads}{A \code{GRanges} object. The reads used in the segmentation process.}
\item{chr}{A character string for the chromosome name.}
\item{gen}{An optional character string for the genome name.}
\item{from}{An \code{integer}, the first base of the plot}
\item{to}{An \code{integer}, the last base of the plot}
\item{FragmentLength}{An \code{integer}, the length of XSET profile extension in the \code{CoverageTrack}. This argument will be ignored if the reads are paired-end sequencing data.}
%\item{GRT}{A \code{logical}. If set to TRUE, Gviz will attempt to create a \code{\link[Gviz]{BiomartGeneRegionTrack}} using the provided genome and chromosome arguments.}
\item{title}{An optional character string that can be prepend to the automatically generated title.}
\item{scoreThreshold}{A \code{numeric}. Removes nucleosome with a score inferior to the threshold from the plot.}
%\item{PE}{A \code{logical}. Set to true for paired-end sequencing data, where the read length is already known.}
}

\value{
This function returns an invisible \code{list} of tracks. The tracks can be used in other Gviz plots.
}

\author{
Renan Sauteraud
}

\seealso{
\code{\link{CoverageTrack}}
\code{\link{RawReadsTrack}}
\code{\link{NucleosomeTrack}}
\code{\link[Gviz]{plotTracks}}
}

\keyword{functions}
