\name{segmentPING} 
\alias{segmentPING}


\title{Segment the genome into candidate regions} 
\description{ 
    Pre-process bidirectional aligned reads data to detect candidate regions with a minimum number of forward and reverse reads. These candidate regions will then be processed by PING.}
\usage{
	segmentPING(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3,
    jitter=FALSE, maxLregion=1200, minLregion=80, step=NULL, width=NULL,
    islandDepth=3, min_cut=50, max_cut=1000, maxReadsWidth=500, PE=FALSE)
}


\arguments{
  \item{data}{A \code{\linkS4class{GRanges}} object containing the IP reads. See details for more information on how to set up the data.}
  \item{dataC}{A \code{\linkS4class{GRanges}} object containing the control reads. Set to NULL by default, i.e. no control.}
  \item{map}{A `RangedData' object containing the mappability profiles. Set to NULL by default, i.e. no profiles.}
  \item{minReads}{The minimum number of F/R reads to be present in the sliding window.}
  \item{minReadsInRegion}{The minimum number of F/R reads to be present in the region.}
  \item{jitter}{A logical value stating whether some noise should be added to the read locations. This is recommended if the read positions have lots of duplicates.}
  \item{step}{An \code{integer}. The step size for the sliding window.}
  \item{width}{An \code{integer}. The size of the sliding window.}
  %\item{dataType}{Type of experiment.}
  \item{maxLregion}{The maximum length.}
  \item{minLregion}{The minimum length.}
  \item{PE}{A \code{logical}. Set to TRUE for Paired-End sequencing data.}
  \item{islandDepth, min_cut, max_cut, maxReadsWidth}{Parameters used for paired-end sequencing data segmentation.}
	\describe{
		\item{islandDepth}{An \code{integer}. The minimum number of reads to cover a candidate region.}
		\item{min_cut}{An \code{integer}. The minimum length of a candidate region.}
		\item{max_cut}{An \code{integer}. The maximum length of a candidate region.}
		\item{maxReadsWidth}{An \code{integer}. Reads with width superior to this limit will be removed from the data.}
		%\item{chr}{A character string for the name of the chromosome.}
	}

}

\value{
An object of class \code{\linkS4class{segReadsList}} containing the results for all regions pre-processed. 
}
\references{
Xuekui Zhang, Gordon Robertson, Sangsoon Woo, Brad G. Hoffman, and Raphael Gottardo, "Probabilistic Inference for Nucleosome Positioning with MNase-based or Sonicated Short-read Data" PlosONE.
}

\author{ Xuekui Zhang, Arnaud Droit <\email{arnaud.droit@crchuq.ualaval.ca}> and Raphael Gottardo <\email{rgottard@fhcrc.org}>}


\examples{
# Read data
path<-system.file("extdata",package="PING")
dataIP<-read.table(file.path(path,"GSM351492_R4_chr1.bed"),header=TRUE)
dataIP<-as(dataIP,"GRanges")
seg<-segmentPING(dataIP,  minReads=NULL, maxLregion=1200,minLregion=80, jitter=TRUE)
}

\seealso{
  \code{\linkS4class{segReadsList}}
}

\keyword{functions}
