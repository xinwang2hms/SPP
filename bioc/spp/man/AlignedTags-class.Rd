\name{AlignedTags-class}
\docType{class}
\alias{AlignedTags-class}
\alias{AlignedTags}
\alias{subset,AlignedTags-methods}
\alias{remove.tag.anomalies,AlignedTags-methods}
\alias{compute.cross.cor,AlignedTags-methods}
\alias{get.smoothed.density,AlignedTags-methods}
\alias{view.cross.cor,AlignedTags-methods}
\alias{remove.local.tag.anomalies,AlignedTags-methods}
\alias{get.cross.cor,AlignedTags-methods}
\alias{set.cross.cor,AlignedTags-methods}
\title{
A reference class representing aligned tags
}
\description{
This reference class represents aligned tags of a ChIP-Seq experiment
}
\section{Objects from the class}{
Objects of class \code{AlignedTags} should be created from its 
subclasses such as \code{BAMTags} and \code{BowtieTags}.
}
\section{Fields}{
	\describe{
		\item{\code{genome_build}:}{
		a character value specifying the genome build such as `mm9' and `hg19'
		}
		\item{\code{tags}:}{
		a list of vectors of 5' tag coordinates, with negative values
        	corresponding to tags mapped to the negative strand
		}
		\item{\code{quality}:}{
		a list of vectors of number of mismatches
		}
		\item{\code{names}:}{
		a list of vectors of tag names
		}
		\item{\code{file}:}{
		a character value specifying the path to the file, from which 
		the aligned tags are read
		}
		\item{\code{bd_chrtcs}:}{
		a list including cross correlation, peak position and height and 
		optimized window half-size for binding detection
		}
		\item{\code{smoothed_density}:}{
		an object of class \code{smoothedTagDensity}, which represents the 
		smoothed tag density profile
		}
	}
}
\section{Methods}{
	\describe{
		\item{\code{subset}}{
		}
		\item{\code{remove.tag.anomalies}}{
		}
		\item{\code{remove.local.tag.anomalies}}{
		}		
		\item{\code{compute.cross.cor}}{
		}
		\item{\code{get.cross.cor}}{
		}	
		\item{\code{set.corss.cor}}{
		}	
		\item{\code{get.smoothed.density}}{
		}	
		\item{\code{view.cross.cor}}{
		}	
	}
}
\seealso{
	\code{\link[spp:smoothedTagDensity]{smoothedTagDensity}}
	\code{\link[spp:BAMTags]{BAMTags}}
	\code{\link[spp:BowtieTags]{BowtieTags}}
}
\author{Xin Wang \email{xin_wang@hms.harvard.edu}}
\references{
Kharchenko PV, Tolstorukov MY and Park PJ, Design and analysis of 
ChIP-seq experiments for DNA-binding proteins. Nat Biotechnol 2008, 
26(12):1351-9
}

\examples{
\dontrun{

}
}

\keyword{classes}
