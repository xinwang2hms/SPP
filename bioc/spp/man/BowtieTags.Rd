\name{BowtieTags-class}
\docType{class}
\alias{BowtieTags-class}
\alias{BowtieTags}
\alias{read,BowtieTags-method}
\title{
A reference class representing aligned tags in the Bowtie format 
}
\description{
The reference class is a subclass of class \code{\link[spp:AlignedTags]{AlignedTags}}
}
\section{Objects from the class}{
Objects of class \code{BowtieTags} can be created in multiple ways (see examples 
below for details).
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
	Methods inherited from class \code{\link[spp:AlignedTags]{AlignedTags}}:
	\describe{
		\item{\code{subset}}{
		}
		\item{\code{remove.tag.anomalies}}{
		}
		\item{\code{remove.local.tag.anomalies}}{
		}		
		\item{\code{compute.bd.chrtcs}}{
		}
		\item{\code{get.bd.chrtcs}}{
		}	
		\item{\code{set.bd.chrtcs}}{
		}	
		\item{\code{get.smoothed.density}}{
		}	
		\item{\code{view.cross.cor}}{
		}	
	}
	Methods specific to class \code{BowtieTags}:
	\describe{
		\item{\code{read}}{
		}
	}
}
\seealso{
	\code{\link[spp:smoothedTagDensity]{smoothedTagDensity}}
	\code{\link[spp:AlignedTags]{AlignedTags}}
	\code{\link[spp:BAMTags]{BAMTags}}
}
\author{Xin Wang \email{xin_wang@hms.harvard.edu}}
\references{
Kharchenko PV, Tolstorukov MY and Park PJ, Design and analysis of 
ChIP-seq experiments for DNA-binding proteins. Nat Biotechnol 2008, 
26(12):1351-9
}

\examples{
\dontrun{
encH3K4 <- BowtieTags(genome_build="mm9")
encH3K4$read(file=file.path("H3K4me3.btm.bz2"), read_tag_names=F, 
read_tag_qualities=F, fix_chr_names=F)
}
}

\keyword{classes}













