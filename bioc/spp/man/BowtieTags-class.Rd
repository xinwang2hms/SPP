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
	Methods inherited from class \code{\link[spp:AlignedTags]{AlignedTags}}, 
	so please look into class \code{\link[spp:AlignedTags]{AlignedTags}} 
	for more details:
	\describe{
		\item{\code{subset(filter)}}{
		}
		\item{\code{remove.tag.anomalies(bin=1, trim_fraction=1e-3, z=5, 
		zo=3*z, var_base=0.1)}}{
		}
		\item{\code{remove.local.tag.anomalies(window_size=200, eliminate_fold=10, 
		cap_fold=4, z_threshold=3)}}{
		}		
		\item{\code{compute.cross.cor(srange=c(50, 500), bin=5, min_tag_count=1e3,
		acceptance_z_score=3, accept_all_tags=FALSE)}}{
		}
		\item{\code{get.cross.cor()}}{
		}	
		\item{\code{set.cross.cor(value)}}{
		}	
		\item{\code{view.cross.cor()}}{
		}	
	}
	Methods specific to class \code{BowtieTags}:
	\describe{
		\item{\code{read(file=NULL, read_tag_names=F, read_tag_qualities=F, 
		fix_chr_names=F)}}{
		Read in aligned reads from Bowtie file.
		\emph{Arguments} \cr
		\code{file}: The file of bowtie alignment results. \cr
		\code{read_tag_names}: whether the tag names should be read in. \cr
		\code{read_tag_qualities}: not in use. \cr
		\code{fix_chr_names}: Whether to remove ".fa" from the end of the sequence names. \cr
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













