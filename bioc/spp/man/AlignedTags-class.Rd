\name{AlignedTags-class}
\docType{class}
\alias{AlignedTags-class}
\alias{AlignedTags}
\alias{subset}
\alias{remove.tag.anomalies}
\alias{compute.cross.cor}
\alias{get.smoothed.density}
\alias{view.cross.cor}
\alias{remove.local.tag.anomalies}
\alias{get.cross.cor}
\alias{set.cross.cor}
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
		a character value specifying the genome build such as 'mm9' and 'hg19'
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
		\item{\code{subset(filter)}}{
		Filter the object by chromosomes. \cr
		\emph{Arguments} \cr
		\code{filter}: a character vector of chromosome names: e.g. c("chr1", "chr2"). 
		}
		\item{\code{remove.tag.anomalies(bin=1, trim_fraction=1e-3, z=5, 
		zo=3*z, var_base=0.1)}}{
		Remove tag positions that have anomalously high counts on both strands. \cr
		\emph{Arguments} \cr
		\code{bin}: size of bin in basepairs. \cr
		\code{trim_fraction}: fraction of top bins to discard when calculating 
		overall background density. \cr
		\code{z}: z-score used to determine anomalous bins. \cr
		\code{zo}: z used to filter out one-strand matches. \cr
		\code{var_base}: minimal base variability of tag counts for processing 
		of flattened datasets with close to 0 variance. \cr
		}
		\item{\code{remove.local.tag.anomalies(window_size=200, eliminate_fold=10, 
		cap_fold=4, z_threshold=3)}}{
		Cap or remove tag positions that are significantly higher than local background. \cr
		\emph{Arguments} \cr
		\code{window_size}: Size of the window used to assess local density. \cr
		\code{eliminate_fold}: Threshold definining fold-over background density above 
		which the position is considered anomalous and removed completely. \cr
		\code{cap_fold}: Threshold fold-over background density above which the 
		position is capped to the maximum statistically likely given local tag density. \cr
		\code{z_threshold}: Z-score used to assess significance of a given position 
		exceeding either of the two density thresholds. \cr
		}		
		\item{\code{compute.cross.cor(srange=c(50, 500), bin=5, min_tag_count=1e3, 
		acceptance_z_score=3, accept_all_tags=FALSE)}}{
		Calculates strand cross-correlation profile to determine binding peak 
		separation distance and approximate window size that should be used for 
		binding detection. \cr
		\emph{Arguments} \cr
		\code{srange}: A range within which the binding peak separation is expected 
		to fall. Should be larger than probe size to avoid artifacts. \cr
		\code{bin}: Size of bin in basepairs specifying at which cross-corrrelation 
		should be calculated. \cr
		\code{min_tag_count}: minimal number of tags on the chromosome to be 
		considered in the cross-correlation calculations. \cr
		\code{acceptance_z_score}: A Z-score used to determine if a given tag quality 
		bin provides significant improvement to the strand cross-correlation. \cr
		\code{accept_all_tags}: Whether tag alignment quality calculations should be 
		skipped and all available tags should be accepted in the downstream analysis. \cr
		}
		\item{\code{get.cross.cor()}}{
		Get binding characteristics, which is a list including \sQuote{cross.correlation} 
		(cross correlation profile as a $x/$ydata frame), \sQuote{peak} (Position $x 
		and height $y of detected cross-correlation peak) and \sQuote{whs} (Optimized window 
		half-size for binding detection). \cr
		}	
		\item{\code{set.corss.cor(value)}}{
		Set binding characteristics. \code{value} should be binding characteristics, 
		which is a list including \sQuote{cross.correlation}, \sQuote{peak} and 
		\sQuote{whs}. See \code{get.cross.cor} for more details. \cr
		}	
		\item{\code{view.cross.cor()}}{
		Plot cross correlation profile and highlight detected peak. \cr
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
