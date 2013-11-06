\name{broadRegion-class}
\docType{class}
\alias{broadRegion-class}
\alias{broadRegion}
\alias{identify,broadRegion-method}
\alias{view,broadRegion-method}
\alias{write.broadpeak}
\title{
A reference class representing broad region profiles
}
\description{
The reference class represents broad region profiles
}
\section{Objects from the class}{
Objects of class \code{broadRegion} can be created from 
\code{broadRegion(ChIP, Input, param)}.
}
\section{Fields}{
	\describe{
		\item{\code{.ChIP}:}{
		a reference to an object of \code{\link[spp:AlignedTags]{AlignedTags}} class, which 
		represents aligned tags of an ChIP experiment
		}
		\item{\code{.Input}:}{
		a reference to an object of \code{\link[spp:AlignedTags]{AlignedTags}} class, which 
		represents aligned tags of an Input experiment
		}
		\item{\code{.profile}:}{
		a variable for caching profile temporarily
		}
		\item{\code{.param}:}{
		a list of parameters used to compute profile. \cr

		\code{window_size=500}: window size to be used for tag counting. \cr

		\code{tag_shift=NULL}: Distance by which the positive and negative strand tags should 
		be shifted towards eachother. This normally corresponds to the half of the 
		cross-correlation peak position. The parameter should be set at initialization, if 
		cross correlation has already been calculated in the object of class \code{AlignedTags}. \cr

		\code{bg_weight=NULL}:optional weight by which the background density should be 
		multipled for scaling. If not supplied, the weight is calculated based on the ratio of 
		the reduced ChIP to input dataset sizes. The parameter is calculated and set at the 
		first time the profile is computed. \cr

		\code{bg_density_scaling=TRUE}: If TRUE, regions of significant tag enrichment will 
		be masked out when calculating size ratio of the signal to control datasets (to estimate 
		ratio of the background tag density). If FALSE, the dataset ratio will be equal to the 
		ratio of the number of tags in each dataset. \cr

		\code{chrl=NULL}: a vector of chromosomes. \cr

		\code{rngl=NULL}: a list of chromosomal regions. \cr

		}
		\item{\code{.param.updated}:}{
		a logical value indicating whether the parameter has been updated
		}
	}
}
\section{Methods}{
	Methods inherited from \code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}}, so please look into 
	class \code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}} for more details:
	\describe{
		\item{\code{set.ChIP(ChIP)}}{
		}		
		\item{\code{set.Input(Input)}}{
		}			
		\item{\code{get.param(what)}}{
		}			
	}
	Methods overridden from \code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}}:
	\describe{
		\item{\code{set.param(...)}}{
		}
		\item{\code{get.profile(sort_by="chr")}}{
		get the profile in a data frame. \cr
		\emph{Arguments} \cr
		\code{sort_by}: how to sort the data frame, by chromosome (\dQuote{chr}) or 
		by fold enrichment score (\dQuote{score}). \cr
		}
	}
	Methods specific to \code{broadRegion}:
	\describe{
		\item{\code{identify()}}{
		identify broad regions of enrichment. \cr
		}
		\item{\code{view(chr=NULL, start=NULL, end=NULL, col_sig="red",
		col_bg="green", ...)}}{
		preview identified broad regions in specified genomic region. If \code{start} and 
		\code{end} are both ignored, the method will display profile in the whole chromosome. \cr
		Arguments are the same as the root class \code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}}. \cr
		}
		\item{\code{write.broadpeak(file)}}{
		save identified broad regions in the broadPeak format. \cr
		\emph{Arguments} \cr
		\code{file}: the name of the file to save the identified broad regions. \cr
		}
	}
}

\seealso{
	\code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}}
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













