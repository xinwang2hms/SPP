\name{smoothedEnrich-class}
\docType{class}
\alias{smoothedEnrich-class}
\alias{smoothedEnrich}
\title{
A reference class representing smoothed enrichment profiles
}
\description{
The reference class is designed to represent smoothed enrichment profiles
}
\section{Objects from the class}{
Objects of class \code{smoothedEnrich} can be created from 
\code{smoothedEnrich(ChIP, Input, param)}.
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
		a list of parameters used to compute profile. Parameters (and default settings) are: \cr
		\code{tag_shift=NULL}: Distance by which the positive and negative strand tags should 
		be shifted towards eachother. This normally corresponds to the half of the 
		cross-correlation peak position. The parameter should be set at initialization, if 
		cross correlation has already been calculated in the object of class \code{AlignedTags}. \cr
		\code{bg_weight=NULL}:optional weight by which the background density should be 
		multipled for scaling. If not supplied, the weight is calculated based on the ratio of 
		the reduced ChIP to input dataset sizes. The parameter is calculated and set at the 
		first time the profile is computed. \cr
		\code{bandwidth=150}: standard deviation of the Gaussian kernel. \cr
		\code{step=50}: The distance between the regularly spaced points for which the values 
		should be calculated. The parameter can be tuned to achieve a balance between the 
		resolution and memory consumption when visualizing the profile using method \code{view}. \cr
		\code{bg_density_scaling=TRUE}: If TRUE, regions of significant tag enrichment will 
		be masked out when calculating size ratio of the signal to control datasets (to estimate 
		ratio of the background tag density). If FALSE, the dataset ratio will be equal to the 
		ratio of the number of tags in each dataset. \cr
		\code{chrl=NULL}: a vector of chromosomes. \cr
		\code{rngl=NULL}: a list of chromosomal regions. \cr
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
		\item{\code{write.wig(file, name, feature, threshold=5, zip=F)}}{
		}			
		\item{\code{write.tdf(file, name, feature, save_wig=F, zip_wig=T)}}{
		}			
		\item{\code{view(chr=NULL, start=NULL, end=NULL, col_sig="red", col_bg="green")}}{
		}
	}
	Methods overridden from \code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}}:
	\describe{
		\item{\code{set.param(...)}}{
		}			
		\item{\code{get.profile()}}{
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













