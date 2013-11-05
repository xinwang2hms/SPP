\name{ChIPSeqProfile-class}
\docType{class}
\alias{ChIPSeqProfile-class}
\alias{ChIPSeqProfile}
\alias{set.ChIP,ChIPSeqProfile-method}
\alias{set.Input,ChIPSeqProfile-method}
\alias{set.param,ChIPSeqProfile-method}
\alias{get.param,ChIPSeqProfile-method}
\alias{get.profile,ChIPSeqProfile-method}
\alias{write.wig,ChIPSeqProfile-method}
\alias{write.tdf,ChIPSeqProfile-method}
\alias{view,ChIPSeqProfile-method}
\title{
A reference class representing ChIP-Seq profiles
}
\description{
The reference class represents various sorts of ChIP-Seq profiles
}
\section{Objects from the class}{
Objects of class \code{ChIPSeqProfile} should be created from its 
subclasses such as \code{\link[spp:smoothedEnrich]{smoothedEnrich}}.
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
		a list of parameters used to compute profile
		}
	}
}
\section{Methods}{
	\describe{
		\item{\code{set.ChIP(ChIP)}}{
		Assign a reference to a ChIP sample. \cr
		\emph{Arguments} \cr
		\code{ChIP}: an object of \code{AlignedTags} representing a ChIP sample. \cr
		}		
		\item{\code{set.Input(Input)}}{
		Assign a reference to an Input sample. \cr
		\emph{Arguments} \cr
		\code{Input}: an object of \code{AlignedTags} representing a ChIP sample. \cr
		}			
		\item{\code{set.param(...)}}{
		Set parameters for computing a ChIP-seq profile. \cr
		\emph{Arguments} \cr
		\code{...}: parameters allowed by field \code{.param} (see details in 
		the description of field \code{.param} in subclasses inheriting this class). \cr
		}			
		\item{\code{get.param(what)}}{
		Get parameter(s). \cr
		\emph{Arguments} \cr
		\code{what}: a character specifying which parameter to get. If \code{what="all"}, 
		it returns all available parameters as a list. \cr
		}			
		\item{\code{write.wig(file, name, feature, threshold=5, zip=F)}}{
		Write profile to a wig file. \cr
		\emph{Arguments} \cr
		\code{file}: name of the file to save. \cr
		\code{name}: label of the wig file. \cr
		\code{feature}: description of the wig file. \cr
		\code{zip}: whether to compress the wig file or not. \cr
		}			
		\item{\code{write.tdf(file, name, feature, save_wig=F, zip_wig=T)}}{
		Write profile in tdf format which can be loaded to IGV browser. \cr
		\emph{Arguments} \cr
		\code{file}: name of the file to save. \cr
		\code{name}: label of the tdf file. \cr
		\code{feature}: description of the tdf file. \cr
		\code{save_wig}: whether or not to save the wig file as well. \cr
		\code{zip_wig}: If \code{save_wig=TRUE}, whether or not to zip the wig file. \cr
		}			
		\item{\code{view(chr=NULL, start=0, end=Inf, col_sig="red", col_bg="green", ...)}}{
		Preview the profile in a specified genomic region. If start and end are both ignored, 
		the method will display profile in the whole chromosome. \cr
		\emph{Arguments} \cr
		\code{chr}: chromosome name, e.g. "chr1". \cr
		\code{start}: start of a genomic region. If \code{start=0}, the 5' end where tags are 
		available will be used. \cr
		\code{end}: end of a genomic region. If \code{end=Inf}, the 3' end where tags are 
		available will be used. \cr
		\code{col_sig}: color used to highlight enrichment when signal is greater 
		than control. \cr
		\code{col_bg}: color used to highlight enrichment when signal is lower than 
		control. \cr
		\code{...}: additional arguments accepted by \code{plot}. 
		}
		\item{\code{get.profile()}}{
		Return the profile computed. 
		}
	}
}
\seealso{
	\code{\link[spp:smoothedTagDensity]{smoothedTagDensity}}
	\code{\link[spp:smoothedEnrich]{smoothedEnrich}}
	\code{\link[spp:conservEnrich]{conservEnrich}}
	\code{\link[spp:bindingPos]{bindingPos}}
	\code{\link[spp:broadRegion]{broadRegion}}
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













