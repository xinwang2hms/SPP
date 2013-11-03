\name{broadRegion-class}
\docType{class}
\alias{broadRegion-class}
\alias{broadRegion}
\alias{set.param,broadRegion-method}
\alias{identify,broadRegion-method}
\alias{view,broadRegion-method}
\alias{get.profile,broadRegion-method}
\alias{write.broadpeak,broadRegion-method}
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
		a list of parameters used to compute profile
		}
		\item{\code{.param.updated}:}{
		a logical value indicating whether the parameter has been updated
		}
	}
}
\section{Methods}{
	Methods inherited from \code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}}:
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
		}
	}
	Methods specific to \code{broadRegion}:
	\describe{
		\item{\code{identify()}}{
		}
		\item{\code{view(chr=NULL, start=NULL, end=NULL, col_sig="red",
		col_bg="green")}}{
		}
		\item{\code{write.broadpeak(file)}}{
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













