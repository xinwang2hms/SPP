\name{broadRegion-class}
\docType{class}
\alias{broadRegion-class}
\alias{broadRegion}
\alias{identify,broadRegion-method}
\alias{get,broadRegion-method}
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
	}
}
\section{Methods}{
	Methods inherited from \code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}}:
	\describe{
		\item{\code{set.ChIP}}{
		}		
		\item{\code{set.Input}}{
		}			
		\item{\code{set.param}}{
		}			
		\item{\code{get.param}}{
		}			
		\item{\code{write.wig}}{
		}			
		\item{\code{write.tdf}}{
		}			
		\item{\code{view}}{
		}
	}
	Methods specific to \code{broadRegion}:
	\describe{
		\item{\code{identify}}{
		}
		\item{\code{write.broadpeak}}{
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













