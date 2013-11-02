\name{ChIPSeqProfile-class}
\docType{class}
\alias{ChIPSeqProfile-class}
\alias{ChIPSeqProfile}
\alias{set.ChIP,ChIPSeqProfile-method}
\alias{set.Input,ChIPSeqProfile-method}
\alias{set.param,ChIPSeqProfile-method}
\alias{get.profile,ChIPSeqProfile-method}
\title{
A reference class representing ChIP-Seq profiles
}
\description{
The reference class represents various kinds of ChIP-Seq profiles
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
		\item{\code{get.profile}}{
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













