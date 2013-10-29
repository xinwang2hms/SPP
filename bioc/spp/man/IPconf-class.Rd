\name{IPconf-class}
\docType{class}
\alias{IPconf-class}
\alias{IPconf}
\title{
A reference class representing IP configuration 
}
\description{
The reference class is designed to represent an IP configuration for a 
ChIP-Seq experiment.
}
\section{Objects from the class}{
Objects of class \code{IPconf} can be created from \code{IPconf(ChIP,
Input)}. 
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
		\item{\code{smoothed_enrichment}:}{
		an object of class \code{\link[spp:smoothedEnrich]{smoothedEnrich}}, which represents 
		smoothed enrichment profile
		}
		\item{\code{conserved_enrichment}:}{
		an object of class \code{\link[spp:conservEnrich]{conservEnrich}}, which represents 
		conserved enrichment profile
		}
		\item{\code{broad_region}:}{
		an object of class \code{\link[spp:broadRegion]{broadRegion}}, which represents broad 
		region profile
		}
		\item{\code{binding_position}:}{
		an object of class \code{\link[spp:bindingPos]{bindingPos}}, which represents binding 
		position profile
		}
	}
}
\section{Methods}{
	The class contains member variables that are objects of specific ChIP-seq 
profiles. Check classes representing specific profiles for available methods.
}
\seealso{
	\code{\link[spp:ChIPSeqProfile]{ChIPSeqProfile}}
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













