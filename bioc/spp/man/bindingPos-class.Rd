\name{bindingPos-class}
\docType{class}
\alias{bindingPos-class}
\alias{bindingPos}
\alias{identify,bindingPos-method}
\alias{get.bp}
\alias{write.bp}
\alias{view.bp}
\alias{get.narrowpeak}
\alias{write.narrowpeak}
\title{
A reference class representing point binding position profiles
}
\description{
The reference class is represents point binding position profiles
}
\section{Objects from the class}{
Objects of class \code{bindingPos} can be created from 
\code{bindingPos(ChIP, Input, param)}. 
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

		\emph{position stringency criteria} \cr
		
		\code{e_value=NULL}: E-value defining the desired statistical significance of binding positions. \cr
		
		\code{fdr=NULL}: FDR defining statistical significance of binding positions. \cr
		
		\code{topN=NULL}: instead of determining statistical significance thresholds, return the 
		specified number of highest-scoring positions. \cr
		
		\emph{other parameters} \cr
		
		\code{whs=NULL}: window half-sized that should be used for binding detection. If not specified, 
		\code{whs} will be determined from cross-correlation profiles of \code{.ChIP}. \cr
		
		\code{masked_data=NULL}: optional set of coordinates that should be masked (e.g. known 
		non-unique regions). \cr
		
		\code{min_dist=200}: minimal distance that must separate detected binding positions. In case 
		multiple binding positions are detected within such distance, the position with the 
		highest score is returned. \cr
		
		\code{window_size=4e7}: size of the window used to segment the chromosome during calculations 
		to reduce memory usage. \cr
		
		\code{min_thr=2}: minimal score requirement for a peak. \cr
		
		\code{bg_density_scaling=TRUE}: If TRUE, regions of significant tag enrichment will be 
		masked out when calculating size ratio of the signal to control datasets (to estimate ratio 
		of the background tag density). If FALSE, the dataset ratio will be equal to the ratio of 
		the number of tags in each dataset. \cr
		
		\emph{randomized controls} \cr		
		
		\code{nrand=3}: number of tag randomziations that should be performed (when the control data 
		is not provided). \cr
		
		\code{use_rand_ctr=FALSE}: use randomized tag control, even if \code{control.data} is supplied. \cr
		
		\code{shuffle_window=1}: during tag randomizations, tags will be split into groups of 
		\code{shuffle_window} and will be maintained together throughout the randomization. \cr
		
		\emph{fold-enrichment confidence intervals} \cr
		
		\code{tag_count_whs=100}: half-size of a window used to assess fold enrichment of a 
		binding position. \cr
		
		\code{enrichment_z=2}: z-score used to define the significance level of the fold-enrichment 
		confidence intervals. \cr
		
		\code{enrichment_scale_down_ctr=FALSE}: if TRUE, sample down control to be the same size as 
		true signal.data. \cr

		\code{enrichment_bg_scales=c(1,5,10)}: in estimating the peak fold-enrichment confidence 
		intervals, the background tag density is estimated based on windows with half-sizes of 
		\code{2*tag_count_whs*enrichment_bg_scales}. \cr
		
		\code{method="tag.wtd"}: either \dQuote{tag.wtd} for WTD method, or \dQuote{tag.lwcc} for MTC method. \cr
		
		\code{mle_filter=FALSE}: if TRUE, will exclude predicted positions whose MLE enrichment ratio 
		(for any of the background scales) is below a specified \code{min_mle_threshold}. \cr
		
		\code{min_mle_threshold=1}: MLE enrichment ratio threshold that each predicted position must 
		exceed if \code{mle_filter} is TRUE.
		
		\emph{masking regions of significant control enrichment} \cr
		
		\code{tec_filter=TRUE}: Whether to mask out the regions exhibiting significant enrichment 
		in the control data in doing other calculations. The regions are identified using Poisson 
		statistics within sliding windows, either relative to the scaled signal (tec_z), or 
		relative to randomly-distributed expectation (tec_poisson_z). \cr
		
		\code{tec_window_size=1e4}: size of the window used to determine significantly enrichent 
		control regions. \cr
		
		\code{tec_masking_window_size=1e4}: size of the window used to mask the area around 
		significantly enrichent control regions. \cr
		
		\code{tec_z=5}: z-score defining statistical stringency by which a given window is determined 
		to be significantly higher in the input than in the signal, and masked if that is the case. \cr
		
		\code{tec_poisson_z=5}: z-score defining statistical stringency by which a given window is 
		determined to be significantly higher than the tec_poisson_ratio above the expected uniform 
		input background. \cr
		
		\code{tec_poisson_ratio=5}: Fold ratio by which input must exceed the level expected from 
		the uniform distribution. \cr

		\emph{broader regions of enrichment} \cr
		
		\code{add_broad_peak_reg=FALSE}: whether or not to add broader regions of enrichment associated 
		with determined peaks. \cr

		\code{broad_peak_reg_window_size=500}: window size to be used in calculating enrichment. \cr

		\code{broad_peak_reg_z_thr=2}: z-score corresponding to the Poisson ratio threshold used 
		to flag significantly enriched windows. \cr
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
		\item{\code{get.profile()}}{
		}
	}
	Methods specific to \code{bindingPos}:
	\describe{
		\item{\code{identify()}}{
		identify binding positions.
		}
		\item{\code{get.bp(sort_by="chr")}}{
		get binding positions in a table. If binding positions are not yet available, the method will 
		invoke \code{identify} to find binding positions first. \cr
		\emph{Arguments}
		\code{sort_by}: how should the data frame of binding positions be sorted, by 
		chromosome (\dQuote{chr}), fold enrichment score (\dQuote{score}), E-value 
		(\dQuote{Evalue}) or FDR (\dQuote{FDR}). 
		}
		\item{\code{write.bp(file, sort_by="chr")}}{
		save identified binding positions as a table. If binding positions are not yet available, 
		the method will invoke \code{identify} to find binding positions first. \cr
		\emph{Arguments}
		\code{file}: the name of the file to save. \cr
		\code{sort_by}: save as \code{get.bp}. \cr
		}
		\item{\code{view.bp(chr=NULL, start=NULL, end=NULL, col="red")}}{
		Preview binding positions in a specified genomic region. If start and end are both ignored, 
		the method will display profile in the whole chromosome. If binding positions are not yet available, 
		the method will invoke \code{identify} to find binding positions first. \cr
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
		\item{\code{get.narrowpeak(sort_by="chr")}}{
		get binding position profile in narrowPeak format. If binding positions are not yet available, 
		the method will invoke \code{identify} to find binding positions first. \cr
		\emph{Arguments} \cr
		\code{sort_by}: how should the table be sorted, by chromosome (\dQuote{chr}), 
		fold enrichment score (\dQuote{score}) or FDR (\dQuote{FDR}). \cr
		}
		\item{\code{write.narrowpeak(file)}}{
		save binding position profile in narrowPeak format. If binding positions are not yet available, 
		the method will invoke \code{identify} to find binding positions first. \cr
		\emph{Arguments} \cr
		\code{file}: the name of the file to save. \cr
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













