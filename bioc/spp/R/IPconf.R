setClassUnion("smoothedEnrich_Or_NULL", c("smoothedEnrich", "NULL"))
setClassUnion("conservEnrich_Or_NULL", c("conservEnrich", "NULL"))
setClassUnion("broadRegion_Or_NULL", c("broadRegion", "NULL"))
setClassUnion("bindingPos_Or_NULL", c("bindingPos", "NULL"))
IPconf = setRefClass(
	Class = "IPconf",
	fields = list(
		.ChIP = "AlignedTags_Or_NULL", 	##AlignedTags object for ChIP
		.Input = "AlignedTags_Or_NULL", 	##AlignedTags object for Input	
		smoothed_enrichment = "smoothedEnrich_Or_NULL", 
		conserved_enrichment = "conservEnrich_Or_NULL", 
		broad_region = "broadRegion_Or_NULL", 
		binding_position = "bindingPos_Or_NULL"
	)
)
IPconf$methods(
	initialize = function(..., ChIP=NULL, Input=NULL) {
		callSuper(...)

		##check input args
		if(!is(ChIP, "AlignedTags_Or_NULL") || 
			!is(Input, "AlignedTags_Or_NULL"))
			stop("'ChIP' and 'Input' should be AlingedTags objects or 'NULL'!")
		##assign ChIP and Input
		##new object
		if(is.null(.ChIP)) {
			if(!is.null(ChIP))
				.ChIP <<- ChIP
		} 
		##otherwise, copied object
		##new Input object
		if(is.null(.Input)) {
			if(!is.null(Input))	
				.Input <<- Input
		}
		##otherwise, copied object		
		if(is.null(smoothed_enrichment))			##new
			smoothed_enrichment <<- smoothedEnrich(ChIP=ChIP, Input=Input)
		else 										##copy										
			smoothed_enrichment <<- smoothedEnrich(.self$smoothed_enrichment)
		if(is.null(conserved_enrichment))			##new
			conserved_enrichment <<- conservEnrich(ChIP=ChIP, Input=Input)
		else 										##copy
			conserved_enrichment <<- conservEnrich(.self$smoothed_enrichment)
		if(is.null(broad_region))					##new
			broad_region <<- broadRegion(ChIP=ChIP, Input=Input)
		else 										##copy
			broad_region <<- broadRegion(.self$broad_region)			
		if(is.null(binding_position))				##new
			binding_position <<- bindingPos(ChIP=ChIP, Input=Input)
		else 										##copy
			binding_position <<- bindingPos(.self$binding_position)	
	}
)
IPconf$methods(
	show = function(...) {	

		##message about Input
		if(!is.null(.Input)) {
			cat("Input:\n")
			if(!is.null(.Input$tags)) {
				nchrs <- length(.Input$tags)
				ntags <- sum(unlist(lapply(.Input$tags, length)))
				cat(paste("  ", ntags, " fragments", " across ", nchrs, 
					" chromosome(s)", "\n", sep=""))			
			}		
		}
		##message about ChIP
		if(!is.null(.ChIP)) {
			cat("ChIP:\n")
			if(!is.null(.ChIP$tags)) {
				nchrs <- length(.ChIP$tags)
				ntags <- sum(unlist(lapply(.ChIP$tags, length)))
				cat(paste("  ", ntags, " fragments", " across ", nchrs, 
					" chromosome(s)", "\n", sep=""))			
			}
			##show binding characteristics
			if(!is.null(.ChIP$bd_chrtcs)) {
				##
				cat("Binding characteristics:\n")
				cat(paste("  cross correlation peak: Position=", 
					.ChIP$bd_chrtcs$peak$x, ", Height=", 
					round(.ChIP$bd_chrtcs$peak$y, 3), "\n", sep=""))
				cat(paste("  optimized window half-size: ", 
					.ChIP$bd_chrtcs$whs, "\n", sep=""))
			}							
		}
	}
)
