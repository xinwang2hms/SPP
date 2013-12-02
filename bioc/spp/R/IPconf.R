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
		conserv_enrichment = "conservEnrich_Or_NULL", 
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
##			smoothed_enrichment <<- smoothedEnrich(.self$smoothed_enrichment)
			smoothed_enrichment <<- smoothedEnrich(smoothed_enrichment)
		if(is.null(conserv_enrichment))			##new
			conserv_enrichment <<- conservEnrich(ChIP=ChIP, Input=Input)
		else 										##copy
##			conserv_enrichment <<- conservEnrich(.self$smoothed_enrichment)
			conserv_enrichment <<- conservEnrich(smoothed_enrichment)
		if(is.null(broad_region))					##new
			broad_region <<- broadRegion(ChIP=ChIP, Input=Input)
		else 										##copy
##			broad_region <<- broadRegion(.self$broad_region)			
			broad_region <<- broadRegion(.self$broad_region)			
		if(is.null(binding_position))				##new
			binding_position <<- bindingPos(ChIP=ChIP, Input=Input)
		else 										##copy
##			binding_position <<- bindingPos(.self$binding_position)	
			binding_position <<- bindingPos(binding_position)	
	}
)

##save object of IPconf to rds file
IPconf$methods(
	save = function(file, what="config", ...) {
		if(missing(file))
			stop("Please specify file name!")
		if(! what %in% c("config", "all"))
			stop("Please specify what to save: only configuration or including AlignedTags objects?")
		if(what=="config") {
			.ChIP <<- NULL
			.Input <<- NULL
			smoothed_enrichment$.ChIP <<- NULL
			smoothed_enrichment$.Input <<- NULL
			conserv_enrichment$.ChIP <<- NULL
			conserv_enrichment$.Input <<- NULL
			binding_position$.ChIP <<- NULL
			binding_position$.Input <<- NULL
			broad_region$.ChIP <<- NULL
			broad_region$.Input <<- NULL
		}
		conf <- .self
		saveRDS(conf, file=file)
}
)


##set ChIP and Input
##!should set recursively for each attribute
##!if more attributes are added in the future, this part 
##should be revised
IPconf$methods(
	set.ChIP = function(ChIP) {
		.ChIP <<- ChIP
		smoothed_enrichment$.ChIP <<- ChIP
		conserv_enrichment$.ChIP <<- ChIP
		binding_position$.ChIP <<- ChIP
		broad_region$.ChIP <<- ChIP
	}
)
IPconf$methods(
	set.Input = function(Input) {
		.Input <<- Input
		smoothed_enrichment$.Input <<- Input
		conserv_enrichment$.Input <<- Input
		binding_position$.Input <<- Input
		broad_region$.Input <<- Input
	}
)




IPconf$methods(
	show = function(...) {	
		##General message
		cat("~~IP configuration~~\n")
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
