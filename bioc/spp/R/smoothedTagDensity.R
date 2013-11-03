##class for smoothed enrichment profile
smoothedTagDensity = setRefClass(
	Class = "smoothedTagDensity",
##	fields = list(
##
##	), 
	contains = "ChIPSeqProfile"
)
smoothedTagDensity$methods(
	initialize = function(..., ChIP=NULL, Input=NULL, param=NULL) {
		callSuper(..., ChIP=ChIP, Input=Input)
#		callSuper(...)
#if(!is.null(.ChIP) && !is(.ChIP, "AlignedTags"))
#	.ChIP <<- NULL
#if(!is.null(.Input) && !is(.Input, "AlignedTags"))
#	.Input <<- NULL
#
#		if(is.null(.ChIP)) {
#			if(!is.null(ChIP))
#				.ChIP <<- ChIP
#		}
#		if(is.null(.Input)) {
#			if(!is.null(Input))
#				.Input <<- Input
#		}
#		if(is.null(.profile))
#			.profile <<- NULL

		##3. set .param		
		if(!is(param, "list_Or_NULL"))
			stop("'param' should be a list of parameters or 'NULL'!")
		##--new object
		if(is.null(.param))	{
			##initialize default .param
			.param <<- list(tag_shift=NULL, bandwidth=150, step=50, 
				rngl=NULL, chrl=NULL, bg_density_scaling=TRUE, 
				scale_by_dataset_size=FALSE)
			##overwrite paramters according to input argument param
			if(!is.null(param)) {
				known_params <- intersect(names(param), names(.param))
				unknown_params <- setdiff(names(param), names(.param))
				if(length(unknown_params) > 0)
					warning("These parameters are unknown: ", 
						paste(unknown_params, collapse=', '))
				if(length(known_params) > 0) 
					.param[known_params] <<- param[known_params]
			}			
		if(!is.null(.ChIP)) {
			##
			if(is.null(.param$tag_shift)) {
##				cat("setting tag shift: tag_shift ")
				.param$tag_shift <<- .ChIP$bd_chrtcs$peak$x/2
##				cat("=", .param$tag_shift, " [done]\n")
			}
			##		
			##chrl		
			##determine common range
			if(is.null(.param$chrl) || (!is.list(.param$chrl) && .param$chrl=="all")) {
##				cat("setting chromosomes: chrl ")
				if(is.null(.Input))
					.param$chrl <<- names(.ChIP$tags)
				else 
					.param$chrl <<- intersect(names(.ChIP$tags), names(.Input$tags))
				names(.param$chrl) <<- .param$chrl
##				cat("[done]\n")
			}
			##rngl
			if(is.null(.param$rngl) || (!is.list(.param$rngl) && .param$rngl=="all")) {
##				cat("setting chromosome ranges: rngl ")
				if(is.null(.Input))
					.param$rngl <<- lapply(.param$chrl,function(chr) 
						range(abs(.ChIP$tags[[chr]]+.param$tag_shift)))				
				else
					.param$rngl <<- lapply(.param$chrl,function(chr) 
						range(c(range(abs(.ChIP$tags[[chr]]+.param$tag_shift)), 
						range(abs(.Input$tags[[chr]]+.param$tag_shift)))))
##				cat("[done]\n")
			} else {
				.param$chrl <<- names(.param$rngl)
				names(.param$chrl) <<- .param$chrl
			}	
		}						
		} 
		##--copied object		
	}
)
smoothedTagDensity$methods(
	set.param = function(..., verbose=TRUE) {
		callSuper(..., verbose=TRUE)
		##		
		##chrl		
		##determine common range
		if(is.null(.param$chrl) || (!is.list(.param$chrl) && .param$chrl=="all")) {
##			cat("setting chromosomes: chrl ")
			if(is.null(.Input))
				.param$chrl <<- names(.ChIP$tags)
			else 
				.param$chrl <<- intersect(names(.ChIP$tags), names(.Input$tags))
			names(.param$chrl) <<- .param$chrl
##			cat("[done]\n")
		}
		##rngl
		if(is.null(.param$rngl) || (!is.list(.param$rngl) && .param$rngl=="all")) {
##			cat("setting chromosome ranges: rngl ")
			if(is.null(.Input))
				.param$rngl <<- lapply(.param$chrl,function(chr) 
					range(abs(.ChIP$tags[[chr]]+.param$tag_shift)))				
			else
				.param$rngl <<- lapply(.param$chrl,function(chr) 
					range(c(range(abs(.ChIP$tags[[chr]]+.param$tag_shift)), 
					range(abs(.Input$tags[[chr]]+.param$tag_shift)))))
##			cat("[done]\n")
		} else {
			.param$chrl <<- names(.param$rngl)
			names(.param$chrl) <<- .param$chrl
		}	
	}
)
##4. compute smoothed tag densities

smoothedTagDensity$methods(
	get.profile = function() {


		profile <- 
			get.smoothed.tag.density(signal.tags=.ChIP$tags, 
				bandwidth=.param$bandwidth, bg.weight=NULL, 
				tag.shift=.param$tag_shift, step=.param$step, 
				background.density.scaling=.param$bg_density_scaling, 
				rngl=.param$rngl, chrl=.param$chrl, 
				scale.by.dataset.size=.param$scale_by_dataset_size)	

		return(profile)
	}
)

smoothedTagDensity$methods(
	show = function(...) {
##		callSuper(...)
##		if(!is.null(.profile)) {
			cat("Smoothed tag density profile\n")
##			cat("  ", sum(unlist(lapply(.profile, nrow))), " data points\n")
##		}
	}
)




















