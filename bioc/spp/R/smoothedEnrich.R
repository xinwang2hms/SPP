##class representing smoothed enrichment profile
smoothedEnrich = setRefClass(
	Class = "smoothedEnrich",
##	fields = list(
##
##	), 
	contains = "ChIPSeqProfile"
)
smoothedEnrich$methods(
	initialize = function(..., ChIP=NULL, Input=NULL, param=NULL) {
		##check ChIP and Input
		if(!is(ChIP, "AlignedTags_Or_NULL") || 
			!is(Input, "AlignedTags_Or_NULL"))
			stop("'ChIP' and 'Input' should be AlingedTags objects or 'NULL'!")	
		callSuper(..., ChIP=ChIP, Input=Input)
		##3. set .param		
		if(!is(param, "list_Or_NULL"))
			stop("'param' should be a list of parameters or 'NULL'!")
		##--new object
		if(is.null(.param))	{
			##initialize default .param
			##!to delete: scale_by_dataset_size and npseudo
			.param <<- list(bandwidth=150, bg_weight=NULL, tag_shift=NULL, 
				step=50, bg_density_scaling=TRUE, rngl=NULL, chrl=NULL, 
				scale_by_dataset_size=F, npseudo=1)
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
				.param$chrl <<- intersect(names(.ChIP$tags), names(.Input$tags))
				names(.param$chrl) <<- .param$chrl
##				cat("[done]\n")
			}
			##rngl
			if(is.null(.param$rngl) || (!is.list(.param$rngl) && .param$rngl=="all")) {
##				cat("setting chromosome ranges: rngl ")
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
smoothedEnrich$methods(
	set.param = function(..., verbose=TRUE) {
		callSuper(..., verbose=TRUE)
	}
)
##4. compute smoothed enrichment by MLE

smoothedEnrich$methods(
	get.profile = function() {
		##bg_weight
		if(is.null(.param$bg_weight)) {
			cat("estimating background weight: bg_weight ")
			.param$bg_weight <<- dataset.density.ratio(.ChIP$tags, 
				.Input$tags, background.density.scaling=.param$bg_density_scaling)
			cat("=", round(.param$bg_weight, 2), "[done]\n")
		}
							
		ssd <- get.smoothed.tag.density(.ChIP$tags, 
			bandwidth=.param$bandwidth, bg.weight=.param$bg_weight, 
			tag.shift=.param$tag_shift, step=.param$step, 
			background.density.scaling=.param$bg_density_scaling, 
			rngl=.param$rngl, chrl=.param$chrl, 
			scale.by.dataset.size=.param$scale_by_dataset_size)
		csd <- get.smoothed.tag.density(.Input$tags, 
			bandwidth=.param$bandwidth, bg.weight=.param$bg_weight, 
			tag.shift=.param$tag_shift, step=.param$step, 
			background.density.scaling=.param$bg_density_scaling, 
			rngl=.param$rngl, chrl=.param$chrl, 
			scale.by.dataset.size=.param$scale_by_dataset_size)

	
		profile <- lapply(.param$chrl, function(chr) {
			d <- ssd[[chr]]
			d$y <- log2(d$y+.param$npseudo)-log2(csd[[chr]]$y+ 
				.param$npseudo)-log2(.param$bg_weight)
			return(d)
		})	
		return(profile)
	}
)

smoothedEnrich$methods(
	show = function(...) {
		callSuper(...)
		if(!is.null(.profile)) {
			cat("Smoothed enrichment profile:\n")
##			cat("  ", sum(unlist(lapply(.profile, nrow))), " data points\n")
		}
	}
)




















