##class for conservative enrichment profile
conservEnrich = setRefClass(
	Class = "conservEnrich",
##	fields = list(
##
##	), 
	contains = "ChIPSeqProfile"
)
	
conservEnrich$methods(
	initialize = function(..., ChIP=NULL, Input=NULL, param=NULL) {
		##check ChIP and Input
		if(!is(ChIP, "AlignedTags_Or_NULL") || 
			!is(Input, "AlignedTags_Or_NULL"))
			stop("'ChIP' and 'Input' should be AlingedTags objects or 'NULL'!")	
		callSuper(..., ChIP=ChIP, Input=Input)
#		callSuper(...)
#if(!is.null(.ChIP) && !is(.ChIP, "AlignedTags"))
#	.ChIP <<- NULL
#if(!is.null(.Input) && !is(.Input, "AlignedTags"))
#	.Input <<- NULL
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
			.param <<- list(fws=500, bwsl=c(1, 5, 25, 50)*500, step=100, 
				tag_shift=NULL, alpha=0.05, use_most_informative_scale=FALSE, 
				quick_cal=TRUE, bg_density_scaling=TRUE, bg_weight=NULL, 
				chrl=NULL, rngl=NULL, posl=NULL, return_mle=FALSE)
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
conservEnrich$methods(
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
##compute conservative fold enrichment
##!filtering step of tags should be moved to initialization
conservEnrich$methods(
	get.profile = function() {

		##bg_weight
		if(is.null(.param$bg_weight)) {
			cat("estimating background weight: bg_weight ")
			.param$bg_weight <<- dataset.density.ratio(.ChIP$tags, 
				.Input$tags, background.density.scaling=.param$bg_density_scaling)
			cat("=", round(.param$bg_weight, 2), "[done]\n")
		}
							
		profile <- lapply(.param$chrl, function(chr) {
			if(is.null(.Input$tags[[chr]])) 
				bt <- c() 
			else 
				bt <- abs(.Input$tags[[chr]] + .param$tag_shift)
			if(is.null(.param$posl))
				x <- mbs.enrichment.bounds(abs(.ChIP$tags[[chr]]+.param$tag_shift), 
					bt, fws=.param$fws, bwsl=.param$bwsl, step=.param$step, 
					rng = .param$rngl[[chr]], 
					calculate.upper.bound=T, bg.weight=.param$bg_weight, 
					use.most.informative.scale=.param$use_most_informative_scale, 
					quick.calculation=.param$quick_cal, alpha=.param$alpha)
			else 
				x <- mbs.enrichment.bounds(abs(.ChIP$tags[[chr]]+.param$tag_shift), 
					bt, fws=.param$fws, bwsl=.param$bwsl, step=.param$step, 
					rng = .param$rngl[[chr]], 
					calculate.upper.bound=T, bg.weight=.param$bg_weight, 
					use.most.informative.scale=.param$use_most_informative_scale, 
					quick.calculation=.param$quick_cal, alpha=.param$alpha, pos=posl[[chr]])
			
			# compose profile showing lower bound for enriched, upper bound for depleted regions
			ps <- rep(1,length(x$mle))
			vi <- which(!is.na(x$lb) & x$lb>1)
			ps[vi] <- x$lb[vi]
			vi <- which(!is.na(x$ub) & x$ub<1)
			ps[vi] <- x$ub[vi]
			ps <- log2(ps)
			if(is.null(.param$posl)) {
				if(.param$return_mle) {
					return(data.frame(x=seq(x$x$s,x$x$e,by=x$x$step), 
						y=ps, mle=log2(x$mle), lb=log2(x$lb), ub=log2(x$ub)))
				} else {
					return(data.frame(x=seq(x$x$s,x$x$e,by=x$x$step), y=ps))
				}
			} else {
				if(.param$return_mle) {
					return(data.frame(x=posl[[chr]], y=ps, mle=log2(x$mle), 
						lb=log2(x$lb), ub=log2(x$ub)))
				} else {
					return(data.frame(x=posl[[chr]], y=ps))
				}
			}

		})
		return(profile)			
	}
)


conservEnrich$methods(
	show = function(...) {
##		callSuper(...)
##		if(!is.null(.profile)) {
			cat("Conservative fold enrichment profile\n")
##			cat("  ", sum(unlist(lapply(.profile, nrow))), " data points\n")
##		}
	}
)
















