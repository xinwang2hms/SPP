##class for point binding position profile
bindingPos = setRefClass(
	Class = "bindingPos",
	fields = list(
		.param.updated = "logical"
	), 
	contains = "ChIPSeqProfile"
)
bindingPos$methods(
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
			.param <<- list(
				f=1, e_value=NULL, fdr=NULL, masked_data=NULL, whs=NULL, 
				min_dist=200, window_size=4e7, nrand=3, 
				shuffle_window=1, min_thr=2, topN=NULL, tag_count_whs=100, 
				enrichment_z=2, method="tag.wtd", tec_filter=TRUE, 
				tec_window_size=1e4, tec_z=5, tec_masking_window_size=1e4,
				tec_poisson_z=5, tec_poisson_ratio=5, tec=NULL, 
				n_ctr_samples=1, enrichment_scale_down_ctr=FALSE, 
				enrichment_bg_scales=c(1,5,10), use_rand_ctr=FALSE, 
				bg_density_scaling=TRUE, mle_filter=FALSE, 
				min_mle_threshold=1, add_broad_peak_reg=FALSE, 
				broad_peak_reg_window_size=500, broad_peak_reg_z_thr=2, 
				chrl=NULL)
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
			##chrl		
			##determine common range
			if(is.null(.param$chrl) || (!is.list(.param$chrl) && .param$chrl=="all")) {
##				cat("setting chromosomes: chrl ")
				.param$chrl <<- intersect(names(.ChIP$tags), names(.Input$tags))
				names(.param$chrl) <<- .param$chrl
##				cat("[done]\n")
			}
			##!
			if(is.null(.param$whs)) {
##				cat("setting window half size: bg_weight ")
				.param$whs <<- .ChIP$bd_chrtcs$whs
##				cat("=", round(.param$whs, 2), "[done]\n")
			}
		}
			.param.updated <<- TRUE
		##!							
		} 
		##--copied object
	}
)

bindingPos$methods(
	set.param = function(..., verbose=TRUE) {
		callSuper(..., verbose=TRUE)
		##method
		if(! .param$method%in%c("tag.wtd", "tag.lwcc"))
			stop("'method' should be either 'tag.wtd' or 'tag.lwcc'!")
#		##rngl
#		if(is.null(.param$rngl) || (!is.list(.param$rngl) && .param$rngl=="all")) {
#			cat("setting chromosome ranges: rngl ")
#			.param$rngl <<- lapply(.param$chrl,function(chr) 
#				range(c(range(abs(.ChIP$tags[[chr]]+.param$tag_shift)), 
#				range(abs(.Input$tags[[chr]]+.param$tag_shift)))))
#			cat("[done]\n")
#		} else {
#			.param$chrl <<- names(.param$rngl)
#			names(.param$chrl) <<- .param$chrl
#		}
		##bg_weight
#		if(is.null(.param$bg_weight)) {
#			cat("estimating background weight: bg_weight ")
#			.param$bg_weight <<- dataset.density.ratio(.ChIP$tags, 
#				.Input$tags, background.density.scaling=.param$bg_density_scaling)
#			cat("=", round(.param$bg_weight, 2), "[done]\n")
#		}		
		.param.updated <<- TRUE
	}
)


##6. identify point binding positions
##!n.control.samples=1 not found in any low-level function
bindingPos$methods(
	identify = function(..., verbose=FALSE) {
	
		#!waste memory, need to improve
		signal.data <- .ChIP$tags
		if(!is.null(.Input))
			control.data <- .Input$tags
		else 
			control.data <- NULL
		#!

		if(.param$f<1) {
##			if(verbose) 
				cat("subsampling signal ... ")
			signal.data <- lapply(signal.data, function(x) 
				sample(x, length(x) * .param$f))
##			if(verbose) 
				cat("[done]\n")
		}	
		if(!is.null(control.data) && !.param$use_rand_ctr) {
			# limit both control and signal data to a common set of chromosomes
##			chrl <- intersect(names(signal.data), names(control.data))
			signal.data <- signal.data[.param$chrl]
			control.data <- control.data[.param$chrl]
			control <- list(control.data)
		} else {
			control <- NULL
		}	
		#method
		if(.param$method=="tag.wtd")
			bp.method <- tag.wtd		
		else if(.param$method=="tag.lwcc")
			bp.method <- tag.lwcc
		cat("Identifying point binding positions ... ")
		.profile <<- lwcc.prediction(signal.data, min.dist=.param$min_dist, 
			whs=.param$whs, window.size=.param$window_size, 
			e.value=.param$e_value, fdr=.param$fdr, debug=verbose, 
			n.randomizations=.param$nrand, shuffle.window=.param$shuffle_window, 
			min.thr=.param$min_thr, 
			method=bp.method, tag.shift=.param$tag_shift, 
			bg.tl=control.data, mask.tl=.param$masked_data, topN=.param$topN, 
			control=control, tec.filter=.param$tec_filter, tec.z=.param$tec_z, 
			tec.window.size=.param$tec_window_size, tec.masking.window.size=
			.param$tec_masking_window_size, tec.poisson.z=.param$tec_poisson_z, 
			tec.poisson.ratio=.param$tec_poisson_ratio, 
			background.density.scaling=.param$bg_density_scaling)
		
		# add tag counts
##		chrl <- names(profile$npl)
##		names(chrl) <- chrl
		.profile$npl <<- lapply(.param$chrl,function(chr) {
			pd <- .profile$npl[[chr]]
			pd$nt <- points.within(abs(signal.data[[chr]]), 
				pd$x-.param$tag_count_whs, pd$x+.param$tag_count_whs, 
				return.point.counts=T)
			return(pd)
		})
		.profile$f <<- .param$f
		.profile$n <<- sum(unlist(lapply(signal.data,length)))
		if(!is.null(control.data)) {
			.profile$n.bg <<- sum(unlist(lapply(control.data,length)))
		}
		# calculate enrichment ratios
		.profile <<- calculate.enrichment.estimates(.profile, signal.data, 
			control.data=control.data, fraction=1, 
			tag.count.whs=.param$tag_count_whs, z=.param$enrichment_z, 
			scale.down.control=.param$enrichment_scale_down_ctr, 
			background.scales=.param$enrichment_bg_scales)

		if(.param$mle_filter) {
			if(!is.null(.profile$npl)) {
				if(length(.profile$npl)>1) {
					mle.columns <- grep("enr.mle",colnames(.profile$npl[[1]]))
					if(length(mle.columns)>1) {
						.profile$npl <<- lapply(.profile$npl,function(d) 
							d[apply(d[,mle.columns],1,function(x) 
							all(x>min_mle_threshold)),])
					}
				}
			}
		}

		.profile$whs <<- .param$whs
		cat("[done]\n")

		if(.param$add_broad_peak_reg) {
			cat("Adding broad regions to peaks ... ")
			se <- find.significantly.enriched.regions(signal.data, 
				control.data, window.size=.param$broad_peak_reg_window_size, 
				z.thr=.param$broad_peak_reg_z_thr, poisson.z=0, 
				poisson.ratio=0, either=F, tag.shift=.param$tag_shift)
##			chrl <- names(profile$npl)
##			names(chrl) <- chrl
			.profile$bnpl <<- lapply(.param$chrl,function(chr) {
				npl <- .profile$npl[[chr]]
				if(is.null(npl) || dim(npl)[1]<1) {
					return(npl)
				}
				pi <- points.within(npl$x, se[[chr]]$s, se[[chr]]$e, 
					return.list=T)
				pm <- do.call(rbind, lapply(pi, function(rl) {
					if(length(rl)>0) {
						return(range(c(se[[chr]]$s[rl], se[[chr]]$e[rl])))
					} else {
						return(c(NA, NA))
					}
				}))
				npl$rs <- pm[, 1]
				npl$re <- pm[, 2]
				return(npl)
			})		
			cat("[done]\n")		
		}
		##add parameter list to .profile
##		.profile$param <<- .param
		.param.updated <<- FALSE	
	}
)

bindingPos$methods(
	get.profile = function(..., verbose=FALSE) {
		##call 'identify' if the profile does not exist
		if(is.null(.profile) || .param.updated)
			.self$identify(verbose=verbose)
		return(.profile)
	}
)

##get binding positions
bindingPos$methods(
	get.bp = function(sort_by="chr", ...) {

		##check sort.by
		if(! sort_by %in% c("chr", "score", "Evalue", "FDR"))
			stop("'sort_by' should be one of 'chr', 'score', 'Evalue' and 'FDR'!")
		##call 'get' -> 'identify' to get the profile if not exist
		if(is.null(.profile) || .param.updated)
			.self$identify()
		chrl <- names(.profile$npl) 
		names(chrl) <- chrl
		x <- do.call(rbind, lapply(chrl, function(chr) {
			d <- .profile$npl[[chr]]
			if(dim(d)[1]>0) {
				if(.profile$thr$type=="topN") {
					od <- cbind(rep(chr,dim(d)[1]), 
						subset(d,select=c(x,y,enr,enr.mle)))
				} else {
					od <- cbind(rep(chr,dim(d)[1]), 
						subset(d,select=c(x,y,evalue,fdr,enr,enr.mle)))
				}
				return(od)
			}			
		}))
		colnames(x) <- c("chr", "pos", "score", "Evalue", "FDR", 
			"enrichment.lb", "enrichment.mle")
		if(sort_by=="score")
			x <- x[sort.list(x[, "score"], decreasing=T), ]
		if(sort_by=="Evalue")
			x <- x[sort.list(x[, "Evalue"], decreasing=T), ]
		if(sort_by=="FDR")
			x <- x[sort.list(x[, "FDR"], decreasing=F), ]			
		if(sort_by=="chr") {
			x <- x[sort.list(x[, "pos"], decreasing=F), ]
			x <- x[sort.list(as.character(x[, "chr"])), ]
		}
		rownames(x) <- NULL
		return(x)	
	}
)
##write binding positions
bindingPos$methods(
	write.bp = function(file, sort_by="chr", ...) {
	
		x <- .self$get.bp(sort_by=sort_by, ...)
		write.table(x, file=file, col.names=T, row.names=F, 
			sep="\t", append=F, quote=F)		
	}
)
bindingPos$methods(
        view.bp = function(chr, start=0, end=Inf, col="red", ...) {    
        ##check 'chr'
		if(missing(chr)) 
			stop("Please specify 'chr'!")
		if(is.null(.Input))
			tmp.chrl <- names(.ChIP$tags)
		else 
			tmp.chrl <- intersect(names(.ChIP$tags), names(.Input$tags))			
		if(! (chr %in% tmp.chrl))
			stop(paste("No tags in '", chr, "'",  sep=""))
		##check 'start' and 'end'
		if(!is.numeric(start) || !is.numeric(end))
			stop("'start' and 'end' should be numeric, or 'Inf' indicating the end of chromosome!")			
		if(start > end)
			stop("start should be < stop!")
		##compute profile if not yet exist
		if(is.null(.profile) || .param.updated)
			.self$identify()
		##decode input start and end
		if(is.null(.Input))
			tmp.rng <- range(abs(.ChIP$tags[[chr]]+.param$tag_shift))				
		else
			tmp.rng <- range(c(range(abs(.ChIP$tags[[chr]]+.param$tag_shift)), 
				range(abs(.Input$tags[[chr]]+.param$tag_shift))))
		rng.start <- max(start, tmp.rng[1])
		rng.end <- min(end, tmp.rng[2])
		
		if(! chr %in% names(.profile$npl))
			stop(paste("No significant binding positions in ", chr, 
				"!", sep=""))
		##retrieve from .profile
		inds <- which(points.within(.profile$npl[[chr]][, "x"], fs=rng.start, fe=rng.end)==1)
		if(length(inds)==0)
			stop("No significant binding positions in specified region!")		

		##limit maximal number of points by random sampling
		##?maybe fix 'randomness' by set.seed?
		density.max.points <- getOption("density.max.points")
		if(is.null(density.max.points))
			density.max.points <- 1e4
		if(length(inds) > density.max.points)	
##			inds <- sample(inds, density.max.points, replace=FALSE)
			inds <- sample(inds, density.max.points, replace=FALSE, 
				prob=.profile$npl[[chr]][inds, "y"])
		##finally, plot
		par(mar=c(4, 2.5, 1, 1))
		if(start==0)
			start <- rng.start
		if(end==Inf)
			end <- rng.end
		plot(.profile$npl[[chr]][inds, "x"], .profile$npl[[chr]][inds, "y"], 
			type='h', col = col, xlab=paste(chr, ":", 
			format(start, scientific=F), "-", format(end, scientific=F), 
			sep=""), ylab="", xlim=c(start, end), ...)
})
##get table in narrowPeak format 
bindingPos$methods(
	get.narrowpeak = function(sort_by="chr", ...) {
		##check sort.by
		if(! sort_by %in% c("chr", "score", "FDR"))
			stop("'sort_by' should be one of 'chr', 'score', and 'FDR'!")
		if(is.null(.profile) || .param.updated) {
			if(!.param$add_broad_peak_reg)
				.param$add_broad_peak_reg <<- TRUE
			.self$identify()	
		}	
		if(is.null(.param$whs)) 
			margin <- 50
		else 
			margin <- .param$whs
		chrl <- names(.profile$npl)
		names(chrl) <- chrl
##		chrl <- sort(chrl)

		md <- do.call(rbind, lapply(chrl,function(chr) {
			df <- .profile$npl[[chr]]
			x <- df$x
			rs <- df$rs
			if(is.null(rs)) 
				rs <- rep(NA,length(x))
			re <- df$re 
			if(is.null(re)) 
				re <- rep(NA,length(x)) 
			ivi <- which(is.na(rs))
			if(any(ivi)) 
				rs[ivi] <- x[ivi]-margin
			ivi <- which(is.na(re)) 
			if(any(ivi)) 
				re[ivi] <- x[ivi]+margin
			data.frame(chr=chr, start=rs, end=re, score=df$y, FDR=df$fdr)
##			cbind(chr, rs, re, ".", "0", ".", df$y, -1, 
##				format(df$fdr, scientific=T, digits=3), x-rs)
		}))
		if(sort_by=="chr") {
			md <- md[sort.list(md[, "end"], decreasing=F),]
			md <- md[sort.list(md[, "start"], decreasing=F),]
			md <- md[sort.list(md[, "chr"], decreasing=F),]
		} else if(sort_by=="score") {
			md <- md[sort.list(md[, "score"], decreasing=T), ]
		} else if(sort_by=="FDR") {
			md <- md[sort.list(md[, "FDR"], decreasing=F), ]
		}
		rownames(md) <- NULL
		return(md)		
	}

)
##write to a table in narrowPeak format
bindingPos$methods(
	write.narrowpeak = function(file, ...) {
##		md <- .self$get.narrowpeak(...)
		if(is.null(.profile) || .param.updated) {
			if(!.param$add_broad_peak_reg)
				.param$add_broad_peak_reg <<- TRUE
			.self$identify()	
		}	
		if(is.null(.param$whs)) 
			margin <- 50
		else 
			margin <- .param$whs
		chrl <- names(.profile$npl)
		names(chrl) <- chrl
		md <- do.call(rbind, lapply(chrl,function(chr) {
			df <- .profile$npl[[chr]]
			x <- df$x
			rs <- df$rs
			if(is.null(rs)) 
				rs <- rep(NA,length(x))
			re <- df$re 
			if(is.null(re)) 
				re <- rep(NA,length(x)) 
			ivi <- which(is.na(rs))
			if(any(ivi)) 
				rs[ivi] <- x[ivi]-margin
			ivi <- which(is.na(re)) 
			if(any(ivi)) 
				re[ivi] <- x[ivi]+margin
			cbind(chr, rs, re, ".", "0", ".", df$y, -1, 
				format(df$fdr, scientific=T, digits=3), x-rs)
		}))
		md <- md[order(as.numeric(md[,7]),decreasing=T),]


		write.table(md, file=file, col.names=F, row.names=F, quote=F, 
			sep="\t", append=F)
	}
)


bindingPos$methods(
	show = function(...) {
		##message about Input
##		callSuper(...)
##		if(!is.null(.profile)) {
##			if(!is.null(profile$npl)) {
##				cat("Binding position profile:\n")
##				cat("  ", sum(unlist(lapply(.profile$npl, nrow))), 
##					" significant binding positions\n")
##			}
##			if(!is.null(.profile$bnpl)) {
##				cat("  ", sum(unlist(lapply(.profile$bnpl, function(x) 
##					sum(!is.na(x[, "rs"]))))), "broad peak regions\n")
##			}
##		}
		cat("Point binding position profile")
	}
)























