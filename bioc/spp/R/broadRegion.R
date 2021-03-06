##class for broad region profile
broadRegion = setRefClass(
	Class = "broadRegion",
	fields = list(
		.param.updated = "logical"
	), 
	contains = "ChIPSeqProfile"
)
broadRegion$methods(
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
			##!to delete: multiplier, z_thr, mcs, mask_window_size, poisson_z, 
			##poisson_ratio, either
			.param <<- list(window_size=500, tag_shift=NULL, bg_weight=NULL, 
				bg_density_scaling=TRUE, multiplier=1, z_thr=3, mcs=0, 
				mask_window_size=500, poisson_z=0, rngl=NULL, 
				chrl=NULL, poisson_ratio=4, either=FALSE)
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
		if(!is.null(.param)) {
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
			.param.updated <<- TRUE
		} 
		##--copied object			
	}
)

broadRegion$methods(
	update = function() {
		invisible(broadRegion(.self))
	}
)
broadRegion$methods(
	size = function() {
		s <- object.size(.param) + object.size(.profile) + object.size(.param.updated)
		return(s)
	}
)


broadRegion$methods(
	set.param = function(..., verbose=TRUE) {
		callSuper(..., verbose=TRUE)
		.param.updated <<- TRUE
	}
)

broadRegion$methods(
	write.broadpeak = function(file, ...) {
	
	if(is.null(.profile) || .param.updated)
		.self$identify()
	chrl <- names(.profile)
	names(chrl) <- chrl
	chrl <- chrl[unlist(lapply(.profile,function(d) length(d$s)))>0]
	md <- do.call(rbind,lapply(chrl,function(chr) {
		df <- .profile[[chr]]
		cbind(chr, df$s, df$e, ".", "0", ".", df$rv, -1, -1)
	}))
	md <- md[order(as.numeric(md[,7]),decreasing=T),]
	write.table(md,file=file,col.names=F,row.names=F,quote=F,sep="\t",append=F)
		
	}
)

broadRegion$methods(
        view = function(chr, start=0, end=Inf, col_sig="red", 
			col_bg="green", ...) {
		##check 'chr'	
		if(missing(chr)) 
			stop("Please specify 'chr'!")
		if(is.null(.Input))
			tmp.chrl <- names(.ChIP$tags)
		else 
			tmp.chrl <- intersect(names(.ChIP$tags), names(.Input$tags))			
		if(! (chr %in% tmp.chrl))
			stop(paste("No tags in '", chr, "'",  sep=""))
##		if(! (chr %in% .param$chrl))
##			stop(paste("No tags in '", chr, "'", sep=""))
		##check 'start' and 'end'
		if(!is.numeric(start) || !is.numeric(end))
			stop("'start' and 'end' should be numeric, or 'Inf' indicating the end of chromosome!")			
		if(start > end)
			stop("start should be < stop!")
		##!dirty code
		temp.chrl <- .param$chrl
		temp.rngl <- .param$rngl
##		.param$chrl <<- chr
##		.param$rngl <<- list(c(start, end))
		rngl <- list(c(start, end))
		names(rngl) <- chr
		.self$set.param(rngl=rngl)
		start <- .param$rngl[[1]][1]
		end <- .param$rngl[[1]][2]
##		names(.param$rngl)[1] <<- names(.param$chrl) <<- chr
		if(is.null(.profile) || .param.updated)
			.self$identify()
		.param$chrl <<- temp.chrl
		.param$rngl <<- temp.rngl
		temp_profile <- .profile
		##limit maximal number of rectangles to plot by random sampling	
		density.max.rects <- getOption("density.max.rects")
		if(is.null(density.max.rects))
			density.max.rects <- 500
		inds <- 1:nrow(temp_profile[[chr]])
		if(nrow(temp_profile[[chr]])>density.max.rects)
			inds <- sample(inds, density.max.rects, replace=FALSE, 
			prob=temp_profile[[chr]][, "e"]-temp_profile[[chr]][, "s"]+1)
		
##		dev.new(width=16, height=2.5)
		par(mar=c(4, 2.5, 1, 1))
		plot(temp_profile[[chr]][inds, 1], temp_profile[[chr]][inds, 3], type='n', 
			xlim=c(start, end), xlab=paste(chr, ":", format(start, scientific=F),               
                        "-", format(end, scientific=F), sep=""), ylab="", 
			ylim=c(min(temp_profile[[chr]][, 3], 0), max(temp_profile[[chr]][, 3])))
		trash <- sapply(inds, function(x) {
			rect(temp_profile[[chr]][x, 1], 0, temp_profile[[chr]][x, 2], 
				temp_profile[[chr]][x, 3], lwd=1, 
				border=ifelse(temp_profile[[chr]][x, 3]>0, col_sig, col_bg))
##				 density=2, angle=60,  border=TRUE)
		})
		abline(h=0, col=col_sig)
})
broadRegion$methods(
	identify = function(..., verbose=FALSE) {
		## find significantly enriched clusters
		##!
##		if(is.null(.param$tag_shift))
##		.param$tag_shift <<- .ChIP$bd_chrtcs$peak$x/2
		##!window size??		
##		.param$bg_weight <<- dataset.density.ratio(.ChIP$tags, 
##			.Input$tags, background.density.scaling=.param$bg_density_scaling)


		##bg_weight
		if(is.null(.param$bg_weight)) {
			cat("estimating background weight: bg_weight ")
			.param$bg_weight <<- dataset.density.ratio(.ChIP$tags, 
				.Input$tags, background.density.scaling=.param$bg_density_scaling)
			cat("=", round(.param$bg_weight, 2), "[done]\n")
		}							
		##filter tags
		rngl <- .param$rngl
		
		signal.data <- sapply(1:length(rngl), function(ichr) {
			chr <- names(rngl)[ichr]
			from <- rngl[[chr]][1]
			to <- rngl[[chr]][2]
			inds.ChIP <- which(points.within(abs(.ChIP$tags[[chr]]), from - .param$window_size/2, 
				to + .param$window_size/2 )==1)
			.ChIP$tags[[chr]][inds.ChIP]
		}, simplify=F)
		control.data <- sapply(1:length(rngl), function(ichr) {
			chr <- names(rngl)[ichr]
			from <- rngl[[chr]][1]
			to <- rngl[[chr]][2]
			inds.Input <- which(points.within(abs(.Input$tags[[chr]]), from - .param$window_size/2, 
				to + .param$window_size/2 )==1)
			.Input$tags[[chr]][inds.Input]
		}, simplify=F)
		names(signal.data) <- names(control.data) <- names(rngl)
##		se <- find.significantly.enriched.regions(.ChIP$tags, .Input$tags, 
		se <- find.significantly.enriched.regions(signal.data, control.data, 
			window.size=.param$window_size, multiplier=.param$multiplier, 
			z.thr=.param$z_thr, mcs=.param$mcs, debug=verbose, 
			background.density.scaling=.param$bg_density_scaling, 
			masking.window.size=.param$mask_window_size, 
			poisson.z=.param$poisson_z, poisson.ratio=.param$poisson_ratio, 
			either=.param$either, tag.shift=.param$tag_shift, 
			bg.weight=.param$bg_weight)	
		chrl <- names(se)
		names(chrl) <- chrl
	
		.profile <<- lapply(chrl,function(chr) {
			d <- se[[chr]]
			if(length(d$s>1)) {
				d <- regionset.intersection.c(list(d,d),do.union=T)
				sc <- points.within(abs(.ChIP$tags[[chr]]+.param$tag_shift), 
					d$s, d$e, return.point.counts=T)
				cc <- points.within(abs(.Input$tags[[chr]]+.param$tag_shift), 
					d$s, d$e, return.point.counts=T)
				d$rv <- log2((sc+1)/(cc+1)/.param$bg_weight)
				return(d)
			} else {
				return(d)
			}
		})		
##		.profile$param <<- .param
		.param.updated <<- FALSE
	}
)

broadRegion$methods(
	get.profile = function(sort_by="chr", ...) {
	
	if(! sort_by %in% c("chr", "score"))
		stop("'sort_by' should be one of 'chr' and 'score'")
	if(is.null(.profile)  || .param.updated)
		.self$identify()
	chrl <- names(.profile)
	names(chrl) <- chrl
	chrl <- chrl[unlist(lapply(.profile,function(d) length(d$s)))>0]
	md <- do.call(rbind, lapply(chrl, function(chr) {
		data.frame(chr=rep(chr, nrow(.profile[[chr]])), .profile[[chr]])
	}))
	colnames(md) <- c("chr", "start", "end", "score")
	if(sort_by=="score")
		md <- md[order(as.numeric(md[,4]), decreasing=T), ]
	else if(sort_by=="chr") {
		md <- md[sort.list(md[, 2], decreasing=F), ]	
		md <- md[sort.list(as.character(md[, 1]), decreasing=F), ]	
	}
	rownames(md) <- NULL

	return(md)
}
)
broadRegion$methods(
	show = function(...) {
##		callSuper(...)
##		if(!is.null(.profile)) {
			cat("Broad region profile\n")
##			cat("  ", sum(unlist(lapply(.profile, nrow))), 
##				" broad enrichment clusters\n")
##		}
	}
)






































