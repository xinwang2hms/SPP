##A root class representing ChIP-seq profile

ChIPSeqProfile = setRefClass(
	Class = "ChIPSeqProfile",
	fields = list(
		.ChIP = "ANY", 
		.Input = "ANY",
		.profile = "list_Or_NULL",
		.param = "list_Or_NULL"			##parameters
	)
)
ChIPSeqProfile$methods(
	initialize = function(..., ChIP=NULL, Input=NULL) {
		.ChIP <<- NULL
		.Input <<- NULL
		callSuper(...)
	
		##1. set .ChIP
		##--new object
		if(is.null(.ChIP)) {
			if(!is.null(ChIP))
				.ChIP <<- ChIP
		} 
		##--otherwise, copied object
		##2. set .Input
		##--new object
		if(is.null(.Input)) {
			if(!is.null(Input))	
				.Input <<- Input
		}
		##--otherwise, copied object	
		##3. set .param
		if(is.null(.param))
			.param <<- NULL
		##4. set .profile
		if(is.null(.profile))
			.profile <<- NULL			
	}
)
##set reference to ChIP
ChIPSeqProfile$methods(
	set.ChIP = function(ChIP) {
		if(missing(ChIP))
			stop("'ChIP' should be an AlingedTags object")
		.ChIP <<- ChIP
	}
)
##set reference to Input
ChIPSeqProfile$methods(
	set.Input = function(Input) {
		if(missing(Input))
			stop("If 'Input' is provided, it should be an AlingedTags object")
		.Input <<- Input
	}
)
##set parameters
ChIPSeqProfile$methods(
	set.param = function(..., verbose=TRUE) {
		param <- list(...)
		known_params <- setdiff(intersect(names(param), names(.param)), c("rngl", "chrl"))
		unknown_params <- setdiff(names(param), names(.param))
		if(length(unknown_params) > 0)
			warning("These parameters are unknown: ", 
				paste(unknown_params, collapse=', '))
		if(length(known_params) > 0) 
			.param[known_params] <<- param[known_params]
		##set rngl and chrl if available
		if(all(c("chrl", "rngl") %in% names(param)))
			stop("Please do not set 'chrl' and 'rngl' at the same time, which might cause confusion!")
		if("rngl" %in% names(param)) {
			rngl <- param[["rngl"]]
			if(length(rngl)==0)
				stop("Empty 'rngl'!")
			##rngl
			if(is.null(.Input))
				tmp.chrl <- names(.ChIP$tags)
			else 
				tmp.chrl <- intersect(names(.ChIP$tags), names(.Input$tags))
			names(tmp.chrl) <- tmp.chrl
			if(!is.list(rngl) && rngl!="all")
				stop("'rngl' should be either a list of genomic ranges or 'all' indicating all regions!")	
			if(!is.list(rngl) && rngl=="all") {
				if(is.null(.Input))
					tmp.rngl <- lapply(tmp.chrl,function(chr) 
						range(abs(.ChIP$tags[[chr]]+.param$tag_shift)))				
				else
					tmp.rngl <- lapply(tmp.chrl,function(chr) 
						range(c(range(abs(.ChIP$tags[[chr]]+.param$tag_shift)), 
						range(abs(.Input$tags[[chr]]+.param$tag_shift)))))
			} else if(is.list(rngl)) {
				tmp.chrl <- intersect(names(rngl), tmp.chrl)
				names(tmp.chrl) <- tmp.chrl
				if(is.null(.Input))
					tmp.rngl <- lapply(tmp.chrl,function(chr) 
						range(abs(.ChIP$tags[[chr]]+.param$tag_shift)))				
				else
					tmp.rngl <- lapply(tmp.chrl,function(chr) 
						range(c(range(abs(.ChIP$tags[[chr]]+.param$tag_shift)), 
						range(abs(.Input$tags[[chr]]+.param$tag_shift)))))
				tmp.rngl <- lapply(tmp.chrl, function(chr) {
					if(length(rngl[[chr]])!=2)
						stop("Please provide and only provide 'start' and 'end' in 'rngl'")
					if(rngl[[chr]][1]>=rngl[[chr]][2])
						stop("'start' should always be smaller than 'end' in 'rngl'")
					c(max(tmp.rngl[[chr]][1], rngl[[chr]][1]), 
					min(tmp.rngl[[chr]][2], rngl[[chr]][2]))
				})
			} 
			.param$rngl <<- tmp.rngl
			.param$chrl <<- tmp.chrl
		} else if("chrl" %in% names(param)) {
			chrl <- param[["chrl"]]
			##chrl	
			if(is.null(.Input))
				tmp.chrl <- names(.ChIP$tags)
			else 
				tmp.chrl <- intersect(names(.ChIP$tags), names(.Input$tags))
			if(!is.character(chrl))	
				stop("'chrl' should be a vector of chromosome names or 'all'!")
			if(length(chrl)>1) {
				tmp.chrl <- intersect(tmp.chrl, chrl)
			}
			if(length(tmp.chrl)==0)
				stop("'chrl' not available in the data!") 
			names(tmp.chrl) <- tmp.chrl
			##infer rngl from chrl
			if(is.null(.Input))
				tmp.rngl <- lapply(tmp.chrl,function(chr) 
					range(abs(.ChIP$tags[[chr]]+.param$tag_shift)))				
			else
				tmp.rngl <- lapply(tmp.chrl,function(chr) 
					range(c(range(abs(.ChIP$tags[[chr]]+.param$tag_shift)), 
					range(abs(.Input$tags[[chr]]+.param$tag_shift)))))
			.param$rngl <<- tmp.rngl
			.param$chrl <<- tmp.chrl
		} 
	}
)
ChIPSeqProfile$methods(
	get.param = function(what) {
		if(!(what %in% c(names(.param), "all")))
			stop("No such parameter")
		if(what=="all")
			return(.param)
		else 
			return(.param[[what]])
	}
)


ChIPSeqProfile$methods(
	write.wig = function(file, name, feature, threshold=5, zip=F) {
		if(!is(.self, "smoothedEnrich") && !is(.self, "conservEnrich")
			&& !is(.self, "smoothedTagDensity"))
		stop("This function only supports objects of class 'smoothedEnrich', 'conservEnrich' and 'smoothedTagDensity!'")
		dat <- .self$get.profile()
		invisible(lapply(.param$chrl,function(chr) {
			bdiff <- dat[[chr]]
			ind <- seq(1,length(bdiff$x))
			ind <- ind[!is.na(bdiff$y[ind])]
			header <- chr==.param$chrl[1]
			write.probe.wig(chr,bdiff$x[ind],bdiff$y[ind],file,append=!header,
				feature=feature, header=header, name=name)
		}))
		if(zip) {
			zf <- paste(file,"zip",sep=".")
			tryCatch(system(paste("zip \"",zf,"\" \"",file,"\"",sep="")), 
				error=function(e) {print(e)})
			tryCatch(system(paste("rm \"",file,"\"",sep="")), 
				error=function(e) {print(e)})
			cat(zf, "\n")
		} else {
			cat(file, "\n")
		}
	}
)

ChIPSeqProfile$methods(
	write.tdf = function(file, name, feature, save_wig=F, zip_wig=T) {
		if(!is(.self, "smoothedEnrich") && !is(.self, "conservEnrich")
			&& !is(.self, "smoothedTagDensity"))
		stop("This function only supports objects of class 'smoothedEnrich', 'conservEnrich' and 'smoothedTagDensity!'")
		igvtools <- getOption("igvtools")
		if(is.null(igvtools))
			stop("Please specify 'igvtools', which is the path of IGVTools")
		tfile <- paste(file, "wig", sep=".")
		.self$write.wig(tfile, feature=feature, name=name, zip=F)
		paste(igvtools, "-f min,max,mean", tfile, file)
		paste(igvtools,"toTDF -f min,max,mean",tfile,file,.ChIP$genome_build)
		tryCatch(system(paste(igvtools,"toTDF -f min,max,mean",tfile,file,.ChIP$genome_build)))
		if(!save_wig) {
			tryCatch(system(paste("rm",tfile)), error=function(e) {print(e)})
		} else if(zip_wig) {
			zf <- paste(file,"zip",sep=".")
                        tryCatch(system(paste("zip \"",zf,"\" \"",file,"\"",sep="")),
                                error=function(e) {print(e)})
                        tryCatch(system(paste("rm \"",file,"\"",sep="")),
                                error=function(e) {print(e)})
                        cat(zf, "\n")
		}
	}
)


ChIPSeqProfile$methods(
	show = function(...) {	

#		##message about Input
#		if(!is.null(.Input)) {
#			cat("Input:\n")
#			if(!is.null(.Input$tags)) {
#				nchrs <- length(.Input$tags)
#				ntags <- sum(unlist(lapply(.Input$tags, length)))
#				cat(paste("  ", ntags, " fragments", " across ", nchrs, 
#					" chromosome(s)", "\n", sep=""))			
#			}		
#		}
#		##message about ChIP
#		if(!is.null(.ChIP)) {
#			if(!is.null(.ChIP$tags)) {
#				nchrs <- length(.ChIP$tags)
#				ntags <- sum(unlist(lapply(.ChIP$tags, length)))
#				cat(paste("  ", ntags, " fragments", " across ", nchrs, 
#					" chromosome(s)", "\n", sep=""))			
#			}
#			##show binding characteristics
#			if(!is.null(.ChIP$bd_chrtcs)) {
#				##
#				cat("Binding characteristics:\n")
#				cat(paste("  cross correlation peak: Position=", 
#					.ChIP$bd_chrtcs$peak$x, ", Height=", 
#					round(.ChIP$bd_chrtcs$peak$y, 3), "\n", sep=""))
#				cat(paste("  optimized window half-size: ", 
#					.ChIP$bd_chrtcs$whs, "\n", sep=""))
#			}							
#		}
		##message about parameters
		if(!is.null(.param)) {
##			cat("Parameters:\n")
##			cat(paste("  ", paste(paste(names(.param), .param, sep="="), 
##				collapse=', '), "\n", sep=""))
		}	
	}
)
ChIPSeqProfile$methods(
	get.profile = function(...) {
	return(.profile)
})

ChIPSeqProfile$methods(
        view = function(chr, start=0, end=Inf, col_sig="red", 
			col_bg="green", ...) {
		if(!is(.self, "smoothedEnrich") && !is(.self, "conservEnrich")
			&& !is(.self, "smoothedTagDensity"))
		stop("This function only supports objects of class 'smoothedEnrich', 'conservEnrich' and 'smoothedTagDensity!'")
		if(missing(chr)) 
			stop("Please specify 'chr'!")
		##check start end 
		if(! (chr %in% .param$chrl))
			stop(paste("No tags in '", chr, "'",  sep=""))
		if(!is.numeric(start) || !is.numeric(end))
			stop("'start' and 'end' should be numeric, or 'Inf' indicating the end of chromosome!")
		##!dirty code
		##cache chrl and rngl
		temp.chrl <- .param$chrl
		temp.rngl <- .param$rngl
		##temporarily set chrl and rngl
##		.param$chrl <<- chr
##		.param$rngl <<- list(c(start, end))
		rngl <- list(c(start, end))
		names(rngl) <- chr
		.self$set.param(rngl=rngl)
		start <- .param$rngl[[1]][1]
		end <- .param$rngl[[1]][2]
		density.max.points <- getOption("density.max.points")
		if(is.null(density.max.points))
			density.max.points <- 1e4
		min.step <- round((end-start+1)/density.max.points)
		temp.step <- NULL
		if((end-start+1)/.param$step > density.max.points) {
			temp.step <- .param$step
			.param$step <<- min.step
			warning(paste("set 'step' = ", min.step, " to avoid overuse of memory!", sep=""))
		}

##		names(.param$rngl)[1] <<- names(.param$chrl) <<- chr
		temp_profile <- .self$get.profile()
		##revert chrl and rngl
		.param$chrl <<- temp.chrl
		.param$rngl <<- temp.rngl
		if(!is.null(temp.step))
			.param$step <<- temp.step
##		dev.new(width=16, height=2.5)
		par(mar=c(4, 2.5, 1, 1))
		plot(temp_profile[[chr]][, 1], temp_profile[[chr]][, 2], type='h', 
			col = ifelse(temp_profile[[chr]][, 2] >= 0, col_sig, col_bg), 
			xlab=paste(chr, ":", format(start, scientific=F), 
			"-", format(end, scientific=F), sep=""), ylab="", ...)
})



















