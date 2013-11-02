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
		known_params <- intersect(names(param), names(.param))
		unknown_params <- setdiff(names(param), names(.param))
		if(length(unknown_params) > 0)
			warning("These parameters are unknown: ", 
				paste(unknown_params, collapse=', '))
		if(length(known_params) > 0) 
			.param[known_params] <<- param[known_params]
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
        view = function(chr=NULL, start=NULL, end=NULL, col_sig="red", 
			col_bg="green", ...) {
		if(!is(.self, "smoothedEnrich") && !is(.self, "conservEnrich")
			&& !is(.self, "smoothedTagDensity"))
		stop("This function only supports objects of class 'smoothedEnrich', 'conservEnrich' and 'smoothedTagDensity!'")
		if(is.null(chr) || is.null(start) || is.null(end)) 
			stop("Please specify 'chr', 'start' and 'end'!")
		##!dirty code
		##cache chrl and rngl
		temp.chrl <- .param$chrl
		temp.rngl <- .param$rngl
		##temporarily set chrl and rngl
		.param$chrl <<- chr
		.param$rngl <<- list(c(start, end))
		names(.param$rngl)[1] <<- names(.param$chrl) <<- chr
		temp_profile <- .self$get.profile()
		##revert chrl and rngl
		.param$chrl <<- temp.chrl
		.param$rngl <<- temp.rngl
##		dev.new(width=16, height=2.5)
		par(mar=c(4, 2.5, 1, 1))
		plot(temp_profile[[chr]][, 1], temp_profile[[chr]][, 2], type='h', 
			col = ifelse(temp_profile[[chr]][, 2] >= 0, col_sig, col_bg), 
			xlab=paste(chr, ":", format(start, scientific=F), 
			"-", format(end, scientific=F), sep=""), ylab="", ...)
})



















