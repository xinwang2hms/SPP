##Class BAMTags is inherited from AlignedTags for BAM alignment results
##!issue: may need to set specific slots for BAM results
BAMTags = setRefClass(
	Class = "BAMTags", 
	contains = "AlignedTags"
)
##initialization function for BAMTags

BAMTags$methods(
	read = function(file=NULL, read_tag_names=F, read_tag_qualities=F, fix_chr_names=F, ...) {
		if(is.null(file) | !file.exists(file))
			stop("Please specify 'file', which is the BAM alignment file!")
		else 
			file <<- file
		if(!is.logical(read_tag_names))
			stop("'read_tag_names' should be a logical value!")
		if(!is.logical(fix_chr_names))
			stop("'fix_chr_names' should be a logical value!")	
		if(read_tag_names) 
			rtn <- as.integer(1)
		else
			rtn <- as.integer(0)
		tl <- lapply(.Call("read_bam", path.expand(file), rtn), 
				function(d) {
					xo <- order(abs(d$t))
					d$t <- d$t[xo]
					if(read_tag_qualities)
						d$n <- d$n[xo]
					if(read_tag_names) 
						d$s <- d$s[xo]
					return(d)
				})
		if(fix_chr_names) {
			# remove ".fa"
			names(tl) <- gsub("\\.fa","",names(tl))
		}
		# separate tags and quality
		tags <<- lapply(tl,function(d) d$t)
		if(read_tag_qualities)
			quality <<- lapply(tl,function(d) d$n)
		if(read_tag_names)
			names <<- lapply(tl,function(d) d$s)
	}
)

##BAMTags$methods(
##	show = function(...) {
##		callSuper(...)
##	}
##)






