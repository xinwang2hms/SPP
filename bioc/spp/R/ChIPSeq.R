##Class representing a set of ChIP-seq experiments

ChIPSeq = setRefClass(
	Class = "ChIPSeq", 
	fields = list(
		.sample_info = "df_Or_NULL", 
		.l_aligned_tags = "list_Or_NULL",
		.l_confs = "list_Or_NULL"
	)
)

ChIPSeq$methods(
	show = function(...) {
		cat("~~ChIP-Seq Object~~\n")
		if(!is.null(.l_aligned_tags)) {
			atg.nload <- sum(unlist(sapply(1:nrow(.sample_info), function(x) {
				samp <- as.character(.sample_info[x, "SampleID"])
				if(! samp %in% names(.l_aligned_tags))
					return(FALSE)
				else if(!"atg" %in% names(.l_aligned_tags[[samp]]))
					return(FALSE)
				else if(!is(.l_aligned_tags[[samp]][["atg"]], "AlignedTags"))
					return(FALSE)
				else if(is(.l_aligned_tags[[samp]][["atg"]], "AlignedTags"))
					return(TRUE)
			})))
			cat(atg.nload, " of total ", nrow(.sample_info), " samples loaded\n")
		}
	}
) 

ChIPSeq$methods(
	initialize = function(..., sample_info) {
		##if copy
		callSuper(...)
		##else new
		if(is.null(.sample_info) && !is.null(sample_info)) {
##			if(!file.exists(file))
##				stop("Sample file does not exist!")
##			target <- tryCatch(read.csv(file), error=function(e) print(e))
			if(!all(c("SampleID") %in% colnames(sample_info)))
				stop("Please provide 'SampleID'!")
			if(!any(c("TagFile", "AlignedTagsFile") %in% colnames(sample_info)))
				stop("Please provide 'TagFile' and/or 'AlignedTagsFile'!")
			.sample_info <<- sample_info
			.l_aligned_tags <<- list()
			if("AlignedTagsFile" %in% colnames(sample_info)) {
				sapply(1:nrow(sample_info), function(i) {
					.l_aligned_tags[[as.character(sample_info[i, "SampleID"])]] <<- 
						list(file=as.character(sample_info[i, "AlignedTagsFile"]), 
						atg=NULL)
				})
			} 
##			else {
##				sapply(1:nrow(target), function(i) {
##					.l_aligned_tags[[as.character(target[i, "SampleID"])]] <<- 
##						list(file=NULL, atg=NULL)
##				})
##			}
			
			.l_confs <<- list()
		}
	}
)
##lazy loading
ChIPSeq$methods(
	get.AlignedTags = function(sample, type="BAM", genome_build="mm9", 
		read_tag_names=F, read_tag_qualities=F, fix_chr_names=F) {
		if(! sample %in% names(.l_aligned_tags)) {
		##neither in memory, nor in saved file
			.l_aligned_tags[[sample]] <<- list()
			if(!type %in% c("BAM", "Bowtie"))
				stop("Only 'BAM' and 'Bowtie' tags are supported now!")
			if(type=="BAM") {
				.l_aligned_tags[[sample]][["atg"]] <<- BAMTags(genome_build=genome_build)
			} else if(type=="Bowtie") {
				.l_aligned_tags[[sample]][["atg"]] <<- BowtieTags(genome_build=genome_build)
			}
			.l_aligned_tags[[sample]][["atg"]]$read(read_tag_names=read_tag_names, 
				read_tag_qualities=read_tag_qualities, fix_chr_names=fix_chr_names)
##			cat("AlignedTags object created and not yet saved!\n")
		} else {
			if((! "atg" %in% names(.l_aligned_tags[[sample]])) || is.null(.l_aligned_tags[[sample]][["atg"]])) {
				if(!"file" %in% names(.l_aligned_tags[[sample]]))
					stop("If TagFile is not speficied, AlignedTagsFile should be!")
				if(!file.exists(.l_aligned_tags[[sample]][["file"]]))
					stop(paste("AlignedTags object not exist in ", file, sep=""))
				.l_aligned_tags[[sample]][["atg"]] <<- tryCatch(readRDS(file.path(.l_aligned_tags[[sample]][["file"]])),
					error = function(e) print(e))
			} 
		}
		return(.l_aligned_tags[[sample]][["atg"]])
	}
)

ChIPSeq$methods(
	QC.AlignedTags = function(sample, ...) {
	
	
	} 
)

ChIPSeq$methods(
	save.AlignedTags = function(sample, file) {
		if(! sample %in% names(.l_aligned_tags))
			stop("AlignedTags object not yet loaded in memory!")
		if(missing(file) && is.null(.l_aligned_tags[[sample]][["file"]]))
			stop("Please specify 'file'!")
		if(is.null(.l_aligned_tags[[sample]][["file"]]) && !missing(file))
			.l_aligned_tags[[sample]][["file"]] <<- getAbsolutePath(file)
		if(! is(.l_aligned_tags[[sample]][["atg"]], "AlignedTags") && file.exists(.l_aligned_tags[[sample]][["file"]]))
			stop(paste("AlignedTags object is probably cached already in", 
				.l_aligned_tags[[sample]][["file"]], "?", sep=""))
		ChIP <- .l_aligned_tags[[sample]][["atg"]]
		saveRDS(ChIP, file=.l_aligned_tags[[sample]][["file"]])
		##maybe not needed here
		rm(ChIP)
		.l_aligned_tags[[sample]][["atg"]] <<- NULL
		invisible(gc())
	}
)

ChIPSeq$methods(
	get.IPconf = function(ChIP_sample, Input_sample) {
		if(ChIP_sample %in% .l_confs) {
			if("conf" %in% names(.l_confs[[ChIP_sample]]) && is(.l_confs[[ChIP_sample]][["conf"]], "IPconf")) {
				conf <- .l_confs[[ChIP_sample]][["conf"]]
			} else if("file" %in% names(.l_confs[[ChIP_sample]])) {
				if(!file.exists(.l_confs[[ChIP_sample]])) 
					stop("IPconf file not exist!")
				conf <- readRDS(.l_confs[[ChIP_sample]][[file]])
				##set ChIP and Input
				ChIP <- .self$get.AlignedTags(ChIP_sample)
				Input <- .self$get.AlignedTags(Input_sample)
				conf$set.ChIP(ChIP)
				conf$set.Input(Input)
				.l_confs[[ChIP_sample]][["conf"]] <<- conf
			}
		}
		else {
		##conf not exist
			ChIP <- .self$get.AlignedTags(ChIP_sample)
			if(is.null(ChIP$bd_chrtcs))
				stop("Please run 'QC.AlignedTags' first!")
			Input <- .self$get.AlignedTags(Input_sample)
			if(is.null(Input$bd_chrtcs))
				stop("Please run 'QC.AlignedTags' first!")
			conf <- IPconf(ChIP=ChIP, Input=Input)
			.l_confs[[ChIP_sample]] <<- list(conf=conf)
		} 
		return(conf)
	}
)

##Assuming that ChIP and Input has already been cached 
ChIPSeq$methods(
	save.IPconf = function(ChIP_sample, Input_sample, file, ...) {
		if(missing(ChIP_sample))
			stop("Please specify 'ChIP_sample'!")
		if(missing(Input_sample))
			stop("Please specify 'Input_sample'!")
		if(! ChIP_sample %in% names(.l_confs))
			stop("IPconf object not exist!")
		if(! "conf" %in% names(.l_confs[[ChIP_sample]]))
			stop("IPconf object not exist!")
		if(! is(.l_confs[[ChIP_sample]][["conf"]], "IPconf"))
			stop("No IPconf object found!")
		.l_confs[[ChIP_sample]][["file"]] <<- file
		.l_confs[[ChIP_sample]][["conf"]]$save(file=file)
		.l_confs[[ChIP_sample]][["conf"]] <<- NULL
	}
)


##Pairwise IDR
ChIPSeq$methods(
	IDR.pairwise = function(conf1, conf2, ...) {
		##Get corresponding configuration
##		conf1 <- .self$get.IPconf(ChIP_sample=ChIP_rep1_sample, Input_sample=Input_sample)
##		conf2 <- .self$get.IPconf(ChIP_sample=ChIP_rep2_sample, Input_sample=Input_sample)
	}
)

##IDR for multiple configurations
ChIPSeq$methods(
	IDR = function() {

	}
)




















