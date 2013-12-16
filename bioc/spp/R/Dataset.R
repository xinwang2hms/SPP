##Class representing a set of ChIP-seq experiments

Dataset = setRefClass(
	Class = "Dataset", 
	fields = list(
		.sample_info = "df_Or_NULL", 
		.conf_info = "df_Or_NULL",
		.l_aligned_tags = "list_Or_NULL",
		.l_confs = "list_Or_NULL" 
	)
)

Dataset$methods(
	show = function(...) {
		cat("~~Dataset Object~~\n")
		cat("  total size: ", object.size.format(.self$size()), "\n", sep="")
		if(!is.null(.l_aligned_tags) && length(.l_aligned_tags) > 0) {
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
		} else
			atg.nload <- 0
		cat("  ", atg.nload, "/", nrow(.sample_info), " AlignedTags objects in memory\n", sep="")
		if(!is.null(.l_confs) && length(.l_confs) > 0) {
			conf.nload <- sum(unlist(sapply(1:length(.l_confs), function(x) {
				confID <- names(.l_confs)[x]
				if(! "conf" %in% names(.l_confs[[confID]]))
					return(FALSE)
				else if(! is(.l_confs[[confID]][["conf"]], "IPconf"))
					return(FALSE)
				else if(is(.l_confs[[confID]][["conf"]], "IPconf"))
					return(TRUE)
			})))
		} else 
			conf.nload <- 0
		cat("  ", conf.nload, " IPconf objects in memory\n", sep="")

	}
) 

Dataset$methods(
	size = function() {
		##size of AlignedTags objects
		if(!is.null(.l_aligned_tags)) {
			atg.size <- sum(unlist(lapply(.l_aligned_tags, function(x) {
				if(!is.null(x[["atg"]]) && is(x[["atg"]], "AlignedTags"))
					x[["atg"]]$size()
			})))
		} else atg.size <- 0
		##size of IPconf objects
		if(!is.null(.l_confs)) {
			conf.size <- sum(unlist(lapply(.l_confs, function(x) {
				if(!is.null(x[["conf"]]) && is(x[["conf"]], "IPconf"))
					x[["conf"]]$size()
			})))
		} else conf.size <- 0
		s <- object.size(.sample_info) + object.size(.conf_info) + conf.size + atg.size
		return(s)	
	}
)
##Initialization by existing object or providing 'sample_info' and 'conf_info'
Dataset$methods(
	initialize = function(..., sample_info, conf_info) {
		##if copy
		callSuper(...)
		##else new .sample_info and .l_aligned_tags
		if(is.null(.sample_info)) {
			.l_aligned_tags <<- list()
			if(!is.null(sample_info)) {
			if(!all(c("SampleID") %in% colnames(sample_info)))
				stop("Please provide 'SampleID'!")
			if(!any(c("TagFile", "AlignedTagsFile") %in% colnames(sample_info)))
				stop("Please provide 'TagFile' (BAM or Bowtie alignment result) and/or ", 
					"'AlignedTagsFile' (saved 'AlignedTags' objects in rds format)!")
			rownames(sample_info) <- as.character(sample_info[, "SampleID"])
			.sample_info <<- sample_info
			if("AlignedTagsFile" %in% colnames(sample_info)) {
				sapply(1:nrow(sample_info), function(i) {
					sample_i <- as.character(sample_info[i, "SampleID"])
					.l_aligned_tags[[sample_i]] <<- 
						list(file=as.character(sample_info[i, 
							"AlignedTagsFile"]), atg=NULL)
					if(!file.exists(as.character(sample_info[i, "AlignedTagsFile"])))
						stop("'AlignedTags' rds file does not exist for sample ", 
							sample_i, "!")
				})
			} 
##			else {
##				sapply(1:nrow(target), function(i) {
##					.l_aligned_tags[[as.character(target[i, "SampleID"])]] <<- 
##						list(file=NULL, atg=NULL)
##				})
##			}
			}
		}
		##new .conf_info and .l_confs
		if(is.null(.conf_info)) {
			.l_confs <<- list()
			if(!missing(conf_info)) {
			.conf_info <<- conf_info
			if(!all(c("confID", "ChIP", "Input", "IPconfFile") %in% colnames(conf_info)))
				stop("Please provide 'confID', 'ChIP', 'Input' and 'IPconfFile' columns!")
			if(sum(is.na(conf_info[, "ChIP"]))>0 || sum(is.na(conf_info[, "Input"]))>0)
				stop("Some ChIP and/or Input sample IDs are not provided in 'conf_info'!")
##			if("IPconfFile" %in% colnames(conf_info)) {
				sapply(1:nrow(conf_info), function(i) {
					conf_i <- paste(as.character(conf_info[i, "ChIP"]), 
						as.character(conf_info[i, "Input"], sep="_VS_"))
					.l_confs[[conf_i]] <<- list(ChIP=as.character(conf_info[i, "ChIP"]), 
						Input=as.character(conf_info[i, "Input"]), 
						file=as.character(conf_info[i, "IPconfFile"]))
					if(!file.exists(as.character(conf_info[i, "IPconfFile"])))
						stop("'IPconf' rds file does not exist for configuration ", 
							conf_i, "!")
				})
##			}
			}
		}
	}
)
Dataset$methods(
	update = function() {
		invisible(Dataset(.self))
	}
)

##Select samples

##Get AlignedTags object 
##3 senarios
##a: if object is saved in file, and already loaded in memory; then return reference to the object
##b: if object is saved in file, but not loaded in memory; then load the saved file and return reference
##c: if object is not available; then create this object. Arguments need to be specified in this scenario. 
Dataset$methods(
	get.AlignedTags = function(samples, type="BAM", genome_build="mm9", read_length=50,
		read_tag_names=F, read_tag_qualities=F, fix_chr_names=F, verbose=TRUE) {
		if(missing(samples) || !is.character(samples))
			stop("'samples' should be sample IDs!")
		if(! all(samples%in% as.character(.sample_info[, "SampleID"])))
			stop("Not all 'samples' exist!")
	ref <- NULL
	for(sample in samples) {
		if(! sample %in% names(.l_aligned_tags)) {
		##scenario c:
			.l_aligned_tags[[sample]] <<- list()
			if(!type %in% c("BAM", "Bowtie"))
				stop("Only 'BAM' and 'Bowtie' tags are supported now!")
			if(type=="BAM") {
				.l_aligned_tags[[sample]][["atg"]] <<- BAMTags(genome_build=genome_build, 
					read_length=read_length)
			} else if(type=="Bowtie") {
				.l_aligned_tags[[sample]][["atg"]] <<- BowtieTags(genome_build=genome_build, 
					read_length=read_length)
			}
			.l_aligned_tags[[sample]][["atg"]]$read(file=as.character(.sample_info[sample, "TagFile"]), 
				read_tag_names=read_tag_names, read_tag_qualities=read_tag_qualities, 
				fix_chr_names=fix_chr_names)
			##
			if(verbose)
				cat(ifelse(type=="BAM", "BAMTags", "BowtieTags"), "object created!\n")
		} else {
		##scenario b:
			if((! "atg" %in% names(.l_aligned_tags[[sample]])) || 
				is.null(.l_aligned_tags[[sample]][["atg"]])) {
				if(!"file" %in% names(.l_aligned_tags[[sample]]))
					stop("If 'TagFile' is not speficied in 'sample_info', 
						'AlignedTagsFile' should be!")
				if(!file.exists(.l_aligned_tags[[sample]][["file"]]))
					stop(paste("'AlignedTags' object does not exist in ", 
						.l_aligned_tags[[sample]][["file"]], sep=""))
				.l_aligned_tags[[sample]][["atg"]] <<- 
					tryCatch(readRDS(file.path(.l_aligned_tags[[sample]][["file"]])),
					error = function(e) print(e))
				if(verbose)
					cat("Object loaded for sample", sample, "from saved file!\n") 
			} else {
		##scenario a:
##				if(verbose)
##					cat("Object loaded for sample", sample, "!\n")
			}
		}
		ref <- c(ref, .l_aligned_tags[[sample]][["atg"]])
	}	
		names(ref) <- samples
		invisible(ref)
	}
)
Dataset$methods(
	load.AlignedTags = function(samples, type="BAM", genome_build="mm9", read_length=50, 
		read_tag_names=F, read_tag_qualities=F, fix_chr_names=F, verbose=TRUE) {
		if(missing(samples) || !is.character(samples))
			stop("'samples' should be a vector of sample IDs or 'all'!")
		allsamps <- as.character(.sample_info[, "SampleID"])
		allsamps.saved <- names(.l_aligned_tags)[which(unlist(lapply(.l_aligned_tags, function(x) !is.null(x[["file"]]))))]
##		allsamps.saved <- allsamps[which(unlist(lapply(.l_aligned_tags, function(x) !is.null(x[["file"]]))))]
		if(length(samples)==1 && samples=="all")
			samples <- allsamps.saved
		if(!all(samples%in%allsamps.saved))
			stop("Not all 'samples' were saved!")
		for(sample in samples) {
##			if((! "atg" %in% names(.l_aligned_tags[[sample]])) || 
##				is.null(.l_aligned_tags[[sample]][["atg"]])) {
				
##				if(!"file" %in% names(.l_aligned_tags[[sample]]))
##					stop("If 'TagFile' is not speficied in 'sample_info', 'AlignedTagsFile' should be!")
				if(!file.exists(.l_aligned_tags[[sample]][["file"]]))
					stop(paste("file ", file, " does not exist", sep=""))
				.l_aligned_tags[[sample]][["atg"]] <<- 
					tryCatch(readRDS(file.path(.l_aligned_tags[[sample]][["file"]])),
						error = function(e) print(e))
##			}
		}
		if(verbose)
			cat(length(samples), "object(s) loaded!\n") 
	}
)

Dataset$methods(
	loaded.AlignedTags = function() {
		allsamps.open <- names(.l_aligned_tags)[which(unlist(lapply(.l_aligned_tags, 
			function(x) {"atg" %in% names(x) && is(x[["atg"]], "AlignedTags")})))]
		invisible(allsamps.open)
	}
)

##Close AlignedTags object when not needed in order to save memory
##This method should be used when the user wants to just remove the object from memory 
##without making any change to saved file
Dataset$methods(
	close.AlignedTags = function(samples, verbose=TRUE) {
		if(missing(samples) || !is.character(samples))
			stop("'samples' should be a vector of sample IDs or 'all'!")
		allsamps.open <- .self$loaded.AlignedTags()
		if(length(samples)==1 && samples=="all")
			samples <- allsamps.open
		if(!all(samples%in%allsamps.open))
			stop("Not all 'samples' were loaded, check with 'loaded.AlignedTags'!")
		for(sample in samples) {
##			if(! sample %in% names(.l_aligned_tags) || 
##				is.null(.l_aligned_tags[[sample]][["atg"]]))
##				stop("'AlignedTags' object is not yet loaded in memory!")
			.l_aligned_tags[[sample]][["atg"]] <<- NULL
		}
		if(verbose)
			cat(length(samples), "'AlignedTags' object(s) closed!\n")
	}
)

##Save AlignedTags object
##This method should be used when the user wants to update/save AlignedTags object
##if files, should be full paths to files. 
##if dirname+prefix, each sample will be saved as dirname/prefix.sample.rds
##if both provided, will take files by default
Dataset$methods(
	save.AlignedTags = function(samples, files, dirname, prefix, verbose=TRUE) {
		##check 'samples'
		if(missing(samples) || !is.character(samples))
			stop("'samples' should be a vector of sample IDs or 'all'!")
		allsamps.open <- .self$loaded.AlignedTags()
		if(length(samples)==1 && samples=="all")
			samples <- allsamps.open
		if(!all(samples%in%allsamps.open))
			stop("Not all 'samples' were loaded, check with 'loaded.AlignedTags'!")
		##check 'files', 'dirname', 'prefix'
		if(missing(files) && (!missing(dirname))) {
			if(!file.exists(dirname))
				dir.create(dirname, showWarnings=TRUE, recursive=TRUE)
			if(!missing(prefix) && !is.null(prefix) && prefix!="") 
				files <- file.path(dirname, paste(prefix, samples, "rds", sep="."))
			else
				files <- file.path(dirname, paste(samples, "rds", sep="."))
		} else if(!missing(files)) {
			if(length(files)!=length(samples))
				stop("The No. of 'samples' does not match the No. of 'files'!")
		} else 
			stop("Please specify either 'files' or 'dirname' + 'prefix'(optional)!")
		for(s in 1:length(samples)) {
			sample <- samples[s]
			file <- files[s]
			.l_aligned_tags[[sample]][["file"]] <<- file
			ChIP <- .l_aligned_tags[[sample]][["atg"]]
			saveRDS(ChIP, file=.l_aligned_tags[[sample]][["file"]])
		}
		if(verbose)
			cat(length(samples), "'AlignedTags' objects saved for sample!\n")
	}
)

##QC for AlignedTags object in a sequencial manner
##Processed AlignedTags object is saved and closed to save memory
##a: compute binding characteristics
##b: qc
Dataset$methods(
	QC.AlignedTags = function(samples, verbose=TRUE, 
		param_atg=list(type="BAM", genome_build="mm9", read_length=50, 
			read_tag_names=F, read_tag_qualities=F, fix_chr_names=F), 
		param_cross_cor = list(srange=c(-100,500), bin=5, min_tag_count=1e3, 
			acceptance_z_score=3, accept_all_tags=FALSE), 
		param_QC_NRF=list(adjust=FALSE, size_adj_thresh=10e6, nsamp=100), 
		param_save=list(files=NULL, dirname=".", prefix="")) {

		##check 'samples'
		if(missing(samples) || !is.character(samples))
			stop("'samples' should be sample IDs or 'all'!")
		if(length(samples)==1 && samples=="all")
			samples <- as.character(.sample_info[, "SampleID"])
		##!check param_atg, param_cross_cor, param_QC_NRF
		.check.param <- function(.def, .arg) {
			unknown.args <- setdiff(names(.arg), names(.def))
			if(length(unknown.args)>0)
				warning(paste("Args Unkown: ", unknown.args, collapse=", "))
			def.args <- setdiff(names(.def), names(.arg))
			if(length(def.args)>0 && verbose) {
				cat("Args set by default: ", paste(paste(def.args, 
					unlist(.def[def.args]), sep="="), collapse=", "), sep="")
			}
##			c(.def[def.args], .arg[intersect(names(.def), names(.arg))])	
		}
		.param.atg.def <- list(type="BAM", genome_build="mm9", read_length=50, 
			read_tag_names=F, read_tag_qualities=F, fix_chr_names=F)
		.param.cross.cor <- list(srange=c(-100,500), bin=5, min_tag_count=1e3, 
			acceptance_z_score=3, accept_all_tags=FALSE)
		.param.QC.NRF <- list(adjust=FALSE, size_adj_thresh=10e6, nsamp=100)
		.param.save <- list(files=NULL, dirname=".", prefix="")
##		.check.param(.param.atg.def, param_atg)
##		.check.param(.param.cross.cor, param_cross_cor)
##		.check.param(.param.QC.NRF, param_QC_NRF)
##		.check.param(.param.save, param_save)
		##

		qc.summary <- list(SampleID=samples, NRF=NULL, NRF_nostrand=NULL, 
			NRF_LibSizeadjusted=NULL, NSC=NULL, RSC=NULL, QFlag=NULL)
		for(sample in samples) {
			atg.ref <- do.call(.self$get.AlignedTags, c(list(sample=sample), param_atg))
			atg.ref <- atg.ref[[1]]
			do.call(atg.ref$compute.cross.cor, param_cross_cor)
			do.call(atg.ref$NRF, param_QC_NRF)
			atg.ref$phantomPeak()
			do.call(.self$save.AlignedTags, c(list(samples=sample), param_save))
			##
			qc.summary$NRF <- c(qc.summary$NRF, atg.ref$qc$NRF$NRF)
			qc.summary$NRF_nostrand <- c(qc.summary$NRF_nostrand, atg.ref$qc$NRF$NRF_nostrand)
			qc.summary$NRF_LibSizeadjusted <- c(qc.summary$NRF_LibSizeadjusted, atg.ref$qc$NRF$NRF_LibSizeadjusted)
			qc.summary$NSC <- c(qc.summary$NSC, atg.ref$qc$phantom_peak$NSC)
			qc.summary$RSC <- c(qc.summary$RSC, atg.ref$qc$phantom_peak$RSC)
			qc.summary$QFlag <- c(qc.summary$QFlag, atg.ref$qc$phantom_peak$quality_flag)
			##
			.self$close.AlignedTags(samples=sample)
			atg.ref <- NULL
		}
		##QC summary
		qc.report <- data.frame(SampleID=qc.summary$SampleID, 
			NRF=qc.summary$NRF, NRF_nostrand=qc.summary$NRF_nostrand, 
			NRF_adj= qc.summary$NRF_LibSizeadjusted, NSC=qc.summary$NSC, 
			RSC=qc.summary$RSC, QFlag=qc.summary$QFlag)
		if(nrow(qc.report)>5)
			print(head(qc.report))
		else 
			print(qc.report)
		invisible(qc.report)
	} 
)


Dataset$methods(
	.check.sample.pairing = function(ChIP_samples, Input_samples, sample_pairs) {
		
		if((missing(ChIP_samples) || missing(Input_samples)) && missing(sample_pairs))
			stop("Either ('ChIP_samples' + 'Input_samples') or 'sample_pairs' should be specified!")
		##if ChIP_samples and Input_samples are specified, sample_pairs will be ignored
		if(!missing(ChIP_samples) && !missing(Input_samples)) {
			if(!is(ChIP_samples, "character") || !is(Input_samples, "character"))
				stop("'ChIP_samples' and 'Input_samples' should be SampleIDs!")
			if(length(ChIP_samples)==0 || length(Input_samples)==0)
				stop("'ChIP_samples' and 'Input_samples' should be SampleIDs!")
			##check pairing
			##
			if(length(ChIP_samples)==1 && length(Input_samples)==1) cat("One ChIP, One Input\n")
			else if(length(ChIP_samples)>1 && length(Input_samples)==1) {
				cat("Multiple ChIPs, One Input\n")
				Input_samples <- rep(Input_samples, length(ChIP_samples))
			} else if (length(ChIP_samples)>1 && length(Input_samples)>1) {
				cat("Multiple ChIPs, Multiple Inputs")
				if(length(ChIP_samples)!=length(Input_samples))
					stop("No. of ChIP samples does not match No. of Input samples!")
			}
			return(list(ChIP_samples=ChIP_samples, Input_samples=Input_samples))
		} else if(!missing(sample_pairs)) {
			if(!is(sample_pairs, "data.frame"))
				stop("'sample_pairs' should be a data frame!")
			if(!all(c("ChIP", "Input") %in% colnames(sample_pairs)) || ncol(sample_pairs)!=2)
				stop(paste("'sample_pairs' should only contain columns of ChIP and", 
					" Input sample IDs, and column names should be 'ChIP' and 'Input'!", 
					sep=""))
			if(nrow(sample_pairs)==0)
				stop("No sample pairs found in 'sample_pairs'!")
			return(list(ChIP_samples=as.character(sample_pairs[, "ChIP"]), 
				Input_samples=as.character(sample_pairs[, "Input"])))
		}
	}
)
##load IPconf object(s)
Dataset$methods(
	load.IPconf = function(confIDs, verbose=TRUE, load_atg=FALSE) {
		if(missing(confIDs) || !is.character(confIDs) || length(confIDs)<1)
			stop("'confIDs' should be a character vector of IPconf IDs!")
		if(is.null(.l_confs))
			stop("No saved IPconf files found!")
		
		allconfs.saved <- names(.l_confs)[which(unlist(lapply(.l_confs, 
			function(x) !is.null(x[["file"]]))))]
		if(length(confIDs)==1 && confIDs == "all")
			confIDs <- allconfs.saved
		if(!all(confIDs %in% names(.l_confs)))
			stop("Not all IPconf objects saved in files!")
		for(confID in confIDs) {
			if(!file.exists(.l_confs[[confID]][["file"]]))
				stop(paste("file ", file, " does not exist", sep=""))
			.l_confs[[confID]][["conf"]] <<- 
				tryCatch(readRDS(file.path(.l_confs[[confID]][["conf"]])),
					error = function(e) print(e))
		}
		if(verbose)
			cat(length(confIDs), " IPconf object(s) loaded!\n")
	}
)

##get IPconf object(s)
##a: if already in memory; then return reference
##b: if saved in file; then load to R and return reference; if conf is saved in file, 
##   filenames of corresponding ChIP and Input should also be provided 
##c: if not in memory, nor in file; then stop, call get.AlignedTags first to 
##   get ChIP and Input; 

##Input formats:
##(1) multiple ChIPs, one Input (including one ChIP, one Input)
##(2) multiple paired ChIP-Input
##(3) specify by sample_pairs
Dataset$methods(
	get.IPconf = function(ChIP_samples, Input_samples, sample_pairs, verbose=TRUE, link_atg=TRUE) {
		sampp <- .self$.check.sample.pairing(ChIP_samples, Input_samples, sample_pairs)
		ChIP_samples <- sampp$ChIP_samples
		Input_samples <- sampp$Input_samples
		
		confIDs <- paste(ChIP_samples, Input_samples, sep="_VS_")
		ref <- NULL
		for(s in 1:length(confIDs)) {
			confID <- confIDs[s]
			ChIP_sample <- ChIP_samples[s]
			Input_sample <- Input_samples[s]
		##a:
		if(confID %in% names(.l_confs)) {
			if("conf" %in% names(.l_confs[[confID]]) && 
				is(.l_confs[[confID]][["conf"]], "IPconf")) {

				conf <- .l_confs[[confID]][["conf"]]
			} else if("file" %in% names(.l_confs[[confID]])) {
		##b:
				if(!file.exists(.l_confs[[confID]][["file"]])) 
					stop("IPconf file does not exist!")
				conf <- readRDS(.l_confs[[confID]][["file"]])
				if(! confID%in%names(.l_confs))
					.l_confs[[confID]] <<- list()
				.l_confs[[confID]][["conf"]] <<- conf
				if(verbose)
					cat("'IPconf' object loaded from file ", .l_confs[[confID]][["file"]], 
						" for ChIP sample ", ChIP_sample, " and Input sample ", Input_sample, "\n")
			}
			if(link_atg) {
##				if(!file.exists(.l_aligned_tags[[ChIP_sample]][["file"]]))
##				stop("filename for ", ChIP_sample, " is not available in 'sample_info'!")
##				if(!file.exists(.l_aligned_tags[[Input_sample]][["file"]]))
##				stop("filename for ", Input_sample, " is not available in 'sample_info'!")

				if(! ChIP_sample %in% names(.l_aligned_tags) || 
					(is.null(.l_aligned_tags[[ChIP_sample]][["atg"]]) && 
					 is.null(.l_aligned_tags[[ChIP_sample]][["file"]])))
					stop(paste(ChIP_sample, " 'AlignedTags' object is neither loaded in memory nor saved in file!"), sep="")
				if(! Input_sample %in% names(.l_aligned_tags) || 
					(is.null(.l_aligned_tags[[Input_sample]][["atg"]]) && 
					 is.null(.l_aligned_tags[[Input_sample]][["file"]])))
					stop(paste(Input_sample, " 'AlignedTags' object is neither loaded in memory nor saved in file!"), sep="")
				ChIP <- .self$get.AlignedTags(ChIP_sample)[[1]]
				Input <- .self$get.AlignedTags(Input_sample)[[1]]
				conf$set.ChIP(ChIP)
				conf$set.Input(Input)
			}
		} else {
		##c:
			if(! ChIP_sample %in% names(.l_aligned_tags) || 
				(is.null(.l_aligned_tags[[ChIP_sample]][["atg"]]) && 
				 is.null(.l_aligned_tags[[ChIP_sample]][["file"]])))
				stop(paste(ChIP_sample, " 'AlignedTags' object is neither loaded in memory nor saved in file!"), sep="")
			if(! Input_sample %in% names(.l_aligned_tags) || 
				(is.null(.l_aligned_tags[[Input_sample]][["atg"]]) && 
				 is.null(.l_aligned_tags[[Input_sample]][["file"]])))
				stop(paste(Input_sample, " 'AlignedTags' object is neither loaded in memory nor saved in file!"), sep="")
			ChIP <- .self$get.AlignedTags(ChIP_sample)[[1]]
			Input <- .self$get.AlignedTags(Input_sample)[[1]]
			conf <- IPconf(ChIP=ChIP, Input=Input)
			if(! confID%in%names(.l_confs))
				.l_confs[[confID]] <<- list()
			.l_confs[[confID]][["conf"]] <<- conf
			if(verbose)
				cat("'IPconf' object created for ChIP sample ", ChIP_sample, " and Input sample ", 
					Input_sample, "!\n")
		} 
		.l_confs[[confID]][["ChIP"]] <<- ChIP_sample
		.l_confs[[confID]][["Input"]] <<- Input_sample
		ref <- c(ref, conf)
		}
		names(ref) <- confIDs
		invisible(ref)
	}
)

##Close IPconf without saving
Dataset$methods(
	close.IPconf = function(ChIP_samples, Input_samples, sample_pairs, verbose=TRUE, ...) {
		sampp <- .self$.check.sample.pairing(ChIP_samples, Input_samples, sample_pairs)
		ChIP_samples <- sampp$ChIP_samples
		Input_samples <- sampp$Input_samples
		confIDs <- paste(ChIP_samples, Input_samples, sep="_VS_")
		for(s in 1:length(confIDs)) {
			confID <- confIDs[s]
			ChIP_sample <- ChIP_samples[s]
			Input_sample <- Input_samples[s]
##			if(missing(ChIP_sample))
##				stop("Please specify 'ChIP_sample'!")
##			if(missing(Input_sample))
##				stop("Please specify 'Input_sample'!")
##			confID <- paste(ChIP_sample, Input_sample, sep="_VS_")
			if(! confID %in% names(.l_confs) )
				stop("IPconf object does not exist!")
			if(! ("conf" %in% names(.l_confs[[confID]]) && 
				is(.l_confs[[confID]][["conf"]], "IPconf")))
				stop("IPconf object does not exist!")
			.l_confs[[confID]][["conf"]] <<- NULL
		}
		if(verbose)
			cat(length(confIDs), "IPconf objects closed!\n")
		
	}
)

Dataset$methods(
	loaded.IPconf = function() {
		allconf.open <- names(.l_confs)[which(unlist(lapply(.l_confs, 
			function(x) {"conf" %in% names(x) && is(x[["conf"]], "IPconf")})))]
		if(length(allconf.open) > 0) {
			ChIPs <- unlist(sapply(allconf.open, function(x) .l_confs[[x]][["ChIP"]]))
			Inputs <- unlist(sapply(allconf.open, function(x) .l_confs[[x]][["Input"]]))
			allconf.open <- data.frame(confID=allconf.open, ChIP=ChIPs, Input=Inputs)
		}
		invisible(allconf.open)
	}
)
##Assuming that ChIP and Input has already been cached 
Dataset$methods(
	save.IPconf = function(ChIP_samples, Input_samples, sample_pairs, files, dirname, prefix, verbose=TRUE, ...) {
		sampp <- .self$.check.sample.pairing(ChIP_samples, Input_samples, sample_pairs)
		ChIP_samples <- sampp$ChIP_samples
		Input_samples <- sampp$Input_samples
		confIDs <- paste(ChIP_samples, Input_samples, sep="_VS_")
		##check 'files', 'dirname', 'prefix'
		if(missing(files) && (!missing(dirname))) {
			if(!file.exists(dirname))
				dir.create(dirname, showWarnings=TRUE, recursive=TRUE)
			if(!missing(prefix) && !is.null(prefix) && prefix!="")
				files <- file.path(dirname, paste(prefix, confIDs, "rds", sep="."))
			else
				files <- file.path(dirname, paste(confIDs, "rds", sep="."))
		} else if(!missing(files)) {
			if(length(files)!=length(confIDs))
				stop("The No. of 'IPconf' objects does not match the No. of 'files'!")
		} else
			stop("Please specify either 'files' or 'dirname' + 'prefix'(optional)!")
		for(s in 1:length(confIDs)) {
			confID <- confIDs[s]
			ChIP_sample <- ChIP_samples[s]
			Input_sample <- Input_samples[s]
			file <- files[s]
##		if(missing(ChIP_sample))
##			stop("Please specify 'ChIP_sample'!")
##		if(missing(Input_sample))
##			stop("Please specify 'Input_sample'!")
##		confID <- paste(ChIP_sample, Input_sample, sep="_VS_")
			if(! confID %in% names(.l_confs) )
				stop("IPconf object does not exist!")
			if(! ("conf" %in% names(.l_confs[[confID]]) && 
				is(.l_confs[[confID]][["conf"]], "IPconf")))
				stop("IPconf object does not exist!")
			.l_confs[[confID]][["file"]] <<- file
			.l_confs[[confID]][["conf"]]$save(file=file)
		}
		if(verbose)
			cat(length(confIDs), "IPconf objects saved!\n")
	}
)

##sample pairing
##samples: (1st filtering) used to pre-filter samples
##filter_def: (2nd filtering) used to filter by a list of factors 
##cond_def (optional): definition of conditions. list of factors to stratify samples to condition groups
##ChIP_def: definition of ChIP. used to look for ChIP samples under each condition
##Input_def: definition of Input. used to look for Input samples under each condition


Dataset$methods(
	sample.pairing = function(samples="all", filter_def=list(ChromatinBatch="20100608", Sex="F"), 
		cond_def=c("Age"), ChIP_def=list(IP=c("h3k9me3")), ChIP_samples, Input_global, 
		Input_def=list(Diet="Normal", IP="input")) {
		##check 'sample'
		if(!is.character(samples))
			stop("'samples' should be a character vector of sample IDs or 'all'!")
		allsamps <- as.character(.sample_info[, "SampleID"])
		if(length(samples)==1 && samples=="all")
			samples <- allsamps
		if(!all(samples%in%allsamps))
			stop("Not all samples are found in sample Info!")
		##check def types
		if(!missing(filter_def) && !is.list(filter_def))
			stop("'filter_def' should be a list!")
		if(!missing(cond_def) && !is.character(cond_def))
			stop("'cond_def' should be a character of vector specifying factors defining conditions!")
		if(!missing(ChIP_def) && !is.list(ChIP_def))
			stop("'ChIP_def' should be a list!")
		if(!missing(Input_def) && !is.list(Input_def))
			stop("'Input_def' should be a list!")
		##check Input_def and Input_global
		if(missing(Input_def) && (missing(Input_global) || length(Input_global)!=1))
			stop("Either set a global Input or define Input in Input_def!")
		##check ChIP_def and ChIP_samples
		if(missing(ChIP_def) && (missing(ChIP_samples)))
			stop("Either set a set of ChIP samples or define ChIP in Input_def!")
		if(!missing(ChIP_samples) && !missing(cond_def))
			stop("ChIP_samples and cond_def should not be specified at the same time!")
		##!check def existance

		##step 1: 1st filtering by samples
		samp.df <- .sample_info[samples, ]		
		##step 2: 2nd filtering by filter_def
		if(!missing(filter_def)) {
			for(f in 1:length(filter_def)) {
				samp.df <- samp.df[which(samp.df[, names(filter_def)[f]] %in% filter_def[[f]]), , drop=FALSE]
			}
		}
		if(nrow(samp.df) < 2)
			stop("Less than 2 samples left after filtering!")
		##step 3: split conditions
		if(!missing(cond_def)) {
			samp.df.sp <- split(samp.df, lapply(cond_def, function(x) samp.df[, x]))
		} else {
			samp.df.sp <- list(samp.df)
		}
		##step 4: for each condition, look for ChIP and Input samples to pair
		l.sample.pairs <- list()
		for(cond in 1:length(samp.df.sp)) {
			##look for Input, currently accept only one Input
			if(missing(Input_def))
				Input <- Input_global
			else {
				tmp.df <- samp.df.sp[[cond]]
				for(f in 1:length(Input_def)) {
					tmp.df <- tmp.df[which(tmp.df[, names(Input_def)[f]]%in%Input_def[[f]]), , drop=FALSE]
				}
				if(nrow(tmp.df)!=1)
					stop(paste("Multiple Input samples defined by 'Input_def':", 
						paste(tmp.df[, "SampleID"], collapse=', '), sep=""))
				Input <- as.character(tmp.df[, "SampleID"])
			}
			##look for ChIP
			if(!missing(ChIP_samples))
				ChIP <- ChIP_samples
			else {
				tmp.df <- samp.df.sp[[cond]]
				for(f in 1:length(ChIP_def)) {
					tmp.df <- tmp.df[which(tmp.df[, names(ChIP_def)[f]]%in%ChIP_def[[f]]), , drop=FALSE]
				}
				ChIP <- as.character(tmp.df[, "SampleID"])
			}
			Input <- rep(Input, length(ChIP))
			l.sample.pairs[[cond]] <- data.frame(ChIP=ChIP, Input=Input)
		}
		names(l.sample.pairs) <- names(samp.df.sp)
		l.sample.pairs
	}
)

##Pairwise IDR
Dataset$methods(
	IDR.pairwise = function(conf1, conf2, ...) {
		##Get corresponding configuration
##		conf1 <- .self$get.IPconf(ChIP_sample=ChIP_rep1_sample, Input_sample=Input_sample)
##		conf2 <- .self$get.IPconf(ChIP_sample=ChIP_rep2_sample, Input_sample=Input_sample)
	}
)



##IDR for multiple configurations
Dataset$methods(
	IDR = function(ChIP_samples, Input_samples, sample_pairs, ...) {

	}
)




















