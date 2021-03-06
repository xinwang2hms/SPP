##-1. AlignedTags 
##-The class is designed for representing aligned ChIP-seq tags
##-note: this is a root class, which is inherited by specific classes
setClassUnion("smoothedTagDensity_Or_NULL", c("smoothedTagDensity", "NULL"))
AlignedTags = setRefClass(
	Class = "AlignedTags",
	fields = list(
		file = "char_Or_NULL", 				##read from file later	
		genome_build = "char_Or_NULL",			##e.g. mm9, hg19
		read_length = "numeric_Or_NULL", 		##read length
		tags = "list_Or_NULL", 				##sequence tags
		quality = "list_Or_NULL",			##not in use actually 
		names = "list_Or_NULL", 			
		bd_chrtcs = "list_Or_NULL", 			##binding characteristics
		qc = "list_Or_NULL",				##quality control
		smoothed_density = "smoothedTagDensity_Or_NULL"
	)
)
setClassUnion("AlignedTags_Or_NULL", c("AlignedTags", "NULL"))
AlignedTags$methods(
	initialize = function(..., genome_build=NULL, read_length=NULL) {
		callSuper(...)
		if(is.null(.self$genome_build))
			genome_build <<- genome_build
		if(is.null(.self$read_length))
			read_length <<- read_length
		if(is.null(smoothed_density))
			smoothed_density <<- smoothedTagDensity(ChIP=.self)
		else {
			tmp <- .self$smoothed_density$.profile
			.self$smoothed_density <<- smoothedTagDensity(ChIP=.self, param=.self$smoothed_density$.param)
			.self$smoothed_density$.profile <<- tmp
		}
	}
)
##object size
AlignedTags$methods(
	size = function() {
		s <- object.size(file) + object.size(genome_build) + 
			object.size(read_length) + object.size(tags) + 
			object.size(quality) + object.size(names) + 
			object.size(bd_chrtcs) + smoothed_density$size()
		return(s)
	}
)


##1. display a brief summary of the object
AlignedTags$methods(
	show = function(...) {
		nchrs <- length(tags)
		ntags <- sum(unlist(lapply(tags, length)))
		cat("~~", as.character(class(.self)), "object~~\n")
		cat("  total size: ", object.size.format(.self$size()), "\n", sep="")
		if(!is.null(tags)) {
			cat(paste("  ", ntags, " fragments", " across ", nchrs, 
				" chromosome(s)", "\n", sep=""))			
			cat(paste("  read from file '", file, "'", "\n", sep=""))
		}
		##show binding characteristics
		if(!is.null(bd_chrtcs)) {
			##
			cat("Binding characteristics:\n")
			cat(paste("  cross correlation peak: Position=", 
				bd_chrtcs$peak$x, ", Height=", 
				round(bd_chrtcs$peak$y, 3), "\n", sep=""))
			cat(paste("  optimized window half-size: ", 
				bd_chrtcs$whs, "\n", sep=""))
		}	
		if(!is.null(qc)) {
			cat("Quality control:\n")
			if(!is.null(qc$phantom_peak))
				cat("  NSC=", round(qc$phantom_peak$NSC, 2), ", RSC=", 
					round(qc$phantom_peak$RSC, 2), ", Quality flag: ", 
					qc$phantom_peak$quality_flag, "\n", sep="")
			if(!is.null(qc$NRF))
				cat("  NRF=", round(qc$NRF$NRF, 2), ", NRF (no strand)=", 
					round(qc$NRF$NRF_nostrand, 2), ", NRF (adjusted)=", 
					round(qc$NRF$NRF_LibSizeadjusted, 2), "\n", sep="")
		}
	}
)

##2. get a subset 
AlignedTags$methods(
	subset = function(filter, ...) {
		##check arguments
		CHRs <- names(tags)
		if(is(filter, "character")) {
			chrs <- filter
			if(!(all(chrs %in% CHRs)))
				stop(paste("Some chromosome(s) are not included in the data!\n", 
				"Here are the chromosomes included:\n", 
				paste(CHRs, collapse=", "), sep=""))
			##filtering
			tags[setdiff(CHRs, chrs)] <<- NULL
			quality[setdiff(CHRs, chrs)] <<- NULL
			if(length(names)>0)
				names[setdiff(CHRs, chrs)] <<- NULL			
		} else
			stop("'filter' should either be a vector of chromosomes")
	}
)

##sampling
AlignedTags$methods(
	sampling = function(f, replace=FALSE) {
		if(!is.numeric(f)) 
			stop("f should be numeric!")
		if(f<0)
			stop("f should be > 0!")
		if(is.null(tags))
			stop("Tags not read yet!")

		trash <- sapply(1:length(tags), function(x) {
			inds <- sample(1:length(tags[[x]]), round(f*length(tags[[x]])), 
				replace=replace)
			inds <- inds[order(inds, decreasing=FALSE)]
			tags[[x]] <<- tags[[x]][inds]
			if(!is.null(quality))
				quality[[x]] <<- quality[[x]][inds]
			if(!is.null(names)) {
				if(replace)
					names[[x]] <<- NULL
				else
					names[[x]] <<- names[[x]][inds]
			}
		})
	}
)



##3. remove tag anomalies
AlignedTags$methods(
	remove.tag.anomalies = function(bin=1, trim_fraction=1e-3, z=5, 
		zo=3*z, var_base=0.1, ...) {
		
		t.remove.tag.anomalies <- function(tv, bin=1, trim_fraction=1e-3, 
			z=5, zo=3*z, return_indecies=F) {
			
			tt <- table(floor(tv/bin))
			# trim value
			stt <- sort(as.numeric(tt))
			stt <- stt[1:(length(stt)*(1-trim_fraction))]
			mtc <- mean(stt)
			tcd <- sqrt(var(stt)+var_base)

			thr <- max(1,ceiling(mtc+z*tcd))
			thr.o <- max(1,ceiling(mtc+zo*tcd))
			# filter tt
			tt <- tt[tt>thr]
			# get + and - tags
			tp <- as.numeric(names(tt))
			pti <- tp>0;
			it <- intersect(tp[pti],(-1)*tp[!pti])
			# add one-strand matches
			it <- unique(c(it,tp[tt>thr.o]))
			sit <- c(it,(-1)*it)

			if(bin>1) {
			  sit <- sit*bin;
			  sit <- c(sit,unlist(lapply(1:bin,function(i) sit+i)))
			}
			if(return_indecies) {
			  return(!tv %in% sit)
			} else {
			  return(tv[!tv %in% sit])
			}
		}
		
		vil <- lapply(tags, t.remove.tag.anomalies, 
			return_indecies=T, bin=bin, trim_fraction=trim_fraction, 
			z=z, zo=zo)		
		chrl <- names(tags)
		names(chrl) <- chrl
		tags <<- lapply(chrl, function(chr) 
			tags[[chr]][vil[[chr]]])
		# count tags to remove empty chromosomes
		nt <- unlist(lapply(tags, length))		
		if(any(nt==0)) 
			tags <<- tags[nt!=0]		
		if(!is.null(quality)) {
			quality <<- lapply(chrl, function(chr) 
				quality[[chr]][vil[[chr]]])
			quality <<- quality[nt!=0]
		}
		if(!is.null(names)) {
			names <<- lapply(chrl, function(chr) 
				names[[chr]][vil[[chr]]])
			names <<- names[nt!=0]
		}			
	}
)	

##4. compute binding characteristics
AlignedTags$methods(
	compute.cross.cor = function(srange=c(-100,500), bin=5, 
		min_tag_count=1e3, acceptance_z_score=3, 
		accept_all_tags=FALSE, ...) {
		
		
		# take highest quality tag bin
		if(!is.null(quality) && !accept_all_tags) {
			min.bin <- min(unlist(lapply(quality, min)))
			chrl <- names(tags)
			names(chrl) <- chrl
			otl <- lapply(chrl, function(chr) 
				tags[[chr]][quality[[chr]]==min.bin])
		} else 
			otl <- tags
		# remove empty chromosomes
		otl <- otl[unlist(lapply(otl,length))!=0]
		#check if parallel
		spp.cores <- getOption("spp.cores")
		if(!is.null(spp.cores) && spp.cores>1 
			&& "package:multicore" %in% search() && length(otl) > 1) {
			
			cc <- mclapply(otl, tag.scc, srange=srange, bin=bin, 
				mc.cores=spp.cores, mc.preschedule=F)
		} else {
			cc <- lapply(otl,tag.scc,srange=srange,bin=bin)
		}
		# calculate strand scc
##		if(!is.null(cluster)) {
##			cc <- clusterApplyLB(cluster,otl,tag.scc,srange=srange,bin=bin)
##			names(cc) <- names(otl)
##		} else {
##			cc <- lapply(otl,tag.scc,srange=srange,bin=bin)
##		}
				
		ccl<-list(sample=cc)
		ccl.av <- lapply(names(ccl),t.plotavcc,type='l',ccl=ccl,return.ac=T,ttl=list(sample=otl),plot=F)[[1]]
		ccl.av <- data.frame(x=as.numeric(names(ccl.av)),y=as.numeric(ccl.av))

		# find peak
		pi <- which.max(ccl.av$y)

		# determine width at third-height
		th <- (ccl.av$y[pi]-ccl.av$y[length(ccl.av$y)])/3+ccl.av$y[length(ccl.av$y)]
		whs <- max(ccl.av$x[ccl.av$y>=th])

		# determine acceptance of different quality bins

		# calculates tag scc for the best tags, and combinations of best tag category with every other category
		# for subsequent selection of acceptable categories
		scc.acceptance.calc <- function() {

			qr <- range(unlist(lapply(quality,range)))

			# start with best tags

			# determine half-width for scc calculations
			pi <- which.max(ccl.av$y)

			# determine width at half-height
			th <- (ccl.av$y[pi]-ccl.av$y[length(ccl.av$y)])/2+ccl.av$y[length(ccl.av$y)]
			lwhs <- max(ccl.av$x[ccl.av$y>=th])-ccl.av$x[pi]
			lwhs <- max(c(20,bin*10,lwhs))
			srange <- ccl.av$x[pi]+c(-lwhs,lwhs)
			# calculate chromosome-average scc		
			t.scc <- function(tags) {
				if(!is.null(spp.cores) && spp.cores>1 
					&& "package:multicore" %in% search() && length(tags) > 1) {
					
					cc <- mclapply(tags, tag.scc, srange=srange, 
						bin=bin, mc.cores=spp.cores, mc.preschedule=F)
				} else {
					cc <- lapply(tags,tag.scc,srange=srange,bin=bin)
				}
				return(t.plotavcc(1,type='l',ccl=list(cc),ttl=list(tags),plot=F,return.ac=T))
			}
			
##			t.scc <- function(tags) {
##			  if(is.null(cluster)) {
##				cc <- lapply(tags,tag.scc,srange=srange,bin=bin)
##			  } else {
##				cc <- clusterApplyLB(cluster,tags,tag.scc,srange=srange,bin=bin)
##				names(cc) <- names(tags)
##			  }		  
##			  return(t.plotavcc(1,type='l',ccl=list(cc),ttl=list(tags),plot=F,return.ac=T))
##			}

			# returns info list for a given tag length (lv), mismatch count (nv)
			t.cat <- function(qual) {
			  # construct tag set
			  if(qual==qr[1]) {
				ts <- otl
			  } else {
				nts <- names(otl)
				names(nts) <- nts
				# select tags
				at <- lapply(nts,function(chr) 
					tags[[chr]][quality[[chr]]==qual])
				ntags <- sum(unlist(lapply(at,length)))
				if(ntags<min_tag_count) 
					return(NULL)

				# append to otl
				ts <- lapply(nts,function(nam) 
					c(otl[[nam]],at[[nam]]))
			  }
			  return(t.scc(ts))
			}

			# calculate cross-correlation values for each quality bin
			ql <- sort(unique(unlist(lapply(quality,unique))))
			names(ql) <- ql

			qccl <- lapply(ql,t.cat)

			# acceptance tests
			ac <- c(T,unlist(lapply(qccl[-1],function(d) 
				if(is.null(d)) { 
					return(F) 
				} else { 
					t.test(d-qccl[[as.character(min.bin)]], 
						alternative="greater")$p.value < 
						pnorm(acceptance_z_score,lower.tail=F) 
				})))
			names(ac) <- names(qccl)
			return(list(informative_bins=ac, quality_cc=qccl))
		}

##		if(accept_all_tags | is.null(object@quality)) {
##			return(list(cross.correlation=ccl.av,peak=list(x=ccl.av$x[pi],y=ccl.av$y[pi]),whs=whs))    
##		} else {
##			acc <- scc.acceptance.calc();
##			return(list(cross.correlation=ccl.av,peak=list(x=ccl.av$x[pi],y=ccl.av$y[pi]),whs=whs,quality.bin.acceptance=acc))
##		}			
		if(is.null(bd_chrtcs)) 
			bd_chrtcs <<- list()
		bd_chrtcs$param <<- list(srange=srange, bin=bin, min_tag_count=min_tag_count, 
			acceptance_z_score=acceptance_z_score, accept_all_tags=accept_all_tags)
		bd_chrtcs$cross_cor <<- ccl.av
		bd_chrtcs$peak <<- list(x=ccl.av$x[pi],y=ccl.av$y[pi])
		bd_chrtcs$whs <<- whs
		if(!(accept_all_tags | is.null(quality))) 
			bd_chrtcs$quality_bin_accpt <<- scc.acceptance.calc()
	}
)

##plot cross correlation 
AlignedTags$methods(
	view.cross.cor = function() {
		if(is.null(bd_chrtcs))
			stop("Please run 'compute.cross.cor' first!")
		plot(bd_chrtcs$cross_cor, type='l', xlab="strand shift", 
			ylab="cross-correlation")
		abline(v=bd_chrtcs$peak$x, lty=2, lwd=3, col='red')
		if(!is.null(qc) && !is.null(qc$phantom_peak)) {
			abline(v=qc$phantom_peak$phantom_cc$x, lty=2, lwd=1, col='blue')
		}	
	}
)

##remove local tag anomalies
AlignedTags$methods(
	remove.local.tag.anomalies = function(window_size=200, 
		eliminate_fold=10, cap_fold=4, z_threshold=3, ...) {
	
		tags <<- lapply(tags, filter.singular.positions.by.local.density, 
			window_size=window_size, eliminate_fold=eliminate_fold, 
			cap_fold=cap_fold, z_threshold=z_threshold)
				
	}
)

##retrieve binding characteristics
AlignedTags$methods(
	get.cross.cor = function() {
		return(bd_chrtcs)
	}
)
##set binding characteristics
AlignedTags$methods(
	set.cross.cor = function(value) {
		bd_chrtcs <<- value
	}
)


## Get NRF scores (Non Redundant Fraction)
## if resampling==TRUE, size_adj_thresh and nsamp are used; otherwise not. 
AlignedTags$methods(
        NRF = function(adjust=FALSE, size_adj_thresh=10e6, nsamp=100) {
                # total number of tags
                ALL_TAGS<-sum(sapply(tags, length))

                # total number of unique positions (with strand specificity)
                UNIQUE_TAGS<-sum(sapply(lapply(tags, unique), length))

                # total number of unique positions (without strand specificity)
                UNIQUE_TAGS_nostrand<-sum(sapply(lapply(tags, FUN=function(x) {unique(abs(x))}), length))

                # Non Redundant Fraction
                NRF<-UNIQUE_TAGS/ALL_TAGS
                # Non Redundant Fraction without strand specificity
                NRF_nostrand<-UNIQUE_TAGS_nostrand/ALL_TAGS


                # With very large libsizes the non redundant fraction might decrease due to
                # the sequencing depth being extremely high rather than the library complxity being low
                ## to compensate for lib size differences we try recomputing the NRF with a subset of 10million reads

                # handle the taglist as a vector instead than as a list for uniform sampling across cheomosomes
##		nomi<-rep(names(tags), sapply(tags, length))
##		chip.data<-unlist(tags)
##		names(chip.data)<-NULL
                # use chsomosome names + reads positions (strand specific) for counting unique tags
##		chip.data<-paste(nomi, chip.data, sep="")
		NRF_LibSizeadjusted <- NA
	if(adjust) {
		chr_n_tags <- unlist(lapply(tags, length))
		chr_samp_tags <- round(size_adj_thresh*chr_n_tags/sum(chr_n_tags))

		spp.cores <- getOption("spp.cores")

                # if larger than 10 million do resampling
                if (ALL_TAGS > size_adj_thresh) {
                    # actually compute the mean over 100 random samplings
##			UNIQUE_TAGS_LibSizeadjusted<-round(mean(sapply(1:nsamp, FUN=function(x) {
##				return(length(unique(sample(chip.data, size=size_adj_thresh))))
##               	 	})))
			if(!is.null(spp.cores) && spp.cores>1 
				&& "package:multicore" %in% search() && nsamp > 1) {
				UNIQUE_TAGS_LibSizeadjusted <- round(mean(unlist(mclapply(1:nsamp, function(x) {
					sum(unlist(sapply(1:length(tags), function(xx) 
						length(unique(sample(tags[[xx]], chr_samp_tags[xx]))))))},
					mc.cores=spp.cores, mc.preschedule=F
				))))
			} else {
				UNIQUE_TAGS_LibSizeadjusted <- round(mean(sapply(1:nsamp, function(x) {
				sum(unlist(sapply(1:length(tags), function(xx) 
					length(unique(sample(tags[[xx]], chr_samp_tags[xx]))))))
				})))
			}
                } else {
                # if less than 10 million reads do resampling with replacement...
                ## (this is still under evaluation, it's not good) because the result is smaller than total NRF
                ## one possibility could be to take the best (higher) NRF among this one and the NRF compute on the entire taglist object
##                UNIQUE_TAGS_LibSizeadjusted<-round(mean(sapply(1:nsamp, FUN=function(x) {
##                    return(length(unique(sample(chip.data, size=size_adj_thresh, replace=TRUE))))
##                })))
			if(!is.null(spp.cores) && spp.cores>1
				&& "package:multicore" %in% search() && nsamp > 1) {
				UNIQUE_TAGS_LibSizeadjusted <- round(mean(unlist(mclapply(1:nsamp, function(x) {
					sum(unlist(sapply(1:length(tags), function(xx) 
						length(unique(sample(tags[[xx]], chr_samp_tags[xx], replace=TRUE))))))},
					mc.cores=spp.cores, mc.preschedule=F
				))))
				
			} else {
				UNIQUE_TAGS_LibSizeadjusted <- round(mean(sapply(1:nsamp, function(x) {
				sum(unlist(sapply(1:length(tags), function(xx) 
					length(unique(sample(tags[[xx]], chr_samp_tags[xx], replace=TRUE))))))
				})))
			}
                }

                NRF_LibSizeadjusted<-UNIQUE_TAGS_LibSizeadjusted/size_adj_thresh
	}
                # return a vector with NRF scores
                STATS_NRF<-list(ALL_TAGS=ALL_TAGS, UNIQUE_TAGS=UNIQUE_TAGS,
                	UNIQUE_TAGS_nostrand=UNIQUE_TAGS_nostrand, NRF=NRF,
                	NRF_nostrand=NRF_nostrand, NRF_LibSizeadjusted=NRF_LibSizeadjusted)
		
		##update NRF stats to qc
		if(is.null(qc))
			qc <<- list()
		qc$NRF <<- STATS_NRF
		cat("NRF=", round(STATS_NRF$NRF, 2), ", NRF (no strand)=", round(STATS_NRF$NRF_nostrand, 2), 
			", NRF (adjusted)=", round(STATS_NRF$NRF_LibSizeadjusted, 2), "\n", sep="")
                invisible(STATS_NRF)
        }
)

AlignedTags$methods(
	phantomPeak = function() {
		if(is.null(bd_chrtcs))
			stop("Please run 'get.cross.cor' first!")
		#bd_chrtcs <- AlignedTags$get.cross.cor
		
		# Phantom peak (shift = read_length) of cross correlation
		ph_peakidx <- which( ( bd_chrtcs$cross_cor$x >= ( read_length - round(2*bd_chrtcs$param$bin) ) ) & 
		( bd_chrtcs$cross_cor$x <= ( read_length + round(1.5*bd_chrtcs$param$bin) ) ) )
		ph_peakidx <- ph_peakidx[ which.max(bd_chrtcs$cross_cor$y[ph_peakidx]) ]
		phantom_cc <- bd_chrtcs$cross_cor[ph_peakidx,]
		
		# Minimum value of cross correlation in srange
		min_cc <- bd_chrtcs$cross_cor[ which.min(bd_chrtcs$cross_cor$y) , ]
		
		# Normalized Strand cross-correlation coefficient (NSC)
		NSC <- bd_chrtcs$peak$y / min_cc$y
		
		# Relative Strand Cross correlation Coefficient (RSC)
		RSC <- (bd_chrtcs$peak$y - min_cc$y) / (phantom_cc$y - min_cc$y)
		
		# Quality flag based on RSC value
		qflag <- NA
		if ( (RSC >= 0) & (RSC < 0.25) ) {
			qflag <- -2
		} else if ( (RSC >= 0.25) & (RSC < 0.5) ) {
			qflag <- -1
		} else if ( (RSC >= 0.5) & (RSC < 1) ) {
			qflag <- 0
		} else if ( (RSC >= 1) & (RSC < 1.5) ) {
			qflag <- 1
		} else if ( (RSC >= 1.5) ) {
			qflag <- 2
		}
		##update stats to qc
		if(is.null(qc))
			qc <<- list()
		phantom_peak <- list(phantom_cc=phantom_cc, NSC=NSC, RSC=RSC, quality_flag=qflag)
		qc$phantom_peak <<- phantom_peak
		cat("NSC=", round(NSC, 2), ", RSC=", round(RSC, 2), ", Quality flag: ", qflag, "\n", sep="")
		invisible(phantom_peak)
	}
)

