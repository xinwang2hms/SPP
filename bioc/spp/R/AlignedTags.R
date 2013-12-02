##-1. AlignedTags 
##-The class is designed for representing aligned ChIP-seq tags
##-note: this is a root class, which is inherited by specific classes
setClassUnion("smoothedTagDensity_Or_NULL", c("smoothedTagDensity", "NULL"))
AlignedTags = setRefClass(
	Class = "AlignedTags",
	fields = list(
		file = "char_Or_NULL", 				##read from file later	
		genome_build = "char_Or_NULL",			##e.g. mm9, hg19
##		read_length = "numeric_Or_NULL", 		##read length
		tags = "list_Or_NULL", 				
		quality = "list_Or_NULL",			##not in use actually 
		names = "list_Or_NULL", 			
		bd_chrtcs = "list_Or_NULL", 			##binding characteristics
		smoothed_density = "smoothedTagDensity_Or_NULL"
	)
)
setClassUnion("AlignedTags_Or_NULL", c("AlignedTags", "NULL"))
AlignedTags$methods(
	initialize = function(..., genome_build=NULL) {
		callSuper(...)
		if(is.null(.self$genome_build))
			genome_build <<- genome_build
		if(is.null(smoothed_density))
			smoothed_density <<- smoothedTagDensity(ChIP=.self)
		else {
			tmp <- .self$smoothed_density$.profile
			.self$smoothed_density <<- smoothedTagDensity(ChIP=.self, param=.self$smoothed_density$.param)
			.self$smoothed_density$.profile <<- tmp
		}
	}
)
##1. display a brief summary of the object
AlignedTags$methods(
	show = function(...) {
		nchrs <- length(tags)
		ntags <- sum(unlist(lapply(tags, length)))
		cat(as.character(class(.self)), "object:\n")
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
	compute.cross.cor = function(srange=c(50,500), bin=5, 
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
		abline(v=bd_chrtcs$peak$x, lty=2, col=2)				
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

















































