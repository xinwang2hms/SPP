
# given a tag vector (signed), identify and clean up (either remove or cap) singular positions that exceed local tag density
# vin - tag vector
# cap_fold - maximal fold over enrichment over local density allowed for a single tag position, at which the tag count is capped
# eliminate_fold - max fold enrichment that, when exceeded, results in exclusion of all the tags at that position (e.g. counted as anomaly)
# z_threshold - Z-score used to determine max allowed counts
filter.singular.positions.by.local.density <- function(tags, 
	window_size=200,cap_fold=4,eliminate_fold=10,z_threshold=3) {
  # tabulate tag positions
  if(length(tags)<2) 
	return(tags)
  
  tc <- table(tags)
  pos <- as.numeric(names(tc))
  storage.mode(pos) <- "double"
  tc <- as.integer(tc)
  storage.mode(tc) <- "integer"
  n <- length(pos)

  whs <- floor(window_size/2)
  
  storage.mode(n) <- storage.mode(whs) <- "integer"
  twc <- .Call("cwindow_n_tags_around",pos,tc,pos,whs)
  twc <- (twc-tc+1)/window_size # local density

  pv <- pnorm(z_threshold,lower.tail=F)
  # exclude
  max.counts <- qpois(pv,twc*eliminate_fold,lower.tail=F)
  tc[tc>max.counts] <- 0
  # cap
  max.counts <- qpois(pv,twc*cap_fold,lower.tail=F)
  ivi <- which(tc>max.counts)
  tc[ivi] <- max.counts[ivi]+1

  # reconstruct tag vector
  tv <- rep(pos,tc)
  to <- order(abs(tv))
  tv <- tv[to]
  return(tv)
}

# calculates tag strand cross-correlation for a range of shifts (on positive strand)
tag.scc <- function(tags,srange=c(50,250),bin=1,tt=NULL,llim=10) {
	if(is.null(tt)) {
		tt <- table(sign(tags)*as.integer(floor(abs(tags)/bin+0.5)))
	}
	if(!is.null(llim)) { 
		l <- mean(tt) 
		tt <- tt[tt<llim*l] 
	}
	tc <- as.integer(names(tt))
	tt <- as.numeric(tt)

	pv <- tt 
	pv[tc<0]<-0
	nv <- tt 
	nv[tc>0]<-0

	pti <- which(tc>0)
	nti <- which(tc<0)

	ptc <- tc[pti]
	ntc <- (-1)*tc[nti]

	ptv <- tt[pti]
	ntv <- tt[nti]

	trng <- range(c(range(ptc),range(ntc)))
	l <- diff(trng)+1
	rm(tc,tt)

	mp <- sum(ptv)*bin/l 
	mn <- sum(ntv)*bin/l
	ptv <- ptv-mp 
	ntv <- ntv-mn
	ss <- sqrt((sum(ptv*ptv)+(l-length(ptv))*mp^2) * 
		(sum(ntv*ntv)+(l-length(ntv))*mn^2))

	t.cor <- function(s) {
		smi <- match(ptc+s,ntc)
		return((sum(ptv[!is.na(smi)]*ntv[na.omit(smi)]) -
				mn*sum(ptv[is.na(smi)]) - mp*sum(ntv[-na.omit(smi)]) +
				mp*mn*(l-length(ptv)-length(ntv) + 
				length(which(!is.na(smi)))))/ss)
	}
	shifts <- floor(seq(srange[1],srange[2],by=bin)/bin+0.5)
	scc <- unlist(lapply(shifts,t.cor))
	names(scc) <- shifts*bin
	return(scc)
}

# plot chromosome-acerage cross-correlation 
t.plotavcc <- function(ci, main=paste(ci,"chromosome average"), 
	ccl=tl.cc, return.ac=F, ttl=tl, plot=T, ... ) {
	
	cc <- ccl[[ci]]
	if(length(cc)==1) 
		return(cc[[1]])
	if(length(cc)==0) 
		return(c()) 
	ac <- do.call(rbind,cc)
	# omit NA chromosomes
	ina <- apply(ac,1,function(d) any(is.na(d)))

	tags <- ttl[[ci]]
	avw <- unlist(lapply(tags,length))
	avw <- avw/sum(avw)
	ac <- ac[!ina,]
	avw <- avw[!ina]
	ac <- apply(ac,2,function(x) sum(x*avw))
	if(plot) {
		m <- t.plotcc(ac, main=main, ...)
		if(!return.ac)
		return(m)
	}
	if(return.ac) 
		return(ac)
}
# given tag position vectors, find contigs of significant enrichment of signal over background
# thr - z score threshold
# mcs - minimal cluster size
# bg.weight - fraction by which background counts should be multipled
# min.tag.count.z will impose a poisson constraint based on randomized signal in parallel of background constaint (0 - no constraint)
tag.enrichment.clusters <- function(signal, background, wsize=200, 
	thr=3, mcs=1, bg.weight=1, min.tag.count.z=0, tag.av.den=NULL, 
	min.tag.count.thr=0, min.tag.count.ratio=4, either=F, tag.shift=146/2) {
	
	if(is.null(tag.av.den)) {
		tag.av.den <- length(signal)/diff(range(abs(signal)))
	}
	if(min.tag.count.z>0) {
		min.tag.count.thr <- qpois(pnorm(min.tag.count.z, lower.tail=F), 
			min.tag.count.ratio*tag.av.den*wsize, lower.tail=F)
	} else {
		min.tag.count.thr <- 0
	}
  
	#if(bg.weight!=1) {
	#  background <- sample(background,length(background)*(bg.weight),replace=T)
	#}
	# make up combined position, flag vectors
	pv <- abs(c(signal, background)+tag.shift)
	fv <- c(rep(1, length(signal)), rep(0, length(background)))
	po <- order(pv)
	pv <- pv[po]
	fv <- fv[po]

	#thr <- pnorm(thr,lower.tail=F)

	storage.mode(wsize) <- storage.mode(mcs) <- storage.mode(fv) <- "integer"
	storage.mode(thr) <- storage.mode(pv) <- "double"
	storage.mode(bg.weight) <- "double"
	storage.mode(min.tag.count.thr) <- "double"
	either <- as.integer(either)
	storage.mode(either) <- "integer"

	z <- .Call("find_poisson_enrichment_clusters", pv, fv, wsize, thr, 
		mcs, bg.weight, min.tag.count.thr, either)
	return(z)
}

# determine d1/d2 dataset size ratio. If background.density.scaling=F, the ratio of tag counts is returned.
# if background.density.scaling=T, regions of significant tag enrichment are masked prior to ratio calculation.
# bug? lacking argument tag.shift
dataset.density.ratio <- function(d1, d2, min.tag.count.z=4.3, wsize=1e3, 
	mcs=0, background.density.scaling=T, return.proportion=F, tag.shift=146/2) {
	
	if(!background.density.scaling) {
		return(sum(unlist(lapply(d1,length)))/sum(unlist(lapply(d2,length))))
	}

	chrl <- intersect(names(d1),names(d2))
	ntc <- do.call(rbind,lapply(chrl,function(chr) {
		x1 <- tag.enrichment.clusters(abs(d1[[chr]]), c(), wsize=wsize, 
			bg.weight=0, min.tag.count.z=min.tag.count.z, mcs=mcs, 
			either=F, tag.shift=tag.shift)
		x2 <- tag.enrichment.clusters(abs(d2[[chr]]), c(), wsize=wsize, 
			bg.weight=0, min.tag.count.z=min.tag.count.z, mcs=mcs, 
			either=F, tag.shift=tag.shift)
		return(c(length(which(points.within(abs(d1[[chr]]), 
			c(x1$s,x2$s)-wsize/2, c(x1$e,x2$e)+wsize/2)==-1)), 
			length(which(points.within(abs(d2[[chr]]), 
			c(x1$s,x2$s)-wsize/2, c(x1$e,x2$e)+wsize/2)==-1))))
	}))
	ntcs <- apply(ntc,2,sum)
	#print(ntcs/c(sum(unlist(lapply(d1,length))),sum(unlist(lapply(d2,length)))))
	return(ntcs[1]/ntcs[2])
}
# determine membership of points in fragments
points.within <- function(x, fs, fe, return.list=F, return.unique=F, 
	sorted=F, return.point.counts=F) {
	
	if(is.null(x) | length(x) < 1) 
		return(c()) 
	if(!sorted) {
		ox <- rank(x,ties="first")
		x <- sort(x)
	}

	se <- c(fs,fe)
	fi <- seq(1:length(fs))
	fi <- c(fi,-1*fi)

	fi <- fi[order(se)]
	se <- sort(se)

	storage.mode(x) <- storage.mode(fi) <- storage.mode(se) <- "integer"
	if(return.unique) 
		iu <- 1
	else 
		iu <- 0
	if(return.list) 
		il <- 1
	else 
		il <- 0
	if(return.point.counts) 
		rpc <- 1 
	else 
		rpc <- 0
	storage.mode(iu) <- storage.mode(il) <- storage.mode(rpc) <- "integer"
	result <- .Call("points_within", x, se, fi, il, iu, rpc)
	if(!sorted & !return.point.counts) {
		result <- result[ox]
	}
	return(result)
}
# returns effective size of the dataset based on the same logic as dataset.density.ratio
# bug? lacking argument tag.shift
dataset.density.size <- function(d1, min.tag.count.z=4.3, wsize=1e3, 
	mcs=0, background.density.scaling=T, tag.shift=146/2) {
	
	if(!background.density.scaling) {
		return(sum(unlist(lapply(d1,length))))
	}

	chrl <- names(d1)
	ntc <- lapply(chrl,function(chr) {
		x1 <- tag.enrichment.clusters(abs(d1[[chr]]), c(), wsize=wsize, 
			bg.weight=0, min.tag.count.z=min.tag.count.z, mcs=mcs, 
			either=F, tag.shift=tag.shift)
		return(length(which(points.within(abs(d1[[chr]]), x1$s-wsize/2, 
			x1$e+wsize/2)==-1)))
	})
	return(sum(unlist(ntc)))
}
# calculate cumulative density based on sum of scaled gaussian curves
# (by Michael Tolstorukov)
#
# vin - input vector; bw -- standard deviation, dw-gaussina cutoff in stdev; dout - output "density")
# output - if return.x=F vector of cumulative density values corresponding to integer positions described by range(vin)
# output - if return.x=T a data structure with $x and $y corresponding to the cumulative density
# optional match.wt.f is a function that will return weights for a tag vector
densum <- function(vin, bw=5, dw=3, match.wt.f=NULL, return.x=T, 
	from=min(vin), to=max(vin), step=1)    {
	# construct vector of unique tags and their counts
	tc <- table(vin[vin>=from & vin<=to])
	pos <- as.numeric(names(tc))
	storage.mode(pos) <- "double"
	tc <- as.numeric(tc)
	storage.mode(tc) <- "double"
	n <- length(pos)
	# weight counts
	if(!is.null(match.wt.f)) {
		tc <- tc*match.wt.f(pos)
	}

	rng <- c(from,to)
	if(rng[1]<0) 
		stop("range extends into negative values") 
	if(range(pos)[1]<0) 
		stop("position vector contains negative values") 

	storage.mode(n) <- storage.mode(rng) <- storage.mode(bw) <- 
		storage.mode(dw) <- storage.mode(step) <- "integer"

	spos <- rng[1]
	storage.mode(spos) <- "double"

	dlength <- floor((rng[2] - rng[1])/step) + 1 # length of output array
	if(dlength<1) { stop("zero data range") }
	dout <- numeric(dlength); storage.mode(dout) <- "double"
	storage.mode(dlength) <- "integer"
	.C("cdensum",n,pos,tc,spos,bw,dw,dlength,step,dout,DUP=F)

	if(return.x) {
		return(list(x=c(rng[1],rng[1]+step*(dlength-1)), y=dout, 
			step=step))
	} else {
		return(dout)
	}
}
# calculate smoothed tag density, optionally subtracting the background
get.smoothed.tag.density <- function(signal.tags,control.tags=NULL,
	bandwidth=150,bg.weight=NULL,tag.shift=146/2,step=round(bandwidth/3),
	background.density.scaling=T,rngl=NULL,chrl=NULL,scale.by.dataset.size=F) {
 
    if(is.null(chrl)) {	
		chrl <- names(signal.tags)
		names(chrl) <- chrl
    }
    if(!is.null(control.tags)) {
		if(is.null(bg.weight))
			bg.weight <- dataset.density.ratio(signal.tags, control.tags, 
				background.density.scaling=background.density.scaling)
    }
    
    if(scale.by.dataset.size) {
        den.scaling <- 1/(dataset.density.size(signal.tags, 
			background.density.scaling=background.density.scaling)/1e6)
    } else {
        den.scaling <- 1
    }
##
	.internal.get.tag.density <- function(chr) {
        ad <- abs(signal.tags[[chr]]+tag.shift)
        rng <- NULL
        if(!is.null(rngl)) {
            rng <- rngl[[chr]]
        }
        if(is.null(rng)) {
            rng <- range(ad)
        }      
        ds <- densum(ad, bw=bandwidth, from=rng[1], to=rng[2], 
			return.x=T, step=step)
        if(!is.null(control.tags)) {
            if(!is.null(control.tags[[chr]])) {
                bsd <- densum(abs(control.tags[[chr]]+tag.shift), 
					bw=bandwidth, from=rng[1], to=rng[2], return.x=F, 
					step=step)
                ds$y <- ds$y-bsd*bg.weight
            }
        }
        return(data.frame(x=seq(ds$x[1], ds$x[2], by=step), 
			y=den.scaling*ds$y))
	}		
	#check if parallel
	spp.cores <- getOption("spp.cores")	
	if(!is.null(spp.cores) && spp.cores>1 
		&& "package:multicore" %in% search() && length(chrl) > 1) {
		ldens <- mclapply(chrl, .internal.get.tag.density, 
			mc.cores=spp.cores, mc.preschedule=F)
	} else 
		ldens <- lapply(chrl, .internal.get.tag.density)
    return(ldens)
}
# count tags within sliding window of a specified size
# vin - tag vector (postive values, pre-shifted)
# window.size/window.step - window characteristics
# tv - optional, pre-sorted, pre-trimmed tag vector
window.tag.count <- function(vin,window.size,window.step=1,return.x=T,from=min(vin)+floor(window.size/2),to=max(vin)-floor(window.size/2),tv=NULL) {
  whs <- floor(window.size/2);
  # select tags with margins
  if(is.null(tv)) {
    tv <- sort(vin[vin>=from-whs-1 & vin<=to+whs+1])
  }
  storage.mode(tv) <- "double";
  n <- length(tv)
  nsteps <- ceiling((to-from)/window.step);
  
  storage.mode(n) <- storage.mode(nsteps) <- storage.mode(window.size) <- storage.mode(window.step) <- "integer";
  
  spos <- from; storage.mode(spos) <- "double";

  if(nsteps<1) { stop("zero data range") }
  #dout <- integer(nsteps); storage.mode(dout) <- "integer";
  #.C("window_n_tags",n,tv,spos,window.size,window.step,nsteps,dout,DUP=F);
  dout <- .Call("cwindow_n_tags",tv,spos,window.size,window.step,nsteps);
  
  if(return.x) {
    return(list(x=c(from,from+(nsteps-1)*window.step),y=dout,step=window.step))
  } else {
    return(dout)
  }
}
# count tags in windows around specified positions (pos)
window.tag.count.around <- function(vin,window.size,pos,return.x=T,tc=NULL,sorted=F) {
  if(is.null(tc)) {
    tc <- table(vin);
  }
  if(!sorted) {
    op <- rank(pos);
    pos <- sort(pos);
  }
  storage.mode(pos) <- "double";
  tpos <- as.integer(names(tc)); storage.mode(tpos) <- "double";
  tc <- as.integer(tc); storage.mode(tc) <- "integer";
  
  whs <- floor(window.size/2);
  
  storage.mode(whs) <- "integer";
  twc <- .Call("cwindow_n_tags_around",tpos,tc,pos,whs);
  if(return.x) {
    if(sorted) {
      return(data.frame(x=pos,y=twc));
    } else {
      return(data.frame(x=pos[op],y=twc[op]));
    }
  } else {
    if(sorted) {
      return(twc);
    } else {
      return(twc[op]);
    }
  }
}
# calculates enrichment bounds using multiple background scales
# ft - foreground tags (pre-shifted, positive)
# bt - background tags
# fws - foreground window size
# bwsl - background window size list
# step - window step
# rng - from/to coordinates (to will be adjusted according to step)
#
# returns: a list with $x ($s $e $step), $lb vector and $mle vector ($ub if calculate.upper.bound=T)
mbs.enrichment.bounds <- function(ft,bt,fws,bwsl,step=1,rng=NULL,alpha=0.05,calculate.upper.bound=F,bg.weight=length(ft)/length(bt),use.most.informative.scale=F,quick.calculation=F,pos=NULL) {
  # determine range
  if(is.null(rng)) {
    rng <- range(range(ft));
  }
  # foreground counts
  if(is.null(pos)) {
    fwc <- window.tag.count(ft,fws,window.step=step,from=rng[1],to=rng[2],return.x=T);
  } else {
    fwc <- window.tag.count.around(ft,fws,pos,return.x=T)
  }
  fwc$y <- fwc$y+0.5;

  zal <- qnorm(alpha/2,lower.tail=F);

  # background counts
  bt <- sort(bt);
  if(!is.null(pos)) {
    tc <- table(bt);
  }
  bgcm <- lapply(bwsl,function(bgws) {
    if(is.null(pos)) {
      window.tag.count(bt,bgws,window.step=step,from=rng[1],to=rng[2],return.x=F,tv=bt)+0.5;
    } else {
      window.tag.count.around(bt,bgws,pos,return.x=F,tc=tc)+0.5
    }
  })
  if(!is.null(pos)) {
    rm(tc);
  }

  # pick most informative scale
  if(use.most.informative.scale) {
    bgcm <- t(do.call(cbind,bgcm))
    isi <- max.col(t((bgcm)/(bwsl/fws))) # add pseudo-counts to select lowest scale in case of a tie

    bgc <- c(bgcm)[isi+dim(bgcm)[1]*(c(1:length(isi))-1)]

    if(quick.calculation) {
      rte <- fwc$y+bgc-0.25*zal*zal; rte[rte<0] <- 0;
      dn <- bgc - 0.25*zal*zal;
      lbm=(sqrt(fwc$y*bgc) - 0.5*zal*sqrt(rte))/dn;
      ivi <- which(lbm<0);
      lbm <- lbm*lbm*bwsl[isi]/fws/bg.weight;
      lbm[rte<=0] <- 1;
      lbm[dn<=0] <- 1;
      lbm[ivi] <- 1;
    } else {
      lbm <- (fwc$y/bgc)*qf(1-alpha/2,2*fwc$y,2*bgc,lower.tail=F)*bwsl[isi]/fws/bg.weight;
    }
    
    mle <- fwc$y/bgc*bwsl[isi]/fws/bg.weight; mle[is.nan(mle)] <- Inf; mle[is.na(mle)] <- Inf;
    
    rl <- list(x=list(s=fwc$x[1],e=fwc$x[2],step=fwc$step),lb=lbm,mle=mle);
    
    if(calculate.upper.bound) {
      isi <- max.col(t((-bgcm)/(bwsl/fws))) # add pseudo-counts to select highest scale in case of a tie
      bgc <- c(bgcm)[isi+dim(bgcm)[1]*(c(1:length(isi))-1)]

      if(quick.calculation) {
        ubm=(sqrt(fwc$y*bgc) + 0.5*zal*sqrt(rte))/dn;
        ivi <- which(ubm<0);
        ubm <- ubm*ubm*bwsl[isi]/fws/bg.weight;
        ubm[rte<=0] <- 1;
        ubm[ivi] <- 1;
        lbm[dn<=0] <- 1;
      } else {
        ubm <- (fwc$y/bgc)*qf(alpha/2,2*fwc$y,2*bgc,lower.tail=F)*bwsl[isi]/fws/bg.weight;
      }
      rl <- c(rl,list(ub=ubm));
    }
    return(rl);
    
  } else {
    # determine lower bounds
    lbm <- lapply(c(1:length(bgcm)),function(i) {
      nbg <- bgcm[[i]];
      if(quick.calculation) {
        rte <- fwc$y+nbg-0.25*zal*zal; rte[rte<0] <- 0;
        dn <- (nbg - 0.25*zal*zal);
        lbm=(sqrt(fwc$y*nbg) - 0.5*zal*sqrt(rte))/dn;
        ivi <- which(lbm<0);  
        lbm <- lbm*lbm*bwsl[i]/fws/bg.weight;
        lbm[rte<=0] <- 1;
        lbm[dn<=0] <- 1;
        lbm[ivi] <- 1;
        return(lbm);
      } else {
        return((fwc$y/nbg)*qf(1-alpha/2,2*fwc$y,2*nbg,lower.tail=F)*bwsl[i]/fws/bg.weight);
      }
    })
    lbm <- do.call(pmin,lbm);

    # calculate mle
    #mle <- do.call(pmin,lapply(bgcm,function(bgc) fwc/bgc))
    mle <- do.call(pmin,lapply(c(1:length(bgcm)),function(i) {
      bgc <- bgcm[[i]];
      x <- fwc$y/bgc*bwsl[i]/fws/bg.weight; x[is.nan(x)] <- Inf; x[is.na(x)] <- Inf; return(x);
    }))

    rl <- list(x=list(s=fwc$x[1],e=fwc$x[2],step=fwc$step),lb=lbm,mle=mle);
    
    if(calculate.upper.bound) {
      # determine upper bound
      ubm <- lapply(c(1:length(bgcm)),function(i) {
        nbg <- bgcm[[i]];
        if(quick.calculation) {
          rte <- fwc$y+nbg-0.25*zal*zal; rte[rte<0] <- 0;
          dn <- (nbg - 0.25*zal*zal);
          ubm=(sqrt(fwc$y*nbg) + 0.5*zal*sqrt(rte))/dn;
          ivi <- which(ubm<0);  
          ubm <- ubm*ubm*bwsl[i]/fws/bg.weight;
          ubm[rte<=0] <- 1;
          ubm[dn<=0] <- 1;
          ubm[ivi] <- 1;
          return(ubm);
        } else {
          return((fwc$y/nbg)*qf(alpha/2,2*fwc$y,2*nbg,lower.tail=F)*bwsl[i]/fws/bg.weight);
        }
      })
      ubm <- do.call(pmax,ubm);
      rl <- c(rl,list(ub=ubm));
    }

    return(rl);
  }
}
# find regions of significant tag enrichment
find.significantly.enriched.regions <- function(signal.data,control.data, 
	window.size=500,multiplier=1,z.thr=3,mcs=0,debug=F,
	background.density.scaling=T,masking.window.size=window.size,
	poisson.z=0,poisson.ratio=4,either=F,tag.shift=146/2,bg.weight=NULL) {
	
	if(is.null(bg.weight)) {
		bg.weight <- dataset.density.ratio(signal.data, control.data, 
			background.density.scaling=background.density.scaling)
	}

	if(debug) {
		cat("bg.weight=",bg.weight,"\n")
	}
	chrl <- names(signal.data)
	names(chrl) <- chrl
	.internal.tag.enrichment.clusters <- function(chr) {
		d <- tag.enrichment.clusters(signal.data[[chr]], 
			control.data[[chr]], bg.weight=bg.weight*multiplier, 
			thr=z.thr, wsize=window.size, mcs=mcs, 
			min.tag.count.z=poisson.z, min.tag.count.ratio=poisson.ratio,
			either=either,tag.shift=tag.shift)
		d$s <- d$s-masking.window.size/2
		d$e <- d$e+masking.window.size/2
		return(d)
	}
	#check if parallel
	spp.cores <- getOption("spp.cores")
	if(!is.null(spp.cores) && spp.cores>1 
		&& "package:multicore" %in% search() && length(chrl) > 1) {
		tec <- mclapply(chrl, .internal.tag.enrichment.clusters, 
			mc.cores=spp.cores, mc.preschedule=FALSE)
	} else {		
		tec <- lapply(chrl, .internal.tag.enrichment.clusters)
	}
	return(tec)
}
# returns intersection of multiple region sets
# each regionset needs to contain $s, $e and optional $v column
regionset.intersection.c <- function(rsl,max.val=-1,do.union=F) {
  # translate into position/flag form
  rfl <- lapply(rsl,function(rs) {
    rp <- c(rs$s,rs$e); rf <- c(rep(c(1,-1),each=length(rs$s)));
    
    ro <- order(rp);
    rp <- rp[ro]; rf <- rf[ro];
    if(!is.null(rs$v)) {
      rv <- c(rs$v,rs$v)[ro];
      return(data.frame(p=as.numeric(rp),f=as.integer(rf),v=as.numeric(rv)));
    } else {
      return(data.frame(p=as.numeric(rp),f=as.integer(rf)));
    }
  })
  rfd <- data.frame(do.call(rbind,lapply(1:length(rfl),function(i) {
    d <- rfl[[i]]; d$f <- d$f*i; return(d);
  })))
  rfd <- rfd[order(rfd$p),];
  if(is.null(rfd$v)) { max.val <- 0; }
  if(do.union) { ur <- 1; } else { ur <- 0; }; 
  rl <- .Call("region_intersection",as.integer(length(rfl)),as.numeric(rfd$p),as.integer(rfd$f),as.numeric(rfd$v),as.integer(max.val),as.integer(ur));
  return(data.frame(do.call(cbind,rl)));
}

# estimates threshold, calculates predictions on complete data and randomized data
# input: tvl
# control - a list of control tag datasets
# no randomization is done if control is supplied
# return.rtp - return randomized tag peaks - do not fit thresholds or do actual predictions
# topN - use min threshold to do a run, return topN peaks from entire genome
# threshold - specify a user-defined threshold
lwcc.prediction <- function(tvl,e.value=NULL, fdr=0.01, chrl=names(tvl), 
	min.thr=0, n.randomizations=1, shuffle.window=1, debug=T, 
	predict.on.random=F, shuffle.both.strands=T,strand.shuffle.only=F, 
	return.rtp=F, control=NULL, print.level=0, threshold=NULL, topN=NULL, 
	bg.tl=NULL, tec.filter=T, tec.window.size=1e3,tec.z=3, 
	tec.masking.window.size=tec.window.size, tec.poisson.z=3,
	tec.poisson.ratio=4, bg.reverse=T, return.control.predictions=F, 
	return.core.data=F, background.density.scaling=T, bg.weight=NULL, 
	tag.shift=146/2, ... ) {

	control.predictions <- NULL
	core.data <- list()

	if(!is.null(bg.tl) && tec.filter) {
		if(debug) 
			cat("finding background exclusion regions ... ")
		tec <- find.significantly.enriched.regions(bg.tl, tvl, 
			window.size=tec.window.size, z.thr=tec.z, 
			masking.window.size=tec.masking.window.size, 
			poisson.z=tec.poisson.z, poisson.ratio=tec.poisson.ratio, 
			background.density.scaling=background.density.scaling, tag.shift=tag.shift, 
			either=T)
		if(return.core.data) 
			core.data <- c(core.data,list(tec=tec))
		if(debug) 
			cat("[done]\n")
	}

  
	if(is.null(threshold) && is.null(topN)) { # threshold determination is needed
		# generate control predictions
		if(!is.null(control)) {
			if(debug) 
				cat("determining peaks on provided", length(control), "control datasets:\n")
			##check what to use as background tags
			if(!is.null(bg.tl)) {
				if(bg.reverse) {
					if(debug) 
					cat("using reversed signal for FDR calculations\n")
					rbg.tl <- tvl
				} else {
					if(debug) 
					cat("generating randomized (within chromosome) background ... ")
					rbg.tl <- lapply(bg.tl, function(d) {
						if(length(d)<2) 
							return(d)
						rng <- range(abs(d))
						rd <- round(runif(length(d),rng[1],rng[2]))
						nrd <- sample(1:length(rd),length(which(d<0)))
						rd[nrd] <- rd[nrd]*(-1)
						return(rd)
					})
					if(debug) 
						cat("[done]\n")
				}
			} else {
				rbg.tl <- NULL
			}
			n.randomizations <- length(control)
			#signal.size <- sum(unlist(lapply(tvl,length)));
			rtp <- lapply(control, function(d) {
				# calculate tag.weight
				#tag.weight <- sum(unlist(lapply(tvl,length)))/sum(unlist(lapply(d,length)));
				tag.weight <- dataset.density.ratio(tvl, d, 
					background.density.scaling=background.density.scaling)
				#cat("tag.weight=",tag.weight," ");
				return(window.call.mirror.binding(d, min.thr=min.thr, 
					tag.weight=tag.weight, bg.tl=rbg.tl, debug=debug, 
					round.up=T,background.density.scaling=background.density.scaling, ...))
				#return(window.call.mirror.binding(d,min.thr=min.thr, method=tag.wtd,wsize=200,bg.tl=control.data,window.size=window.size,debug=T,min.dist=min.dist,cluster=cluster))
			})
			if(return.core.data) {
				core.data <- c(core.data, list(rtp.unfiltered=rtp))
			}
			if(tec.filter) {
				if(debug) 
					cat("excluding systematic background anomalies ... ")
				rtp <- lapply(rtp, filter.binding.sites, tec, exclude=T)
				if(debug) 
					cat("[done]\n")
			}
		} else {
			if(debug) 
				cat("determining peaks on ", n.randomizations, "randomized datasets:\n")
			rtp <- lapply(1:n.randomizations, function(i) {
				rd <- generate.randomized.data(tvl, shuffle.window=shuffle.window, 
					shuffle.both.strands=shuffle.both.strands, 
					strand.shuffle.only=strand.shuffle.only)
				return(window.call.mirror.binding(rd, min.thr=min.thr, 
					bg.tl=bg.tl, debug=debug, ...))
				#return(window.call.mirror.binding(rd,min.thr=min.thr, method=tag.wtd,wsize=200,bg.tl=control.data,window.size=window.size,debug=T,min.dist=min.dist))
			})
		}
		if(return.control.predictions) {
			control.predictions <- rtp
		} 
		rtp <- do.call(rbind, lapply(rtp, function(d) do.call(rbind,d)))	# merge tables
		
		# generate real data predictions
		if(debug) 
			cat("determining peaks on real data:\n")
		npl <- window.call.mirror.binding(tvl, min.thr=min.thr, bg.tl=bg.tl, 
			debug=debug, background.density.scaling=background.density.scaling, ...)
		#npl <- window.call.mirror.binding(tvl,min.thr=min.thr, method=tag.wtd,wsize=200,bg.tl=control.data,window.size=window.size,debug=T,min.dist=min.dist,cluster=cluster);
		if(return.core.data) {
			core.data <- c(core.data, list(npl.unfiltered=npl))
		}

		if(!is.null(bg.tl) && tec.filter) {
			if(debug) 
				cat("excluding systematic background anomalies ... ")
			npl <- filter.binding.sites(npl, tec, exclude=T)
			if(debug) 
				cat("[done]\n")
		}

		# calculate E-value and FDRs for all of the peaks
		if(debug) 
			cat("calculating statistical thresholds\n")
		chrl <- names(npl)
		names(chrl) <- chrl
		npld <- do.call(rbind, lapply(names(npl),function(chr) {
			k <- npl[[chr]]
			if(!is.null(k) && dim(k)[1]>0) {
				k$chr <- rep(chr,dim(k)[1])
			}
			return(k)
		}))
		npld <- cbind(npld, get.eval.fdr.vectors(npld$y, rtp$y))
		# correct for n.randomizations
		npld$fdr <- npld$fdr/n.randomizations
		npld$evalue <- npld$evalue/n.randomizations

		if(return.core.data) {
			core.data <- c(core.data,list(npld=npld))
		}

		# determine actual thresholds
		if(is.null(e.value)) {
			if(is.null(fdr)) 
				fdr <- 0.01
			thr <- list(root=min(npld$y[npld$fdr<=fdr]), type="FDR", fdr=fdr)
			if(debug) 
				cat("FDR",fdr,"threshold=",thr$root,"\n")
		} else {
			# determine threshold based on e-value
			thr <- list(root=min(npld$y[npld$evalue<=e.value]), type="Evalue", 
				e.value=e.value)
			if(debug) 
				cat("E-value", e.value, "threshold=", thr$root, "\n")
		}
		npld <- npld[npld$y>=thr$root,]
		if(dim(npld)[1]>0) {
			npl <- tapply(c(1:dim(npld)[1]), as.factor(npld$chr), 
				function(ii) {
					df <- npld[ii,]
					df$chr <- NULL
					return(df) 
				})
		} else {
			npl <- list()
		}
	} else {
		if(is.null(threshold)) {
			thr <- list(root=min.thr, type="minimal")
		} else {
			thr <- list(root=threshold, type="user specified")
		}

		cat("calling binding positions using", thr$type, 
			"threshold (", thr$root, ") :\n")
		npl <- window.call.mirror.binding(tvl=tvl, min.thr=thr$root, 
			bg.tl=bg.tl, debug=debug, ...)
		if(!is.null(bg.tl) && tec.filter) {
			if(debug) 
				cat("excluding systematic background anomalies ... ")
			npl <- filter.binding.sites(npl,tec,exclude=T)
			if(debug) 
				cat("[done]\n")
		}

		if(!is.null(topN)) {
			# determine threshold based on topN peaks
			ay <- unlist(lapply(npl,function(d) d$y))
			if(length(ay)>topN) {
				thr <- list(root=sort(ay,decreasing=T)[topN],type="topN",topN=topN)
				cat(paste("determined topN threshold :",thr$root,"\n"))    
				npl <- lapply(npl,function(d) d[d$y>thr$root,])
			}
		}
	}

	if(return.core.data) {
		return(c(list(npl=npl,thr=thr), core.data))
	}
	if(return.control.predictions && !is.null(control.predictions)) {
		return(list(npl=npl,thr=thr,control.predictions=control.predictions))
	}
	return(list(npl=npl,thr=thr))
}
# determine mirror-based binding positions using sliding window along each chromosome
# extra parameters are passed on to call.nucleosomes()
window.call.mirror.binding <- function(tvl,window.size=4e7, debug=T, 
	bg.tl=NULL, mask.tl=NULL, background.density.scaling=T, ...) {
	
	chrl <- names(tvl)
	# determine bg.weight
	if(!is.null(bg.tl)) {
		bg.weight <- dataset.density.ratio(tvl, bg.tl, 
			background.density.scaling=background.density.scaling)
	} else {
		bg.weight <- NULL
	}
	if(debug) {
		cat("bg.weight=", bg.weight, " ")
	}
  
	names(chrl) <- chrl
	#check if parallel
	spp.cores <- getOption("spp.cores")
	if(!is.null(spp.cores) && spp.cores>1 
			&& "package:multicore" %in% search() && length(chrl) > 1) {
		# add bg.ctv and mask.ctv to parallel call
		tvll <- lapply(chrl,function(chr) {
			bg.ctv <- NULL
			if(!is.null(bg.tl)) {
				bg.ctv <- bg.tl[[chr]]
			}
			mask.ctv <- NULL
			if(!is.null(mask.tl)) {
				mask.ctv <- mask.tl[[chr]]
			}
			return(list(ctv=tvl[[chr]], bg.ctv=bg.ctv, mask.ctv=mask.ctv, chr=chr))
		})
##		bl <- clusterApplyLB(cluster, tvll, window.chr.call.mirror.binding, 
##			window.size=window.size, debug=debug, bg.weight=bg.weight, ...)
		bl <- mclapply(tvll, window.chr.call.mirror.binding, 
			window.size=window.size, debug=debug, bg.weight=bg.weight, ..., 
			mc.preschedule = FALSE, mc.cores=spp.cores)
		names(bl) <- chrl
		return(bl)	
	} else {
		return(lapply(chrl,function(chr) {
			bg.ctv <- NULL
			if(!is.null(bg.tl)) { 
				bg.ctv <- bg.tl[[chr]]
			}
			mask.ctv <- NULL
			if(!is.null(mask.tl)) { 
				mask.ctv <- mask.tl[[chr]]
			}
	  
			window.chr.call.mirror.binding(list(ctv=tvl[[chr]], 
				bg.ctv=bg.ctv, mask.ctv=mask.ctv, chr=chr), 
				window.size=window.size, chr=chr, debug=debug, 
				bg.weight=bg.weight, bg.ctv=bg.ctv, mask.ctv=mask.ctv, ...)
		}))
	}
}

window.chr.call.mirror.binding <- function(ctvl, window.size, debug=T, 
	chr="NA", method=tag.wtd, bg.ctv=NULL, mask.ctv=NULL, ...) {
	chr <- ctvl$chr
	ctv <- ctvl$ctv
	bg.ctv <- ctvl$bg.ctv
	mask.ctv <- ctvl$mask.ctv
	if(is.null(ctv)) 
		return(data.frame(x=c(),y=c())) 
	if(length(ctv)<2) 
		return(data.frame(x=c(),y=c())) 

	dr <- range(unlist(lapply(ctv,function(x) range(abs(x)))))
	n.windows <- ceiling(diff(dr)/window.size)


	pinfo <- c()
	if(debug) {
		cat(paste("processing ",chr," in ",n.windows," steps [",sep=""))
	}
	for(i in 1:n.windows) {
		s <- dr[1]+(i-1)*window.size
		npn <- method(s=s, e=s+window.size,ctv=ctv, return.peaks=T, 
			bg.ctv=bg.ctv, mask.ctv=mask.ctv, ... )
		if(length(npn) > 0) { pinfo <- rbind(pinfo,npn)  }
		if(debug) {
		  cat(".")
		}
	}
	if(debug) {
		cat(paste("] done (",dim(pinfo)[1],"positions)\n"))
	} 
##	else {
##		cat(".")
##	}
	return(data.frame(x=pinfo[,1],y=pinfo[,2]))
}

tag.wtd <- function(ctv,s,e,return.peaks=T, bg.ctv=NULL,  mask.ctv=NULL, ...) {
  x <- ctv[ctv>=s & ctv<=e];
  y <- (-1)*ctv[ctv<=-s & ctv>=-e];

  if(!is.null(bg.ctv)) {
    bg.x <- bg.ctv[bg.ctv>=s & bg.ctv<=e];
    bg.y <- (-1)*bg.ctv[bg.ctv<=-s & bg.ctv>=-e];
  } else {
    bg.x <- bg.y <- NULL;
  }

  if(!is.null(mask.ctv)) {
    mask.x <- mask.ctv[mask.ctv>=s & mask.ctv<=e];
    mask.y <- (-1)*mask.ctv[mask.ctv<=-s & mask.ctv>=-e];
  } else {
    mask.x <- mask.y <- NULL;
  }

  if(length(x)==0 | length(y) ==0) {
    if(return.peaks) {
      return(data.frame(x=c(),y=c()));
    } else {
      rx <- range(c(x,y));
      return(list(x=rx,y=numeric(diff(rx)+1)));
    }
  } else {
    return(wtd(x,y,s,e,return.peaks=return.peaks,  bg.x=bg.x,bg.y=bg.y, mask.x=mask.x,mask.y=mask.y, ...))
  }
}

# window tag difference method
wtd <- function(x,y,s,e,whs=200,return.peaks=T,min.thr=5,min.dist=200,step=1,direct.count=F,tag.weight=1,bg.x=NULL,bg.y=NULL,bg.weight=1,mask.x=NULL,mask.y=NULL,ignore.masking=F, bg.whs=whs, round.up=F, ...) {
  ignore.masking <- ignore.masking | (is.null(mask.x) & is.null(mask.y));
  if(step>1) {
    x <- floor(x/step+0.5); y <- floor(y/step+0.5)
    
    if(!is.null(bg.x)) {
      bg.x <- floor(bg.x/step+0.5); bg.y <- floor(bg.y/step+0.5)  
    }
    
    if(!is.null(mask.x)) {
      mask.x <- floor(mask.x/step+0.5); mask.y <- floor(mask.y/step+0.5)  
    }

    
    whs <- floor(whs/step+0.5);
    bg.whs <- floor(bg.whs/step+0.5);
    min.dist <- floor(min.dist/step +0.5);
    s <- floor(s/step+0.5)
    e <- floor(e/step+0.5)
  }

  # scale bg.weight, since within calculation they are considered independent
  bg.weight <- bg.weight*tag.weight;

  rx <- c(s-whs,e+whs);

  # compile tag vectors
  xt <- table(x);
  xh <- integer(diff(rx)+1);
  xh[as.integer(names(xt))-rx[1]+1] <- as.integer(xt);

  yt <- table(y);
  yh <- integer(diff(rx)+1);
  yh[as.integer(names(yt))-rx[1]+1] <- as.integer(yt);

  # compile background vectors
  if(!is.null(bg.x) && length(bg.x)>0) {
    bg.subtract <- 1;

    bg.xt <- table(bg.x);
    bg.xh <- integer(diff(rx)+1);
    bg.xh[as.integer(names(bg.xt))-rx[1]+1] <- as.integer(bg.xt);
    rm(bg.xt);

    bg.yt <- table(bg.y);
    bg.yh <- integer(diff(rx)+1);
    bg.yh[as.integer(names(bg.yt))-rx[1]+1] <- as.integer(bg.yt);
    rm(bg.yt);

    # adjust bg.weight according to bg.whs
    if(bg.whs!=whs) {
      bg.weight <- bg.weight*whs/bg.whs;
    }
  } else {
    bg.subtract <- 0;
    bg.xh <- bg.yh <- c();
  }

  # record masked positions
  if(!ignore.masking) {
    if(!is.null(mask.x) && length(mask.x)>0) {
      mvx <- unique(mask.x); mvx <- setdiff(mvx,as.numeric(names(xt)));
      mvx <- mvx[mvx>=rx[1] & mvx<=rx[2]];
      xh[mvx-rx[1]+1] <- -1;
    }

    if(!is.null(mask.y) && length(mask.y)>0) {
      mvy <- unique(mask.y); mvy <- setdiff(mvy,as.numeric(names(yt)));
      mvy <- mvy[mvy>=rx[1] & mvy<=rx[2]];
      yh[mvy-rx[1]+1] <- -1;
    }
  }

  rm(xt,yt);

  if(round.up) { round.up <- 1; } else { round.up <- 0; }
  
  storage.mode(xh) <- storage.mode(yh) <- "integer";
  storage.mode(bg.xh) <- storage.mode(bg.yh) <- "integer";
  nx <- length(xh);   storage.mode(nx) <- storage.mode(whs) <- storage.mode(bg.whs) <- "integer";
  rp <- as.integer(return.peaks);
  dcon <- as.integer(direct.count);
  storage.mode(rp) <- storage.mode(min.dist) <- "integer";
  storage.mode(min.thr) <- "double";
  storage.mode(dcon) <- "integer";
  storage.mode(tag.weight) <- "double";
  storage.mode(bg.weight) <- "double";
  storage.mode(bg.subtract) <- "integer";
  storage.mode(round.up) <- "integer";
  im <- as.integer(ignore.masking);
  storage.mode(im) <- "integer";
  z <- .Call("wtd",xh,yh,whs,rp,min.dist,min.thr,dcon,tag.weight,im,bg.subtract,bg.xh,bg.yh,bg.whs,bg.weight,round.up);
  if(return.peaks) {
    return(data.frame(x=(z$x+rx[1])*step,y=z$v));
  } else {
    return(list(x=rx*step,y=z));
  }
}

# calculate window cross-correlation
lwcc <- function(x,y,s,e,whs=100,isize=20,return.peaks=T,min.thr=1,min.dist=100,step=1,tag.weight=1,bg.x=NULL,bg.y=NULL,bg.weight=NULL,mask.x=NULL,mask.y=NULL,bg.whs=whs,round.up=F) {
  if(step>1) {
    x <- floor(x/step+0.5); y <- floor(y/step+0.5)
    
    if(!is.null(bg.x)) {
      bg.x <- floor(bg.x/step+0.5); bg.y <- floor(bg.y/step+0.5)  
    }
    
    if(!is.null(mask.x)) {
      mask.x <- floor(mask.x/step+0.5); mask.y <- floor(mask.y/step+0.5)  
    }

    whs <- floor(whs/step+0.5);
    bg.whs <- floor(bg.whs/step+0.5);
    isize <- floor(isize/step+0.5);
    min.dist <- floor(min.dist/step +0.5);
    s <- floor(s/step+0.5)
    e <- floor(e/step+0.5)
  }

  # scale bg.weight, since within calculation they are considered independent
  bg.weight <- bg.weight*tag.weight;

  
  rx <- c(s-whs,e+whs);
  xt <- table(x);
  xh <- integer(diff(rx)+1);
  xh[as.integer(names(xt))-rx[1]+1] <- as.integer(xt);

  yt <- table(y);
  
  yh <- integer(diff(rx)+1);
  yh[as.integer(names(yt))-rx[1]+1] <- as.integer(yt);

  # compile background vectors
  if(!is.null(bg.x) && length(bg.x)>0) {
    bg.subtract <- 1;

    bg.xt <- table(bg.x);
    bg.xh <- integer(diff(rx)+1);
    bg.xh[as.integer(names(bg.xt))-rx[1]+1] <- as.integer(bg.xt);
    rm(bg.xt);

    bg.yt <- table(bg.y);
    bg.yh <- integer(diff(rx)+1);
    bg.yh[as.integer(names(bg.yt))-rx[1]+1] <- as.integer(bg.yt);
    rm(bg.yt);

    # adjust bg.weight according to bg.whs
    bg.weight <- bg.weight*(whs-isize)/bg.whs;
  } else {
    bg.subtract <- 0;
    bg.xh <- bg.yh <- c();
  }

  # record masked positions
  if(!is.null(mask.x) && length(mask.x)>0) {
    mvx <- unique(mask.x); mvx <- setdiff(mvx,as.numeric(names(xt)));
    mvx <- mvx[mvx>=rx[1] & mvx<=rx[2]];
    
    xh[mvx-rx[1]+1] <- -1;
  }

  if(!is.null(mask.y) && length(mask.y)>0) {
    mvy <- unique(mask.y); mvy <- setdiff(mvy,as.numeric(names(yt)));
    mvy <- mvy[mvy>=rx[1] & mvy<=rx[2]];
    yh[mvy-rx[1]+1] <- -1;
  } 
  
  rm(xt,yt);
  if(round.up) { round.up <- 1; } else { round.up <- 0; }
  
  storage.mode(xh) <- storage.mode(yh) <- "integer";
  storage.mode(bg.xh) <- storage.mode(bg.yh) <- "integer";
  nx <- length(xh);   storage.mode(nx) <- storage.mode(whs) <- storage.mode(isize) <- storage.mode(bg.whs) <- "integer";
  rp <- as.integer(return.peaks);
  storage.mode(rp) <- storage.mode(min.dist) <- "integer";
  storage.mode(min.thr) <- "double";
  storage.mode(tag.weight) <- "double";
  storage.mode(bg.weight) <- "double";
  storage.mode(bg.subtract) <- "integer";
  storage.mode(round.up) <- "integer";

  # allocate return arrays
  #cc <- numeric(nx); storage.mode(cc) <- "double";
  z <- .Call("lwcc",xh,yh,whs,isize,rp,min.dist,min.thr,tag.weight,bg.subtract,bg.xh,bg.yh,bg.whs,bg.weight,round.up);
  if(return.peaks) {
    return(data.frame(x=(z$x+rx[1])*step,y=z$v));
  } else {
    return(list(x=rx*step,y=z));
  }
}

tag.lwcc <- function(ctv,s,e,return.peaks=T, bg.ctv=NULL, mask.ctv=NULL, ...) {
  x <- ctv[ctv>=s & ctv<=e];
  y <- (-1)*ctv[ctv<=-s & ctv>=-e];

  if(!is.null(bg.ctv)) {
    bg.x <- bg.ctv[bg.ctv>=s & bg.ctv<=e];
    bg.y <- (-1)*bg.ctv[bg.ctv<=-s & bg.ctv>=-e];
  } else {
    bg.x <- bg.y <- NULL;
  }

  if(!is.null(mask.ctv)) {
    mask.x <- mask.ctv[mask.ctv>=s & mask.ctv<=e];
    mask.y <- (-1)*mask.ctv[mask.ctv<=-s & mask.ctv>=-e];
  } else {
    mask.x <- mask.y <- NULL;
  }
  
  if(length(x)==0 | length(y) ==0) {
    if(return.peaks) {
      return(data.frame(x=c(),y=c()));
    } else {
      rx <- range(c(x,y));
      return(list(x=rx,y=numeric(diff(rx)+1)));
    }
  } else { 
    return(lwcc(x,y, s,e,return.peaks=return.peaks, bg.x=bg.x,bg.y=bg.y,  mask.x=mask.x,mask.y=mask.y, ...))
  }
}

# filter predictions to remove calls failling into the tag enrichment clusters ( chr list of $s/$e dfs)
filter.binding.sites <- function(bd,tec,exclude=F) {
  chrl <- names(bd); names(chrl) <- chrl;
  lapply(chrl,function(chr) {
    cbd <- bd[[chr]];
    if(is.null(cbd)) { return(NULL) };
    if(length(cbd)==0) { return(NULL) };
    if(dim(cbd)[1]>0) {
      ctec <- tec[[chr]];
      if(length(ctec$s)>0) {
        if(exclude) {
          pwi <- which(points.within(cbd$x,ctec$s,ctec$e)== -1);
        } else {
          pwi <- which(points.within(cbd$x,ctec$s,ctec$e)> -1);
        }
        return(cbd[pwi,]);
      } else {
        if(exclude) {
          return(cbd);
        } else {
          return(data.frame(x=c(),y=c()));
        }
      }
    } else {
      return(cbd);
    }
  });  
}

# given list of magnitude values for signal(x) and control (y),
# return a dataframe with $e.val and $fdr
get.eval.fdr.vectors <- function(x,y) {
  nx <- length(x); ny <- length(y);
  if(nx==0) { return(data.frame(evalue=c(),fdr=c())) }
  if(ny==0) { return(data.frame(evalue=rep(0,nx),fdr=rep(1,nx))) }
  ex <- ecdf(x); ey <- ecdf(y);

  evals <- (1-ey(x))*ny;
  yvals <- (1-ex(x))*nx;
  fdr <- (evals+0.5)/(yvals+0.5); # with pseudo-counts
  fdr[yvals==0] <- min(fdr); # correct for undercounts
  # find a min x corresponding to a minimal FDR
  mfdr <- min(fdr);
  mfdrmx <- min(x[fdr==mfdr]);
  # correct
  fdr[x>=mfdrmx] <- mfdr;
  return(data.frame(evalue=(evals+1),fdr=fdr));
}
# estimates enrichment confidence interval based on 2*tag.count.whs window around each position, and a z-score (alpha/2)
# if(multiple.background.scales=T) the enrichment is also estimated using 5- and 10-fold increased background tag window
# adds $enr (lower bound), $enr.ub (upper bound) and $enr.mle fields
calculate.enrichment.estimates <- function(binding.positions, 
	signal.data=NULL, control.data=NULL, fraction=1, tag.count.whs=100, 
	z=2, effective.genome.size=3e9, scale.down.control=F, 
	background.scales=c(1), bg.weight=NULL) {
	
	f <- fraction
	qv <- pnorm(z,lower.tail=F)
	cn <- names(binding.positions$npl)
	names(cn) <- cn

	if(is.null(control.data)) {
		# estimate from gamma distribution
		fg.lambda <- f*sum(unlist(lapply(signal.data,length)))*2*tag.count.whs/effective.genome.size
		binding.positions$npl <- lapply(binding.positions$npl,function(d) {
			d$enr <- qgamma(qv,d$nt,scale=1)/fg.lambda
			d$enr.ub <- qgamma(1-qv,d$nt,scale=1)/fg.lambda
			d$enr.mle <- d$nt/fg.lambda
			return(d)
		})      
	} else {
		# estimate using beta distribution
		if(is.null(bg.weight)) {
			bg.weight <- sum(unlist(lapply(signal.data,length)))/sum(unlist(lapply(control.data,length)))
		}

		if(scale.down.control) {
			# sample down control to be the same size as true signal.data (bg.weight*f)
			control.data <- lapply(control.data,function(d) 
				sample(d,length(d)*bg.weight*f,replace=(f*bg.weight>1)))
			#bg.weight <- sum(unlist(lapply(signal.data,length)))/sum(unlist(lapply(control.data,length)))
			bg.weight <- 1/f
		}

		binding.positions$enrichment.bg.weight <- bg.weight
		binding.positions$enrichment.whs <- tag.count.whs
		binding.positions$enrichment.z <- z

		binding.positions$npl <- lapply(cn, function(chr) {
			d <- binding.positions$npl[[chr]]

			edf <- lapply(background.scales,function(background.width.multiplier) {
				sig.mult <- bg.weight*f/background.width.multiplier
				nbg <- points.within(abs(control.data[[chr]]), 
					d$x-tag.count.whs*background.width.multiplier, 
					d$x+tag.count.whs*background.width.multiplier,
					return.point.counts=T, return.unique=F)
				nfg <- d$nt
				# Poisson ratio Bayesian LB with non-informative prior (Clopper & Pearson 1934)
				nf <- ((nfg+0.5)/(nbg+0.5))*qf(1-qv,2*(nfg+0.5),2*(nbg+0.5),lower.tail=F)
				nf <- nf/sig.mult

				ub <- ((nfg+0.5)/(nbg+0.5))*qf(qv,2*(nfg+0.5),2*(nbg+0.5),lower.tail=F)
				ub <- ub/sig.mult

				mle <- (nfg+0.5)/(nbg+0.5)
				mle <- mle/sig.mult
				if(is.null(nbg)) { nbg <- numeric(0) }
				if(is.null(nf)) { nf <- numeric(0) }
				if(is.null(ub)) { ub <- numeric(0) }
				if(is.null(mle)) { mle <- numeric(0) }
				return(data.frame(nbg=nbg,lb=nf,ub=ub,mle=mle))
			})

			adf <- do.call(cbind,lapply(c(1:length(background.scales)),function(i) {
				df <- edf[[i]]
				cn <- c("nbgt","enr","enr.ub","enr.mle")
				if(background.scales[i]>1) {
					cn <- paste(cn,as.character(background.scales[i]),sep=".")
				}
				names(df) <- cn
				return(df)
			}))
			return(cbind(d,adf))
		})
	} 
	return(binding.positions)
}

write.probe.wig <- function(chr,pos,val,fname,append=F,feature="M",probe.length=35,header=T,name="Bed Format") {
  if(is.null(pos) | length(pos)==0) { return(); }
  #mdat <- paste(rep(chr,length(s)),as.integer(pos),as.integer(pos+probe.length),val,sep="\t")
  mdat <- data.frame(chr=chr,s=as.integer(pos),e=as.integer(pos+probe.length),v=val)
  # make sure segments don't overlap
  mdd <- mdat$s[-1]-mdat$e[-dim(mdat)[1]];
  ivi <- which(mdd<1);
  mdat$e[ivi] <- as.integer(mdat$e[ivi]+mdd[ivi]-1);
  mdat <- data.frame(chr=mdat$chr,s=as.integer(mdat$s),e=as.integer(mdat$e),v=mdat$v);

  if(header) {
    write(paste("track type=wiggle_0 name=\"",name,"\" description=\"",feature,"\" visibility=dense color=200,100,0 altColor=0,100,200 priority=20",sep=""),file=fname,append=append)
    write.table(mdat,file=fname,col.names=F,row.names=F,quote=F,sep=" ",append=T);
  } else {
    write.table(mdat,file=fname,col.names=F,row.names=F,quote=F,sep=" ",append=append);
  }

}




