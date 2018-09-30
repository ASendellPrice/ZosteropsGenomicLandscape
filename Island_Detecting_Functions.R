# Functions required for detecting genomic islands of divergence from windowed FSTs
# Based on script developed by Benjamin Van Doren

######### Functions #########

add.alpha <- function(col, alpha = 1) {
  if (missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb) / 255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha = alpha))
}




# Peak calling


fill.peaks = function(logical.vector,
                      retain.peaks.greater.than,
                      fill.gaps.less.than
                      ) {
  u = rle(logical.vector)
  # find T-F-T pattern
  if (length(u$values) > 2) {
    for (i in 2:(length(u$values) - 1)) {
      if (!NA %in% u$values[(i - 1):(i + 1)]) {
        # browser()
        if (u$values[i - 1] == TRUE &&
            u$values[i] == FALSE && u$values[i + 1] == TRUE) {
          if (u$lengths[i] < fill.gaps.less.than) {
            u$values[i] = T
          }
        }
      }
    }
  }
  # u$values[!is.na(u$values) & u$values == F & u$lengths < fill.gaps.less.than] = T
  u = rle(inverse.rle(u))
  u$values[u$values == T &
             u$lengths <= retain.peaks.greater.than] = F
  return(inverse.rle(u))
}

get.peak.info = function(peak.rle,
                         data,
                         data.col,
                         win.size,
                         type = c("high", "low"),
                         rm.false = T) {
  r = peak.rle
  rdf = data.frame(
    type = type[1],
    len = r$lengths,
    peak = r$values,
    end = cumsum(r$lengths)
  )
  rdf$beg = rdf$end - rdf$len + 1
  rdf = rdf[, c(1, 3, 5, 4, 2)]
  rdf$mean.val = apply(rdf[, 3:4], 1, function(X) {
    mean(data[X[1]:X[2], data.col])
  })
  # rdf$n.bases = apply(rdf[,3:4],1,function(X) { diff(data[c(X[2],X[1]),"mean.pos"])+win.size })
  if (rm.false) {
    rdf = rdf[rdf$peak, ]
  }
  rownames(rdf) = NULL
  if (nrow(rdf)==0) {
    rdf = data.frame(type=type,peak=NA,beg=NA,end=NA,len=NA,mean.val=NA)
  }
  return(rdf)
}

# Output number of regions, average sizes, sd of sizes, total size, 
peak.data.summary = function(pd,win.size.multiplier=50,type.return="both") { # approx. size in Kb
  pd = pd[pd$peak & !is.na(pd$peak),]
  pd$len = pd$len*win.size.multiplier
  res = by(pd,pd$pop,function(df) { 
    by(df,df$type,function(ddf) { 
      list(n.regions = nrow(ddf),median.size = median(ddf$len),max.size=max(ddf$len),
           tot.size = sum(ddf$len))
    })  
  })
  res.sum = lapply(res,function(X) { 
    d = do.call("rbind",X) 
    if (!"high" %in% rownames(d)) {
      d = rbind(rep(NA,4),d)
    } else if (!"low" %in% rownames(d)) {
      d = rbind(d,rep(NA,4))
    }
    return(d)
  })
  rownames(res.sum) = NULL
  res.out = suppressWarnings(data.frame(type=rep(c("High","Low"),length(res.sum)),pop=rep(names(res.sum),each=2),do.call("rbind",res.sum)))
  for (i in 3:6) { res.out[,i]=as.numeric(res.out[,i]) }
  colnames(res.out) = c("Type","Group","Num. Outlier Regions","Median Region Size (Kb)","Max Region Size (Kb)","Total Length (Kb)")
  if (type.return!="both") {
    type.return=toTitleCase(type.return)
    res.out=res.out[res.out$Type==type.return,]
  }
  out = rbind(res.out,c(NA,NA,apply(res.out[,3:6],2,mean)))
  return(out[complete.cases(out),])
}

peak.permutation.test = function(dataset,
                                 data.col,
                                 bandwidth,
                                 N,
                                 win.size,
                                 peak.pop,
                                 plot = T,
                                 treat.Z = c("same", "exclude", "only", "separate"),
                                 retain.peaks.greater.than = 0,
                                 fill.gaps.less.than = 5,
                                 kernel.type="box") {
  # data.col is number or name of column
  y = dataset[, data.col]
  x = dataset$index
  cs = dataset$ficed_chr
  stopifnot(!NA %in% cs)
  if (length(treat.Z) > 1) {
    treat.Z == "same"
  }
  stopifnot(treat.Z %in% c("same", "exclude", "only", "separate"))
  
  onZ = cs %in% c("chr_Z", "chr_Z_un")
  
  y.split.smooth = by(cbind(x, y), cs, function(df) {
    ys = ksmooth(
      df$x,
      df$y,
      kernel = kernel.type,
      bandwidth = bandwidth,
      x.points = df$x
    )$y
    if (nrow(df) != length(ys))
      browser()
    return(cbind(df, ys))
  })
  y.split.smooth = y.split.smooth[unique(cs)]
  
  # permutation loop
  sims = array(dim = c(length(y), N))
  for (i in 1:N) {
    y.perm = rep(NA, length(y))
    sg.y.perm = rep(NA, length(y))
    if (treat.Z == "same") {
      y.perm = sample(y, length(y), replace = F)  # permute everything equally
      sg.y.perm = ksmooth(x, y.perm, kernel = kernel.type, bandwidth = bandwidth)$y
    } else if (treat.Z == "exclude") {
      y.perm[!onZ] = sample(y[!onZ], length(y[!onZ]), replace = F)  # permute non-Z locations
      # leave windows on Z as NA
      sg.y.perm[!onZ] = ksmooth(1:length(x[!onZ]),
                                y.perm[!onZ],
                                kernel = kernel.type,
                                bandwidth = bandwidth)$y
    } else if (treat.Z == "only") {
      y.perm[onZ] = sample(y[onZ], length(y[onZ]), replace = F)  # permute Z locations
      # leave windows not on Z as NA
      sg.y.perm[onZ] = ksmooth(1:length(x[onZ]),
                               y.perm[onZ],
                               kernel = kernel.type,
                               bandwidth = bandwidth)$y
    } else if (treat.Z == "separate") {
      y.perm[!onZ] = sample(y[!onZ], length(y[!onZ]), replace = F)  # permute non-Z locations
      y.perm[onZ] = sample(y[onZ], length(y[onZ]), replace = F)  # permute Z locations
      sg.y.perm[!onZ] = ksmooth(1:length(x[!onZ]),
                                y.perm[!onZ],
                                kernel = kernel.type,
                                bandwidth = bandwidth)$y
      sg.y.perm[onZ] = ksmooth(1:length(x[onZ]),
                               y.perm[onZ],
                               kernel = kernel.type,
                               bandwidth = bandwidth)$y
    }
    sims[, i] = sg.y.perm
  }
  
  # ci = apply(sims, 1, quantile,probs=c(0.001,0.999))
  # ci = rbind(apply(sims, 1, min), apply(sims, 1, max)) # missing values will return as NA
  
  if (treat.Z %in% c("same","exclude","only")) {
    # qq = quantile(sims,probs=c(0.0001,0.9999),na.rm=T)
    # min = rep(qq[1],dim(sims)[1])
    # max = rep(qq[2],dim(sims)[1])
    min = apply(sims, 1, min)
    max = apply(sims, 1, max)
    if (treat.Z == "exclude") {
      min[onZ] = max[onZ] = NA
    } else if (treat.Z == "only") {
      min[!onZ] = max[!onZ] = NA 
    }
  } else if (treat.Z == "separate") {
    # qq.Z = quantile(sims[onZ],probs=c(0.0001,0.9999),na.rm=T)
    # qq.notZ = quantile(sims[!onZ],probs=c(0.0001,0.9999),na.rm=T)
    # min = rep(qq.notZ[1],dim(sims)[1])
    # max = rep(qq.notZ[2],dim(sims)[1])
    # min[onZ] = qq.Z[1]
    # max[onZ] = qq.Z[2]
    min = apply(sims, 1, min)
    max = apply(sims, 1, max)
  }
  ci = rbind(min,max)
  
  # browser()
  # single 
  
  
  if (plot) {
    plot(x, y, col = gray(.4), pch = ".")
    lines(ci[1, ],
          col = 3,
          type = "l",
          lwd = 2)
    lines(ci[2, ],
          col = 3,
          type = "l",
          lwd = 2)
    lapply(y.split.smooth, function(df) {
      lines(df$x, df$ys, col = 2, lwd = 2)
    })
  }
  
  y.combined = unlist(lapply(y.split.smooth, function(df) {
    df$ys
  }))
  
  stopifnot(length(y.combined) == length(y))
  
  over = y.combined > ci[2, ]
  under = y.combined < ci[1, ]
  
  # fill/call peaks separately for each chromosome
  over.per.peak = unlist(tapply(over, cs, fill.peaks, 
                                retain.peaks.greater.than = retain.peaks.greater.than,
                                fill.gaps.less.than = fill.gaps.less.than)[unique(cs)])
  under.per.peak = unlist(tapply(under, cs, fill.peaks, 
                                 retain.peaks.greater.than = retain.peaks.greater.than,
                                 fill.gaps.less.than = fill.gaps.less.than)[unique(cs)])
  
  if (plot) {
    abline(v = which(over), col = add.alpha("orange", .01))
    abline(v = which(under), col = add.alpha("blue", .01))
    
    abline(v = which(over.per.peak),
           col = add.alpha("orange", .05))
    abline(v = which(under.per.peak),
           col = add.alpha("blue", .05))
  }
  
  over.info = get.peak.info(rle(over.per.peak),
                            dataset,
                            data.col,
                            win.size,
                            "high",
                            rm.false = T)
  under.info = get.peak.info(rle(under.per.peak),
                             dataset,
                             data.col,
                             win.size,
                             "low",
                             rm.false = T)
  # browser()
  peak.info = rbind(over.info, under.info)
  peak.info$pop = peak.pop
  return(peak.info)
  
}


plot.peaks.as.lines = function(peak.data,
                               peak.pop,
                               y.low,
                               y.high,
                               col.high = add.alpha("orange", .2),
                               col.low = add.alpha("blue", .2),
                               border.opacity=0.5,
                               lwd = 2) {
  inds.high = peak.data[peak.data$pop == peak.pop &
                          peak.data$type == "high", c("beg", "end")]
  inds.low = peak.data[peak.data$pop == peak.pop &
                         peak.data$type == "low", c("beg", "end")]
  if (nrow(inds.high)>0) {
    apply(inds.high, 1, function(X) {
      # lines(c(X[1],X[2]),rep(y.low,2),lwd=lwd,lend=0,col=col.high)
      polygon(
        c(X[1], X[1], X[2], X[2]),
        c(y.low, y.high, y.high, y.low),
        col = col.high,
        border = add.alpha(col.high,border.opacity),
        lwd=.5
      )
    })
  }
  if (nrow(inds.low)>0) {
    apply(inds.low, 1, function(X) {
      # lines(c(X[1],X[2]),rep(y.low,2),lwd=lwd,lend=0,col=col.low)
      polygon(
        c(X[1], X[1], X[2], X[2]),
        c(y.low, y.high, y.high, y.low),
        col = col.low,
        border = add.alpha(col.low,border.opacity),
        lwd=.5
      )
    })
  }
}

plot.genomic.scan = function(data,
                             chrs,
                             ylim,
                             peak.pop,
                             peak.data = NULL,
                             scaff.names = F,
                             scatter = T,
                             fill.val = 0,
                             span = 0.001,
                             bw = 30,
                             lwd=2,
                             scatter.opacity = 0.2,
                             line.col = "red",
                             xlim=NA,
                             all.chr.names.above=F,
                             add=F,
                             chr.nums=T,
                             col.high = add.alpha("orange", .2),
                             col.low = add.alpha("blue", .2),
                             plot.chrs=T,
                             ...) {
  if (is.na(xlim)) {
    xlim = c(max(1, min(data$index)),
             min(chrs$end.ind[nrow(chrs)], max(data$index)))
  } else if (class(xlim)=="character") {
    data = data[data$ficed_chr %in% xlim,]
    chrs = chrs[chrs$chr %in% xlim,]
    if (!is.null(peak.data)) {
      peak.data = peak.data[peak.data$beg %in% data$index & 
                              peak.data$end %in% data$index,]
    }
    # xlim = c(chrs$start.ind[chrs$chr==xlim],chrs$end.ind[chrs$chr==xlim])
    xlim=NULL
  }
  if (is.null(ylim)) {
    ylim = c(quantile(data$val,.001),quantile(data$val,.999))
  }
  # browser()
  if (scatter) {
    if (add) {
      points(
        data$index,
        data$val,
        pch = ".",
        col = add.alpha(data$color, scatter.opacity),
        ...
      )
    } else {
      plot(
        data$index,
        data$val,
        pch = ".",
        xlab = "",
        xaxt = "n",
        xlim = xlim,
        ylim = ylim,
        col = add.alpha(data$color, scatter.opacity),
        bty = "n",
        ...
      )
    }
    
    if (!is.null(peak.data)) {
      plot.peaks.as.lines(peak.data, peak.pop, ylim[1], ylim[2], col.high = col.high, col.low = col.low)
    }
    
    y.split.smooth = by(data[, c("index", "val")], data$ficed_chr, function(df) {
      ys = ksmooth(
        df$index,
        df$val,
        kernel = "box",
        bandwidth = bw,
        x.points = df$index
      )$y
      return(cbind(df, ys))
    })
    y.split.smooth = y.split.smooth[unique(data$ficed_chr)]
    
    lapply(y.split.smooth, function(df) {
      lines(df$index, df$ys, col = line.col, lwd = lwd)
    })
    # browser()
    
    
  } else {
    inds = 1:length(data$val)
    plot(
      data$index,
      data$val,
      type = "n",
      ylim = ylim,
      bty = "n",
      xlab = "",
      xaxt = "n",
      xlim = xlim,
      ...
    )
    l = loess(data$val ~ data$index, span =
                span, degree = 0)
    lines(data$index[!is.na(data$val)], l$fitted)
    xx = data$index[!is.na(data$val)]
    yy = l$fitted
    yy[1] = yy[length(yy)] = fill.val
    polygon(xx, yy, col = "black")
  }
  
  if (!add) {
    
    if (all.chr.names.above) { chrs$up=0 }
    
    if (plot.chrs) {
      apply(chrs, 1, function(X) {
        lines(X[3:4],
              rep(ylim[2], 2),
              col = X[7],
              lwd = 3,
              lend = 1)
        rect(X[3],
             ylim[1],
             X[4],
             ylim[2],
             col = add.alpha(X[7], .06),
             border = add.alpha("black", .05))
        if (chr.nums) {
          text(
            mean(as.numeric(X[3:4])),
            ylim[2],
            adj = c(0.5, -1 + as.numeric(X[8])),
            labels = X[6],
            cex = .8
          )
        }
      })
    }
    if (scaff.names) {
      scaff.inds = tapply(data$val, data$chr, mean)
      scaff.inds = sort(scaff.inds)
      text(
        scaff.inds[c(T, F, F)],
        ylim[1] + diff(ylim) * .9,
        names(scaff.inds)[c(T, F, F)],
        cex = .2,
        adj = c(0.5, -2)
      )
      text(
        scaff.inds[c(F, T, F)],
        ylim[1] + diff(ylim) * .9,
        names(scaff.inds)[c(F, T, F)],
        cex = .2,
        adj = c(0.5, 0)
      )
      text(
        scaff.inds[c(F, F, T)],
        ylim[1] + diff(ylim) * .9,
        names(scaff.inds)[c(F, F, T)],
        cex = .2,
        adj = c(0.5, 2)
      )
    }
  }
}

add.chr.x.axis <- function(all.chrs, chr,win.size=5e4,inc.size=5) {
  chr.lims= as.numeric(all.chrs[all.chrs$chr==chr,c("start.ind","end.ind")])
  # win.size=5e4
  # inc.size = 5 # Mb
  inc.wins=1e6*inc.size/5e4
  axis(1,at=seq(from=chr.lims[1],to=chr.lims[2],by=inc.wins),
       labels=seq(from=0,to=floor(diff(chr.lims)*win.size/1e6),by=inc.size)) # MB
  mtext(text="Distance Along Draft Chromosome (Mb)",side=1,line=2.5,cex=.65)
}


plot.fst = function(data,
                    chrs,
                    data.col,
                    ylim,
                    peak.pop = NULL,
                    peak.data = NULL,
                    scaff.names = F,
                    scatter = T,
                    fill.val = 0,
                    span = 0.001,
                    bw = 30,
                    scatter.opacity = 0.2,
                    line.col = "red",
                    xlim=NA,
                    all.chr.names.above=F,
                    chr.nums=T,
                    col.high = add.alpha("orange", .2),
                    col.low = add.alpha("blue", .2),
                    ...) {
  data$val = data[, data.col]
  plot.genomic.scan(
    data = data,
    chrs = chrs,
    ylim = ylim,
    peak.pop = peak.pop,
    peak.data = peak.data,
    scaff.names = scaff.names,
    scatter = scatter,
    fill.val = fill.val,
    span = span,
    bw = bw,
    scatter.opacity = scatter.opacity,
    line.col = line.col,
    xlim = xlim,
    all.chr.names.above=all.chr.names.above,
    chr.nums=chr.nums,
    col.high = col.high,
    col.low = col.low,
    ...
  )
}







# p-value for jaccard
# null hypothesis where both pops have the same number of outlier windows of given sizes, 
# but permute their location

# re-draw interval start locations. possible draws are:
# (1 to max-interval_size) .. minus already-taken intervals (minus interval size)
resample.intervals = function(int, max) {
  # browser()
  new.intervals = NULL
  for (i in sample(1:nrow(int),nrow(int))) { # for number of intervals - random 
    new.size = size(int[i,])
    # browser()
    bad.starts = new.intervals
    bad.starts[,1] = bad.starts[,1]-new.size
    # bad.starts = Intervals_full(rbind(bad.starts,c(max-new.size,max)))
    bad.starts = rbind(bad.starts,c(max-new.size,max))
    
    # poss.starts = interval_difference(Intervals_full(c(1,max)),bad.starts)
    poss.starts = which(!(1:max) %in% expand.indices(bad.starts))
    
    # start = sample(expand.indices(poss.starts),1)
    start = sample(poss.starts,1)
    new = c(start,start+new.size)
    new.intervals = rbind(new.intervals,new)
  }
  new.intervals = Intervals(new.intervals)
  return(new.intervals)
}

simulate.null.overlap = function(int1,int2,max,N,hold.constant=NULL) { # to not resample int1, set hold.constant=1, etc.
  jaccards = rep(NA,N)
  n.overlapping = rep(NA,N)
  for (i in 1:N) {
    if (!is.null(hold.constant)) {
      if (hold.constant==1) {
        new.int1 = int1
        new.int2 = resample.intervals(int2,max)
      } else if (hold.constant==2) {
        new.int1 = resample.intervals(int1,max)
        new.int2 = int2
      } 
    } else {
      new.int1 = resample.intervals(int1,max)
      new.int2 = resample.intervals(int2,max)
    }
    new.intersection = interval_intersection(new.int1,new.int2)
    new.jaccard = sum(size(new.intersection)+1)/
      sum(size(interval_union(new.int1,new.int2))+1)
    jaccards[i] = new.jaccard
    # n.overlapping[i] = nrow(new.intersection)
    new.union = interval_union(new.int1,new.int2)
    n.overlapping[i] = sum(unlist(lapply(interval_overlap(new.union,new.int1),function(X) { length(X)>0 })) &
          unlist(lapply(interval_overlap(new.union,new.int2),function(X) { length(X)>0 })))
    # if (i%%100==0) { print(i) }
  }
  res = list(jaccards,n.overlapping)
  names(res) = c("jaccards","n.overlapping")
  return(res)
}

test.jaccard = function(int1,int2,max,N,plot=F,hold.constant=NULL) {
  # browser()
  sim = simulate.null.overlap(int1,int2,max,N,hold.constant = hold.constant)  
  obs.intersection = interval_intersection(int1,int2)
  observed.jacc = sum(size(obs.intersection)+1)/sum(size(interval_union(int1,int2))+1)
  # observed.n.overlapping = nrow(obs.intersection)
  # int1 = Intervals(int1)
  # int2 = Intervals(int2)
  union = interval_union(int1,int2)
  observed.n.overlapping = sum(unlist(lapply(interval_overlap(union,int1),function(X) { length(X)>0 })) &
        unlist(lapply(interval_overlap(union,int2),function(X) { length(X)>0 })))
  p.one.tail.jaccard = sum(sim[["jaccards"]] >= observed.jacc) / length(sim[["jaccards"]])
  p.one.tail.overlap = sum(sim[["n.overlapping"]] >= observed.n.overlapping) / length(sim[["n.overlapping"]])
  if (plot) {
    par(mfrow=c(2,1))
    hist(sim[["jaccards"]],xlim=c(0,1.1*max(observed.jacc,max(sim[["jaccards"]]))),
         main=paste(N,"simulations"))
    abline(v=observed.jacc,col="red",lwd=2)
    hist(sim[["n.overlapping"]],xlim=c(0,1.1*max(observed.n.overlapping,max(sim[["n.overlapping"]]))),
         main=paste(N,"simulations"))
    abline(v=observed.n.overlapping,col="red",lwd=2)
    par(mfrow=c(1,1))
    # browser()
  }
  res= list(p.one.tail.jaccard,p.one.tail.overlap)
  names(res) = c("jaccard","n.overlapping")
  return(res)
}

peak.overlap = function(pd, ind1, ind2, type1, type2, max.length,n.iterations=100,hold.constant=NULL) {
  # max.length is number of windows in whole genome
  if (type1 == "both") {
    t1 = c("low", "high")
  } else {
    t1 = type1
  }
  if (type2 == "both") {
    t2 = c("low", "high")
  } else {
    t2 = type2
  }
  i1 = pd[ind1 & pd$type %in% t1, c("beg", "end")]
  i2 = pd[ind2 & pd$type %in% t2, c("beg", "end")]
  i1 = i1[complete.cases(i1),]
  i2 = i2[complete.cases(i2),]
  # browser()
  if (nrow(i1)>0 & nrow(i2)>0 ) {
    i1 = Intervals(i1)
    i2 = Intervals(i2)
    intersection = interval_intersection(i1, i2)
    # n.overlapping.intervals = nrow(intersection) # number of overlapping peaks/valleys
    union = interval_union(i1, i2)
    shared.1 = unlist(lapply(interval_overlap(union,i1),function(X) { length(X)>0 }))
    shared.2 = unlist(lapply(interval_overlap(union,i2),function(X) { length(X)>0 }))
    # print (shared.1)
    # print(shared.2)
    # browser()
    n.overlapping.intervals = sum(shared.1 & shared.2)
    
    n.1.in.2 = sum(unlist(lapply(interval_overlap(i1,i2),function(X) { length(X)>0 })))
    n.2.in.1 = sum(unlist(lapply(interval_overlap(i2,i1),function(X) { length(X)>0 })))
    
    n.unique.intervals = nrow(union)
    n.windows.shared = sum(size(intersection) + 1) # number of windows in the overlap
    n.windows.both = sum(size(union) + 1) # number of windows in the overlap
    jaccard = n.windows.shared / n.windows.both
    if (n.iterations>1) {
      # browser()
      # set.seed(12345)
      set.seed(5)
      ps = test.jaccard(i1,i2,max.length,n.iterations,hold.constant = hold.constant)  
    } else {
      ps= list(NA,NA)
      names(ps) = c("jaccard","n.overlapping") 
    }
    

    res = c(
     # pop1,
     # pop2,
      type1,
      type2,
      n.1.in.2,
      nrow(i1),
      n.2.in.1,
      nrow(i2),
      n.overlapping.intervals,
      n.unique.intervals,
      n.overlapping.intervals/n.unique.intervals,
      ps[["n.overlapping"]],
      n.windows.shared,
      sum(size(i1) + 1),
      sum(size(i2) + 1),
      n.windows.both,
      jaccard,
      ps[["jaccard"]]
    )
  } else {
    res = c(type1,type2,rep(NA,14))
  } 
  names(res) = c(
  #  "pop1",
  #  "pop2",
    "type1",
    "type2",
    "n.shared.regions.1",
    "n.total.regions.1",
    "n.shared.regions.2",
    "n.total.regions.2",
    "n.shared.regions",
    "n.unique.regions",
    "prop.shared.regions",
    "p.overlapping.regions",
    "n.shared.wins",
    "n.total.wins.1",
    "n.total.wins.2",
    "n.total.wins.both",
    "jaccard",
    "p.jaccard"
  )
  return(res)
}

get.peak.overlap.summary.fst = function(peak.data,max.length,n.iterations=100,type="fst",all.combos=F,
                                        this.one.with.others=NULL,
                                        these.combos=NULL,
                                        hold.constant=NULL,
                                        independent.only=T) {
  peak.data = peak.data[complete.cases(peak.data),]
  comp.data = NULL
  if (!is.null(this.one.with.others)) {
    other.pop = this.one.with.others
    # pops = unique(unlist(indep.comp.names))
    # pops = pops[pops != other.pop]
    # indep.comp.names = lapply(as.list(pops),function(X) { return(c(X,other.pop)) })
    pops = unique(substring(peak.data$pop,5))
    pops = pops[pops!=other.pop]
    indep.comp.names = lapply(as.list(pops),function(X) { return(c(X,other.pop)) })
  }
  if (!independent.only) {
    pops = unique(substring(peak.data$pop,5))
    all.comp.names = lapply(apply(combn(pops,2),2,as.list),unlist)
    indep.comp.names = all.comp.names
  }
  for (comp in indep.comp.names) {
    comp = paste(type,comp,sep='.')
    print(paste("comp",comp,"of",length(indep.comp.names)))
    combos = list(c("high","high"),c("low","low"),c("both","both"))
    if (all.combos) { combos = c(combos,list(c("high","low"),c("low","high"))) }
    if (!is.null(these.combos)) { combos=these.combos }
    r= NULL
    for (combo in combos) {
      r.new = peak.overlap(peak.data, peak.data$pop==comp[1], peak.data$pop==comp[2],
                     combo[1], combo[2],max.length,n.iterations,hold.constant = hold.constant)
      r = data.frame(rbind(r,r.new),stringsAsFactors = F)
    }
    r$pop1 = comp[1]
    r$pop2 = comp[2]
    comp.data = rbind(comp.data, r)
  }
  comp.data = data.frame(comp.data, stringsAsFactors = F)
  comp.data = comp.data[order(comp.data$type1, comp.data$type2), ]
  comp.data$jaccard = as.numeric(comp.data$jaccard)
  comp.data$prop.shared.regions = as.numeric(comp.data$prop.shared.regions)
  return(comp.data)
}

get.peak.overlap.summary.single = function(peak.data,max.length,n.iterations=100,these.combos=NULL) {
  peak.data = peak.data[complete.cases(peak.data),]
  comp.data = NULL
  combs = combn(unique(peak.data$pop), 2)
  for (i in 1:ncol(combs)) {
    combos = list(c("high","high"),c("low","low"),c("both","both"))
    if (!is.null(these.combos)) { combos=these.combos }
    r= NULL
    for (combo in combos) {
      r.new = peak.overlap(peak.data, peak.data$pop==combs[1, i], peak.data$pop==combs[2, i],
                           combo[1], combo[2],max.length,n.iterations)
      r = data.frame(rbind(r,r.new),stringsAsFactors = F)
    }
    r$pop1 = combs[1, i]
    r$pop2 = combs[2, i]
    comp.data = rbind(comp.data, r)
    
    # # browser()
    # print(paste("combination",i,"of",ncol(combs)))
    # r2 = peak.overlap(peak.data, peak.data$pop==combs[1, i], peak.data$pop==combs[2, i], 
    #                   "low", "low",max.length,n.iterations)
    # if (!only.low) {
    #   r1 = peak.overlap(peak.data, peak.data$pop==combs[1, i], peak.data$pop==combs[2, i], 
    #                     "high", "high",max.length,n.iterations)
    #   r3 = peak.overlap(peak.data, peak.data$pop==combs[1, i], peak.data$pop==combs[2, i], 
    #                     "both", "both",max.length,n.iterations)
    #   r = data.frame(rbind(r1,r2,r3),stringsAsFactors = F)
    # } else { r = rbind(r2) }
    # r$pop1 = combs[1, i]
    # r$pop2 = combs[2, i]
    # comp.data = rbind(comp.data, r)
  }
  comp.data = data.frame(comp.data, stringsAsFactors = F)
  comp.data = comp.data[order(comp.data$type1, comp.data$type2), ]
  comp.data$jaccard = as.numeric(comp.data$jaccard)
  comp.data$prop.shared.regions = as.numeric(comp.data$prop.shared.regions)
  return(comp.data)
}

get.peak.overlap.summary.2.metrics.by.pop = function(peak.data1,peak.data2,max.length,
                                                     n.iterations=100,all.combos=F,these.combos=NULL) {
  comp.data = NULL
  peak.data = rbind(peak.data1,peak.data2)
  peak.data$data.type = c(rep(1,nrow(peak.data1)),rep(2,nrow(peak.data2)))
  for (pop in unique(peak.data$pop)) {
    print(paste("pop",pop,"of",length(unique(peak.data$pop))))
    pd.one = peak.data[peak.data$pop==pop,]
    combos = list(c("high","high"),c("low","low"),c("both","both"))
    if (all.combos) { combos = c(combos,list(c("high","low"),c("low","high"))) }
    if (!is.null(these.combos)) { combos=these.combos }
    r= NULL
    for (combo in combos) {
      r.new = peak.overlap(pd.one, pd.one$data.type==1, pd.one$data.type==2,
                           combo[1], combo[2],max.length,n.iterations)
      r = data.frame(rbind(r,r.new),stringsAsFactors = F)
    }
    # browser()
    r$pop = pop
    comp.data = rbind(comp.data, r)
  }
  comp.data = data.frame(comp.data, stringsAsFactors = F)
  comp.data = comp.data[order(comp.data$type1, comp.data$type2), ]
  comp.data$jaccard = as.numeric(comp.data$jaccard)
  comp.data$prop.shared.regions = as.numeric(comp.data$prop.shared.regions)
  return(comp.data)
}

get.peak.overlap.summary.pop.v.fst = function(peak.popgen.data,peak.fst.type.data,
                                              max.length,n.iterations=100,all.combos=F,
                                              these.combos=NULL) {
  pfst = peak.fst.type.data
  pd = peak.popgen.data
  
  # pfst$pop.comb = pfst$pop
  p1 = regmatches(pfst$pop,regexpr("\\..+:",pfst$pop))
  p2 = regmatches(pfst$pop,regexpr(":.*",pfst$pop))
  pfst$pop1 = substring(p1,2,nchar(p1)-1)
  pfst$pop2 = substring(p2,2)
  
  peak.res = NULL
  for (comp in unique(pfst$pop)) {
    pf = pfst[pfst$pop==comp,c("type","peak","beg","end","len","mean.val","pop1")] # POP1
    colnames(pf)[length(colnames(pf))] = "pop"
    pk1 = get.peak.overlap.summary.2.metrics.by.pop(pd,pf,max.length,n.iterations=n.iterations,
                                                    all.combos = all.combos,these.combos = these.combos)  
    
    pf = pfst[pfst$pop==comp,c("type","peak","beg","end","len","mean.val","pop2")] # POP2
    colnames(pf)[length(colnames(pf))] = "pop"
    pk2 = get.peak.overlap.summary.2.metrics.by.pop(pd,pf,max.length,n.iterations=n.iterations,
                                                    all.combos = all.combos,these.combos = these.combos)  
    
    pk = rbind(pk1,pk2)
    pk = pk[!is.na(pk$n.shared.regions),]
    pk$comp = comp
    peak.res = rbind(peak.res,pk)
  }
  rownames(peak.res) = NULL
  return(peak.res)
}

peak.overlap.to.dist = function(peak.overlap.summary, type1, type2) {
  s = peak.overlap.summary
  s = s[s$type1 == type1 & s$type2 == type2, ]
  
  pops = unique(c(s$pop1, s$pop2))
  dist = diag(rep(0, length(pops)))
  dimnames(dist) = list(pops, pops)
  
  for (i in 1:nrow(s)) {
    dist[s$pop1[i], s$pop2[i]] = dist[s$pop2[i], s$pop1[i]] = 1 - as.numeric(s$jaccard[i])
  }
  return(dist)
}

f.dec = function(num,n.dec=2) {
  format(round(num,n.dec),nsmall=n.dec)
}

overlap.table = function(over.data,type1,type2,other.columns=NULL,short.summary=F) {
  o = over.data
  # browser()
  o = o[o$type1 %in% type1 & o$type2 %in% type2,]
  o = o[rev(order(o$prop.shared.regions)),]
  o$prop.shared.regions = round(	  o$prop.shared.regions,2)
  o$prop.shared.1 = round(as.numeric(o$n.shared.regions.1)/as.numeric(o$n.total.regions.1),2)
  o$prop.shared.2 = round(as.numeric(o$n.shared.regions.2)/as.numeric(o$n.total.regions.2),2)
  o$jaccard = round(o$jaccard,2)	  
  # browser()
  if (short.summary) {
    # number (percent) of comparisons that are significant at alpha=0.05
    # mean prop shared regions across all comparisons (sig comparisons)
    # mean prop shared regions 1 across all comparisons (sig comparisons)
    # mean prop shared regions 2 across all comparisons (sig comparisons)
    sig.ind.regions = o$p.overlapping.regions<=0.05
    sig.ind.wins = o$p.jaccard<=0.05
    n.sig.regions = sum(sig.ind.regions,na.rm=T)
    n.sig.wins = sum(sig.ind.wins,na.rm=T)
    mean.prop.shared = mean(o$prop.shared.regions,na.rm=T)
    max.prop.shared = max(o$prop.shared.regions,na.rm=T)
    min.prop.shared = min(o$prop.shared.regions,na.rm=T)
    # mean.prop.shared.sig = mean(o$prop.shared.regions[sig.ind],na.rm=T)
    mean.prop.shared.1 = mean(o$prop.shared.1,na.rm=T)
    max.prop.shared.1 = max(o$prop.shared.1,na.rm=T)
    min.prop.shared.1 = min(o$prop.shared.1,na.rm=T)
    # mean.prop.shared.1.sig = mean(o$prop.shared.1[sig.ind],na.rm=T)
    mean.prop.shared.2 = mean(o$prop.shared.2,na.rm=T)
    max.prop.shared.2 = max(o$prop.shared.2,na.rm=T)
    min.prop.shared.2 = min(o$prop.shared.2,na.rm=T)
    # mean.prop.shared.2.sig = mean(o$prop.shared.2[sig.ind],na.rm=T)
    mean.prop.shared.both = mean(c(o$prop.shared.1,o$prop.shared.2),na.rm=T)
    max.prop.shared.both = max(c(o$prop.shared.1,o$prop.shared.2),na.rm=T)
    min.prop.shared.both = min(c(o$prop.shared.1,o$prop.shared.2),na.rm=T)
    res = c(paste(n.sig.wins,"/",nrow(o),sep=""),
            paste(n.sig.regions,"/",nrow(o),sep=""),
          paste(f.dec(mean.prop.shared)," (",f.dec(min.prop.shared),"-",f.dec(max.prop.shared),")",sep=""),
          paste(f.dec(mean.prop.shared.1)," (",f.dec(min.prop.shared.1),"-",f.dec(max.prop.shared.1),")",sep=""),
          paste(f.dec(mean.prop.shared.2)," (",f.dec(min.prop.shared.2),"-",f.dec(max.prop.shared.2),")",sep=""),
          paste(f.dec(mean.prop.shared.both)," (",f.dec(min.prop.shared.both),"-",f.dec(max.prop.shared.both),")",sep=""))
    names(res) = c("Significant Tests: 50 Kb Windows",
                   "Significant Tests: Outlier Regions",
                   "Prop. Overlap",
                   "Prop. Shared of 1",
                   "Prop. Shared of 2",
                   "Prop. Shared of 1 and 2")
    return(res)
  }
  if ("pop" %in% colnames(over.data)) {
    o = o[,c("pop","n.unique.regions","prop.shared.regions","prop.shared.1","prop.shared.2","p.overlapping.regions",
             "n.total.wins.both","jaccard","p.jaccard",other.columns)]
    colnames(o) = c("Pop","N. Unique Outlier Regions",
                    "Prop. Outlier Region Overlap","Prop. Outlier Regions Shared 1","Prop. Outlier Regions Shared 2",
                    "P-value","N. Unique Outlier Windows","Prop. Outlier Windows Shared","P-value",other.columns)	  			
  } else {
    o = o[,c("pop1","pop2","n.unique.regions","prop.shared.regions","prop.shared.1","prop.shared.2","p.overlapping.regions",
             "n.total.wins.both","jaccard","p.jaccard",other.columns)]
    colnames(o) = c("Pop 1","Pop 2","N. Unique Outlier Regions",
                    "Prop. Outlier Region Overlap","Prop. Outlier Regions Shared 1","Prop. Outlier Regions Shared 2",
                    "P-value","N. Unique Outlier Windows","Prop. Outlier Windows Shared","P-value",other.columns)	  			
  }
  
  rownames(o) = NULL
  return(o)
}


# PBS #

pbs.avg.3 = function(t_ab, t_ac, t_bc) {
  return((t_ab + t_ac - t_bc) / 2)
}

fst2t = function(fst) {
  return(-log10(1 - fst)) # base10 or base e?
}

generate.pbs.comparisons = function(a, b, c, d, cols, data, transform =
                                      F) {
  # B is closer to A than to C/D
  ab = c(paste(a, b, sep = ":"), paste(b, a, sep = ":"))
  ab = ab[ab %in% cols]
  ac = c(paste(a, c, sep = ":"), paste(c, a, sep = ":"))
  ac = ac[ac %in% cols]
  bc = c(paste(b, c, sep = ":"), paste(c, b, sep = ":"))
  bc = bc[bc %in% cols]
  ad = c(paste(a, d, sep = ":"), paste(d, a, sep = ":"))
  ad = ad[ad %in% cols]
  bd = c(paste(b, d, sep = ":"), paste(d, b, sep = ":"))
  bd = bd[bd %in% cols]
  return(pbs.for.4(data[, ab], data[, ac], data[, bc], data[, ad], data[, bd], transform =
                     transform))
}

generate.pbs.comparisons.3 = function(a, b, c, cols, data, transform =
                                      F) {
  # B is closer to A than to C/D
  ab = c(paste(a, b, sep = ":"), paste(b, a, sep = ":"))
  ab = ab[ab %in% cols]
  ac = c(paste(a, c, sep = ":"), paste(c, a, sep = ":"))
  ac = ac[ac %in% cols]
  bc = c(paste(b, c, sep = ":"), paste(c, b, sep = ":"))
  bc = bc[bc %in% cols]
  return(pbs.for.3(data[, ab], data[, ac], data[, bc], transform =
                     transform))
}

pbs.for.3 = function(fst_ab,
                     fst_ac,
                     fst_bc,
                     transform = F) {
  if (transform) {
    # (pbs.avg.3(fst2t(fst_ab), fst2t(fst_ac), fst2t(fst_bc)) + pbs.avg.3(fst2t(fst_ab), fst2t(fst_ad), fst2t(fst_bd))) /
    #  2
    pbs.avg.3(fst2t(fst_ab), fst2t(fst_ac), fst2t(fst_bc))
  } else {
    # (pbs.avg.3((fst_ab), (fst_ac), (fst_bc)) + pbs.avg.3((fst_ab), (fst_ad), (fst_bd))) /
    #  2
    pbs.avg.3((fst_ab), (fst_ac), (fst_bc))
  }
}

pbs.for.4 = function(fst_ab,
                     fst_ac,
                     fst_bc,
                     fst_ad,
                     fst_bd,
                     transform = F) {
  if (transform) {
    # (pbs.avg.3(fst2t(fst_ab), fst2t(fst_ac), fst2t(fst_bc)) + pbs.avg.3(fst2t(fst_ab), fst2t(fst_ad), fst2t(fst_bd))) /
    #  2
    pbs.avg.3(fst2t(fst_ab), fst2t(fst_ac), fst2t(fst_bc)) - pbs.avg.3(fst2t(fst_ab), fst2t(fst_ad), fst2t(fst_bd))
  } else {
    # (pbs.avg.3((fst_ab), (fst_ac), (fst_bc)) + pbs.avg.3((fst_ab), (fst_ad), (fst_bd))) /
    #  2
    pbs.avg.3((fst_ab), (fst_ac), (fst_bc)) - pbs.avg.3((fst_ab), (fst_ad), (fst_bd))
  }
}

pbs.for.5 = function(fst_ab,
                     fst_ac,
                     fst_bc,
                     fst_ad,
                     fst_bd,
                     fst_ae,
                     fst_be) {
  (
    pbs.avg.3(fst2t(fst_ab), fst2t(fst_ac), fst2t(fst_bc)) +
      pbs.avg.3(fst2t(fst_ab), fst2t(fst_ad), fst2t(fst_bd)) +
      pbs.avg.3(fst2t(fst_ab), fst2t(fst_ae), fst2t(fst_be))
  ) / 3
}

# END PBS #


# FIXED SNP DENSITIES #

prepare.fixed.snp.file = function(data, pop.names) {
  f = data
  data.cols = which(grepl("=", f[1,]))
  pop.comps = substring(regmatches(f[1, data.cols], gregexpr(".*=", f[1, data.cols])), 1, 3)
  for (i in 1:length(pop.names)) {
    pop.comps = gsub(i, pop.names[i], pop.comps)
  }
  for (i in data.cols) {
    f[, i] = as.numeric(substring(f[, i], 5))
  }
  colnames(f) = c("chr", "pos", "n.snps", "cov.frac", "min.cov", pop.comps)
  f = f[,-which(colnames(f) %in% c("n.snps", "cov.frac"))]
  f$scaffold.id = gsub("Sc(0)*", "Sc", f$chr)
  f$ficed_chr = synteny$chr[match(f$scaffold.id, synteny$sca)]
  return(f)
}

prepare.fixed.snp.density = function(data, win.len, outline,prop.of.snps.that.are.fixed=F) {
  # calculate fixed snp density following windowing in outline (dxy.data or fst.data)
  # outline must have the same winlen as outline
  stopifnot(win.len == outline$mean.pos[2] - outline$mean.pos[1])
  f = data
  
  # scaffolds = unique(f$chr)
  # scaffolds = as.character(unique(outline$chr))
  data.cols = which(grepl(":", colnames(f)))
  
  # for each row in outline, create new
  outline.data.cols = which(grepl(":", colnames(outline)))
  outline$chr = as.character(outline$chr)
  outline$scaffold.id = as.character(outline$scaffold.id)
  res = outline
  res[, outline.data.cols] = NA
  
  scaffolds = unique(res$chr)
  current.scaffold = NULL
  
  print("main loop...")
  tic = proc.time()
  for (i in 1:nrow(res)) {
    # f.sc = f[scaffold.indices[[outline$chr[i]]],]
    if (is.null(current.scaffold) ||
        res$chr[i] != current.scaffold) {
      current.scaffold = res$chr[i]
      f.sc = f[f$scaffold.id == current.scaffold,]
    }
    
    start = res$mean.pos[i] - (win.len / 2)
    end = res$mean.pos[i] + (win.len / 2)
    
    stopifnot(f.sc$scaffold.id[1] == current.scaffold &&
                f.sc$scaffold.id[nrow(f.sc)] == current.scaffold)
    
    f.rel = f.sc[f.sc$pos > start & f.sc$pos <= end,]
    
    # n.fixed = apply(f.rel[, data.cols], 2, function(X) {
    #   sum(X == 1)
    # })
    
    # n.covered.sites = round(win.len * res$cov.frac[i])

    # mean.fst = apply(f.rel[, data.cols], 2, function(X) {
    #   mean(X)
    # })
    
    if (prop.of.snps.that.are.fixed) {
      # n.greater.than.zero = apply(f.rel[, data.cols], 2, function(X) {
      #   sum(X>0)
      # })
      # n.middle.50 = apply(f.rel[, data.cols], 2, function(X) {
      #   sum(X>0.25 & X<=0.75)
      # })
      mean.fst.greater.than.zero = apply(f.rel[, data.cols], 2, function(X) {
        mean(X[X>0])
      })
      # res[i, outline.data.cols] = n.fixed/n.greater.than.zero
      # res[i, outline.data.cols] = n.middle.50 / n.greater.than.zero
      res[i, outline.data.cols] = mean.fst.greater.than.zero
    } else {
      fixed.density  = n.fixed / n.covered.sites
      stopifnot(length(fixed.density) > 0)
      res[i, outline.data.cols] = fixed.density  
    }
    
    # browser()
    
    if (i %% 100 == 0) {
      elapsed.so.far = (proc.time() - tic)[3]
      to.go = elapsed.so.far / i * (nrow(res) - i)
      print(
        paste(
          round(elapsed.so.far, 1),
          "s elapsed for",
          i,
          "rows.",
          "est.",
          round(to.go, 1),
          "s to go, or",
          round(elapsed.so.far + to.go, 1),
          "sec total."
        )
      )
    }
  }
  toc = proc.time() - tic
  print(paste("done.", round(toc[3], 1), "sec."))
  
  return(res)
}



coding.sequence.density = function(data, replace.col, win.len, outline) {
  # calculate fixed snp density following windowing in outline (dxy.data or fst.data)
  # outline must have the same winlen as outline
  stopifnot(win.len == outline$mean.pos[2] - outline$mean.pos[1])
  f = data
  
  # for each row in outline, create new
  outline.data.col = replace.col
  outline$chr = as.character(outline$chr)
  res = outline
  res[, outline.data.col] = NA
  
  scaffolds = unique(res$chr)
  current.scaffold = NULL
  
  print("main loop...")
  tic = proc.time()
  for (i in 1:nrow(res)) {
    if (is.null(current.scaffold) ||
        res$chr[i] != current.scaffold) {
      current.scaffold = res$chr[i]
      f.sc = f[f$chr == current.scaffold,]
    }
    
    if (nrow(f.sc)>0) {
    
      start = res$mean.pos[i] - (win.len / 2)
      end = res$mean.pos[i] + (win.len / 2)
      
      stopifnot(f.sc$chr[1] == current.scaffold &&
                  f.sc$chr[nrow(f.sc)] == current.scaffold)
      
      f.rel = f.sc[f.sc$end > start & f.sc$start <= end,]
      
      if (nrow(f.rel)>0) {
        coding.length = sum(size(interval_intersection(Intervals(f.rel[,c("start","end")]),Intervals(c(start,end)))))  
        prop.coding = coding.length/win.len
      } else {
        prop.coding=0
      }
    } else {
      prop.coding=0
    }
    
    res[i, outline.data.col] = prop.coding
    
    if (i %% 100 == 0) {
      elapsed.so.far = (proc.time() - tic)[3]
      to.go = elapsed.so.far / i * (nrow(res) - i)
      print(
        paste(
          round(elapsed.so.far, 1),
          "s elapsed for",
          i,
          "rows.",
          "est.",
          round(to.go, 1),
          "s to go, or",
          round(elapsed.so.far + to.go, 1),
          "sec total."
        )
      )
    }
  }
  toc = proc.time() - tic
  print(paste("done.", round(toc[3], 1), "sec."))
  
  return(res)
}



# DISTANCE MATRIX AND TREE #

dist.matrix = function(data.melt, taxa.names=c("Albicollis","Austria", "Irish", "Kenya", "Siberia","Canary"), method = c("mean", "median")) {
  if (method == "mean") {
    m = with(data.melt, tapply(value, variable, mean, na.rm = T))
  } else if (method == "median") {
    m = with(data.melt, tapply(value, variable, median, na.rm = T))
  } else if (class(method) == "numeric") {
    m = with(data.melt,
             tapply(
               value,
               variable,
               quantile,
               probs = method,
               na.rm = T
             ))
  } else {
    stop("Unknown method")
  }
  
  dist = diag(rep(0, length(taxa.names)))
  dimnames(dist) = list(
    taxa.names,
    taxa.names
  )
  # browser()
  comp = data.frame(p1 = unlist(lapply(strsplit(names(
    m
  ), ":"), function(X) {
    return(X[1])
  })),
  p2 = unlist(lapply(strsplit(names(
    m
  ), ":"), function(X) {
    return(X[2])
  })), m)
  rownames(comp) = NULL
  
  # comp = comp[comp$p1!="Canary" & comp$p2!="Canary",]
  
  for (i in 1:nrow(comp)) {
    p1 = as.character(comp$p1[i])
    p2 = as.character(comp$p2[i])
    dist[eval(p1), eval(p2)] = comp$m[i]
    dist[eval(p2), eval(p1)] = comp$m[i]
  }
  return(dist)
}


# OTHER #

expand.indices = function(df) {
  if (nrow(df)>1) {
    df.expanded = as.list(apply(df,1,function(X) { return(X[1]:X[length(X)]) }))
    if (class(df.expanded)=="matrix") {
      df.expanded = lapply(apply(df.expanded,2,list),unlist)
    }
    df.cat = do.call("c",df.expanded)
  } else {
    df.cat = df[1,1]:df[1,2]
  }
  return(df.cat)
}

blast.result.to.window = function(win.scafs,win.start,win.end,s.scafs,s.start,s.end) {
  stopifnot(length(win.start)==mean(c(length(win.end),length(win.scafs))))
  stopifnot(length(s.scafs)==mean(c(length(s.start),length(s.end))))
  bins = rep(NA,length(s.start))
  for (i in 1:length(s.scafs)) {
    # browser()
    bin = which(win.scafs==s.scafs[i] & 
                  ((win.start<as.numeric(s.start[i]) & win.end>as.numeric(s.start[i])) | 
                     (win.start<as.numeric(s.end[i]) & win.end>as.numeric(s.end[i]))))  
    if (length(bin)==0) { bin=NA }
    bins[i]=bin
  }
  return(bins)
}

biplot.pcoa = function (x, Y = NULL, plot.axes = c(1, 2), dir.axis1 = 1, dir.axis2 = 1, 
          rn = NULL, main = "PCoA ordination", text.col="black", ...) 
{
  k <- ncol(x$vectors)
  if (k < 2) 
    stop("There is a single eigenvalue. No plot can be produced.")
  if (k < plot.axes[1]) 
    stop("Axis", plot.axes[1], "does not exist.")
  if (k < plot.axes[2]) 
    stop("Axis", plot.axes[2], "does not exist.")
  if (!is.null(rn)) 
    rownames(x$vectors) <- rn
  labels = colnames(x$vectors[, plot.axes])
  diag.dir <- diag(c(dir.axis1, dir.axis2))
  x$vectors[, plot.axes] <- x$vectors[, plot.axes] %*% diag.dir
  if (is.null(Y)) {
    limits <- apply(x$vectors[, plot.axes], 2, range)
    ran.x <- limits[2, 1] - limits[1, 1]
    ran.y <- limits[2, 2] - limits[1, 2]
    xlim <- c((limits[1, 1] - ran.x/10), (limits[2, 1] + 
                                            ran.x/5))
    ylim <- c((limits[1, 2] - ran.y/10), (limits[2, 2] + 
                                            ran.y/10))
    # par(mai = c(1, 1, 1, 0.5))
    plot(x$vectors[, plot.axes], xlab = labels[1], ylab = labels[2], 
         xlim = xlim, ylim = ylim, asp = 1,main=main,...)
    text(x$vectors[, plot.axes], labels = rownames(x$vectors), col=text.col,
         pos = 4, cex = 1, offset = 0.5)
    # title(main = "PCoA ordination", line = 2.5)
  }
  else {
    n <- nrow(Y)
    points.stand <- scale(x$vectors[, plot.axes])
    S <- cov(Y, points.stand)
    U <- S %*% diag((x$values$Eigenvalues[plot.axes]/(n - 
                                                        1))^(-0.5))
    colnames(U) <- colnames(x$vectors[, plot.axes])
    # par(mai = c(1, 0.5, 1.4, 0))
    biplot(x$vectors[, plot.axes], U, xlab = labels[1], ylab = labels[2],...)
    title(main = c("PCoA biplot", "Response variables projected", 
                   "as in PCA with scaling 1"), line = 4)
  }
  invisible()
}

# peak.data = peak.fst.data
# n.pops = 3
# types="high"
# return peaks shared by at least n.pops comparisons
get.overlapping.peaks = function(peak.data,n.pops,types="high",end.pop=NULL) {
  peak.data = peak.data[complete.cases(peak.data),]
  peak.data = peak.data[peak.data$type %in% types,]
  peak.data$n.overlap = NA
  all.intervals = Intervals(peak.data[,c("beg","end")])
  for (i in 1:nrow(peak.data)) {
    the.int = Intervals(peak.data[i,c("beg","end")])
    n.overlap = sum(unlist(lapply(interval_overlap(all.intervals,the.int),function(X) { length(X)>0 })))
    peak.data$n.overlap[i] = n.overlap
  }
  res.intervals = interval_union(Intervals(peak.data[peak.data$n.overlap>=n.pops,c("beg","end")]))
  if (is.null(end.pop)) {
    end.pop = peak.data$pop[1]
  } 
  res = data.frame(type=types[1],peak=T,beg=res.intervals[,1],end=res.intervals[,2],pop=end.pop)
  return(res)
}

get.overlapping.peaks.intersection = function(peak.data1,peak.data2,types="high",end.pop=NULL) {
  # browser()
  peak.data1 = peak.data1[peak.data1$type %in% types,]
  peak.data2 = peak.data2[peak.data2$type %in% types,]
  res.intervals = interval_intersection(Intervals(peak.data1[,c("beg","end")]),Intervals(peak.data2[,c("beg","end")]))
  if (is.null(end.pop)) {
    end.pop = peak.data1$pop[1]
  } 
  res = data.frame(type=types[1],peak=T,beg=res.intervals[,1],end=res.intervals[,2],pop=end.pop)
  return(res)
}

get.overlapping.peaks.1.in.2 = function(peak.data1,peak.data2,types="high",end.pop=NULL) {
  peak.data1 = peak.data1[peak.data1$type %in% types,]
  peak.data2 = peak.data2[peak.data2$type %in% types,]
  peaks.2.in.1 = unlist(lapply(interval_overlap(Intervals(peak.data1[,c("beg","end")]),Intervals(peak.data2[,c("beg","end")])),function(X) { length(X)>0 }))
  res.intervals = peak.data1[peaks.2.in.1,]
  return(res.intervals)
}

get.genes.in.interval = function(scaf,beg,end) {
  scaf = gsub("Sc(0)*", "Sc", scaf)
  if (scaf %in% gff.genes$scaf) {
    genes.in.scaf = gff.genes[gff.genes$scaf==scaf,]
    gene.intervals = Intervals(genes.in.scaf[,c("start","end")])
    overlapping.genes = which(unlist(lapply(interval_overlap(gene.intervals,Intervals(c(beg,end))),function(X) { length(X)>0 })))
    return(unique(genes.in.scaf$gene.id[overlapping.genes]))  
  }
  return(character(0))
}

get.genes.with.indices = function(indices,df) {
  all.genes = c()
  for (i in indices) {
    row = df[i,]
    genes = get.genes.in.interval(row$scaffold.id,row$mean.pos-2.5e4,row$mean.pos+2.5e4-1)  
    all.genes = c(all.genes,genes)
  }
  return(all.genes)
}

#ref = read.dna("/Users/Benjamin/Documents/Cornell/Senior_Thesis/stonechats/stonechat_genome/FINAL_ASSEMBLY_10JUNE2015/stonechat.haplomerged.v2.fa",format="fasta")
get.sequence.from.genome = function(scaf.name,mean.pos,win.size=5e4) {
  ref.names = substring(names(ref),1,-1+unlist(gregexpr(" ",names(ref))))
  
  ref.names = paste(substring(ref.names,1,unlist(gregexpr("Sc",ref.names))+1),
                    as.numeric(substring(ref.names,unlist(gregexpr("Sc",ref.names))+2)),sep="")
  if (scaf.name=="Sc") { scaf.name="Sc0" }
  if (scaf.name=="xfSc") { scaf.name="xfSc0" }
  pos = mean.pos
  scaf.ind = which(ref.names==scaf.name)
  bot = max(0,pos-win.size/2)
  top = min(pos+win.size/2,length(ref[scaf.ind][[1]]))
  scaf = as.character(ref[scaf.ind])[[1]][bot:top]
  scaf = as.DNAbin(list(scaf))
  # if (output.text) { 
  #   return(paste(as.character(s[[1]]),collapse=""))
  # }
  # return(scaf)
  nm= paste("Saxicola_",scaf.name,"_",bot,"-",top,sep="")
  nm=paste("~/Downloads/",nm,".txt",sep="")
  write.dna(scaf,nm,format="fasta")
  system(paste("open",nm))
            

}

cor.mtest.2 <- function(mat1, mat2, conf.level = 0.95,method="pearson"){
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n1, n2)
  # diag(p.mat) <- 0
  # diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:n1){
    for(j in 1:n2){
      tmp <- cor.test(mat1[,i], mat2[,j], conf.level = conf.level,method=method)
      # p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      p.mat[i,j] <- tmp$p.value
      if (method!="spearman") {
        lowCI.mat[i,j] <- tmp$conf.int[1]
        uppCI.mat[i,j] <- tmp$conf.int[2]
      }
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

cor.mtest <- function(mat, conf.level = 0.95,method="pearson"){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level,method=method)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      if (method!="spearman") {
        lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
        uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
      }
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

# Corrplot related

pairs.to.matrix = function(names1,names2,vals) {
  names = unique(c(names1,names2))
  mat = matrix(nrow=length(names),ncol=length(names))
  dimnames(mat) = list(names,names)
  for (i in 1:length(names1)) {
    mat[names1[i],names2[i]] = mat[names2[i],names1[i]] = vals[i]
  }
  return(mat)
}

corrplot.regions <- function(M, P,type="upper",order="hclust",method="circle",
                             cl.pos=NULL,tl.pos=NULL,pch.col=add.alpha(gray(.5),.7),...) {
  
  # marking independent.comparisons 
  BG = matrix(0,nrow=nrow(M),ncol=ncol(M))
  dimnames(BG) = dimnames(M)
  for (i in 1:nrow(M)) {
    for(j in 1:ncol(M)) {
      n.unique.pops = length(unique(unlist(c(strsplit(colnames(M)[i],"\n"),strsplit(colnames(M)[j],"\n")))))
      if (n.unique.pops==4) {
        BG[i,j] = 1
      }
    }
  }
  pdf(file = NULL); M.order=corrplot(M,order=order); dev.off()
  BG = BG[rownames(M.order),colnames(M.order)]
  diag(M) = 1
  M = M[rownames(M.order),colnames(M.order)]
  if (!is.null(P)) { diag(P) = 0;  P = P[rownames(M.order),colnames(M.order)]}
  corrplot(BG,method="square",tl.pos = tl.pos,tl.col="black",is.corr=F, cl.lim=c(0,1),
           diag=F,col=add.alpha("gold",.1),type=type,cl.pos=cl.pos)
  # if (method=="pie") { insig="blank" } else { insig="pch" }
  corrplot(M,method=method,order=order,tl.pos = "n", addrect = 4,is.corr=F,
           cl.lim=c(0,1),diag=F,bg="#FFFFFF00",
           add=T,outline=add.alpha("darkblue",.3),type=type,cl.pos=cl.pos,pch.col=pch.col,
           p.mat=P,sig.level = 0.05,...)
}




shorten.names.for.corrplot = function(names) {
  # browser()
  names=gsub("fst\\.","",names)
  names=gsub("dxy\\.","",names)
  names = unlist(lapply(strsplit(names,":"),
                        function(X){paste(substring(X,1,3),".",sep="",collapse="\n")}))
  names = gsub("Iri","Ire",names)
  names = gsub("Hyp","*Hyp",names)
  names = gsub("Alb","*Alb",names)
  return(names)
}

get.corrplot.matrices = function(over.table,type,remove.Irish=F,correct.P=T) {
  if (type=="22") {
    M = matrix(pmax(over.table$`Prop. Outlier Regions Shared 1`,
                    over.table$`Prop. Outlier Regions Shared 2`),nrow=1,
               dimnames = list(NULL,over.table[,1]))
    P = matrix(as.numeric(over.table[,grep("P-value",colnames(over.table))[1]]),nrow=1,
               dimnames = list(NULL,over.table[,1]))
    # diag(P)=diag(M)=1
  } else if (type=="11") {
    M = pairs.to.matrix(over.table[,1],over.table[,2],
                        pmax(over.table$`Prop. Outlier Regions Shared 1`,
                             over.table$`Prop. Outlier Regions Shared 2`))
    P = pairs.to.matrix(over.table[,1],over.table[,2],
                        as.numeric(over.table[,grep("P-value",colnames(over.table))[1]]))
    diag(P)=diag(M)=1
  } else if (type=="12") {
    name.split = strsplit(substring(over.table$comp,5),":")
    rn = unlist(lapply(name.split,
                       function(X){paste(substring(X,1,3),".",sep="",collapse="\n")}))
    P=M = matrix(nrow=2,ncol=length(unique(rn)),dimnames=list(c(1,2),unique(rn)))
    for (i in 1:ncol(M)) {
      comp = colnames(M)[i]
      pops = unique(unlist(name.split[which(comp==rn)]))
      M[1,i] = with(over.table[over.table$Pop==pops[1] & rn==comp,],
                    pmax(`Prop. Outlier Regions Shared 1`,`Prop. Outlier Regions Shared 2`))
      M[2,i] = with(over.table[over.table$Pop==pops[2] & rn==comp,],
                    pmax(`Prop. Outlier Regions Shared 1`,`Prop. Outlier Regions Shared 2`))
      P[1,i] = as.numeric(over.table[over.table$Pop==pops[1] & rn==comp,6])
      P[2,i] = as.numeric(over.table[over.table$Pop==pops[2] & rn==comp,6])
    }
    names=colnames(M)
    names=gsub("Iri","Ire",names)
    names = gsub("Hyp","*Hyp",names)
    names = gsub("Alb","*Alb",names)
    colnames(M)=colnames(P)=names
    # rownames(M)=rownames(P)=c("Low dXY\nLow π",
    #                           "Low dXY\nLow π")
    # P=t(P);M=t(M)
  }
  if (remove.Irish) {
    # browser()
    if (nrow(M)>1) {
      rows.M = -matches("^Ire|^Iri|Irish:",vars=rownames(M))
      rows.P = -matches("^Ire|^Iri|Irish:",vars=rownames(P))
    } else {rows.M=rows.P=1}
    M = t(as.matrix(M[rows.M,-matches("^Ire|^Iri|Irish:",vars=colnames(M))]))
    P = t(as.matrix(P[rows.P,-matches("^Ire|^Iri|Irish:",vars=colnames(P))]))
  }
  # browser()
  if (type!="12") {
    colnames(M)=shorten.names.for.corrplot(colnames(M))
    colnames(P)=shorten.names.for.corrplot(colnames(P))
  }
  if (nrow(M)>1) {
    if (type!="12") {
      rownames(M)=shorten.names.for.corrplot(rownames(M))
      rownames(P)=shorten.names.for.corrplot(rownames(P))
    }
    # pop.names.plus[match(pop.names.plus.order,pop.names.plus)]
    M = M[order(rownames(M)),order(colnames(M))]
    P = P[order(rownames(P)),order(colnames(P))]
  } else {
    M = M[,order(colnames(M))]
    P = P[,order(colnames(P))]
  }
  if (correct.P) {
    P.orig = P
    if (type=="22") {
      P = p.adjust(P,method="BH")
    } else if (type=="11") {
      P[upper.tri(P)] = p.adjust(P[upper.tri(P)],method="BH")
      P[lower.tri(P)] = p.adjust(P[lower.tri(P)],method="BH")
    } else if (type=="12") {
      P = p.adjust(P,method="BH")
      P = matrix(P,nrow=nrow(P.orig),ncol=ncol(P.orig))
    }
    dimnames(P) = dimnames(P.orig)
    P.orig = t(P.orig)
  } else { P.orig=NULL }
  return(list(t(M),t(P),P.orig))
}

corrplot.gen = function(corr,p.mat=NULL,
                        method="number",
                        cl.pos="n",
                        tl.col="black",
                        order="original",
                        pch.col=add.alpha(gray(.5),.7),...) {
  corrplot(corr=corr,p.mat=p.mat,method=method,order=order,cl.pos=cl.pos,tl.col=tl.col,pch.col=pch.col,...)
}

plot.outlier.similarities = function(comp,comp.type,corrplot.type,fig.dim,
                                     final=F,na.rm=T,remove.Irish=F,p.correct=T,
                                     correct.P=T,type="full",
                                     ...) {
  over.table = overlap.tables[[comp]][[comp.type]]
  if (na.rm) {
    over.table = over.table[apply(over.table,1,function(x) !all(is.na(x[2:3]))),]
  }
  if (p.correct) {
    over.table[,grep("P-value",colnames(over.table))[1]] = 
      p.adjust(as.numeric(over.table[,grep("P-value",colnames(over.table))[1]]),method="BH")
  }
  mats = get.corrplot.matrices(over.table,type=corrplot.type,remove.Irish=remove.Irish,correct.P=correct.P)
  M = mats[[1]]; P = mats[[2]]; P.orig = mats[[3]]
  if (final) { fig.loc="final" } else { fig.loc="scratch" }
  cairo_pdf(paste("figures/",fig.loc,"/corrplot_",comp,"_",comp.type,".pdf",sep=""),
            width=fig.dim[1],height=fig.dim[2])
  if (comp %in% c("fst","dxy")) {
    # cp = corrplot.regions(corr=M,p.mat=P,order="original",...)  
    cp = corrplot.regions(M=M,P=P,method="number",type=type,order="original",...)
  } else {
    cp = corrplot.gen(corr=M,p.mat=P,type=type,...)  
  }
  dev.off()
  # return(mean(M[upper.tri(M) | lower.tri(M)],na.rm=T))
  write.csv(M,paste("tables/M",comp,comp.type,"csv",sep="."))
  write.csv(P,paste("tables/P",comp,comp.type,"csv",sep="."))
  res = list(M,P,P.orig)
  names(res) = c("vals","p.vals.correct","p.vals.orig")
  return(res)
}



genomic.scatterplot = function(plot.order,sub.plot.list,all.data,title.text.list,
                               col.inds,nm,stat.bquote.x,main.plot.layout,
                               sub.plot.layout,second.stat=NULL,fst.cols=F,abline=F,stat.bquote.y=NULL,
                               col.inds.2=NULL,xlim.equals.ylim=F,lim.probs=c(0,1),condition.string=NULL,
                               span=.4,include.fst=F,stat.pos="lowerright") {
  all.data=as.data.frame(all.data)
  setTimeLimit()
  plotList=list()
  # browser()
  for (j in 1:length(plot.order)) {
    if (is.null(second.stat)) {
      i = plot.order[[j]]
      x = all.data[,i[1]]
      y = all.data[,i[2]]
      df = data.frame(x,y)
    } else {
      if (!second.stat%in%c("fst","dxy","cds")) {
        pair = plot.order[[j]]
        pops = strsplit(unlist(lapply(pair,function(x) { 
          unlist(lapply(strsplit(x,split="\\."),function(y) tail(y,1))) })),":")[[1]]
        x = all.data[,pair]
        y = (all.data[,paste(second.stat,pops[1],sep=".")]+all.data[,paste(second.stat,pops[2],sep=".")])/2
        fst = all.data[,paste("fst",substring(pair,4),sep="")]
        df = data.frame(x,y,fst)  
      } else if (second.stat=="cds") {
        i = plot.order[[j]]
        x = all.data[,"cds"]
        y = all.data[,i]
        df = data.frame(x,y)
      } else {
        i = plot.order[[j]]
        i2 = paste(second.stat,unlist(lapply(i,function(x) { 
          unlist(lapply(strsplit(x,split="\\."),function(y) tail(y,1))) })),sep=".")
        x = all.data[,i]
        y = all.data[,i2]
        df = data.frame(x,y)
      }
      
      if (include.fst) {
        df$fst = all.data[,paste("fst",substring(i,4),sep="")]
      }
      
      
      if (!is.null(condition.string)) {
        df=df[eval(parse(text=condition.string)),]
      }
      
      # browser()
    }
    title.text = title.text.list[[j]]
    spear.cor = suppressWarnings(cor.test(df$x,df$y,method="spearman"))
    if (spear.cor$estimate<0) {
      line.col="dodgerblue2"
    } else {
      line.col="goldenrod1"
    }
    if (spear.cor$p.value<0.05) { # Correct for multiple comparisons?
      spear.col="black"
    } else {
      spear.col="gray"
    }
    if (is.null(stat.bquote.y)) {
      stat.bquote.y=stat.bquote.x
    }
    # browser()
    xrange = quantile(unlist(all.data[,col.inds]),c(lim.probs[1],lim.probs[2])) # range(all.data[,col.inds]) #+.15    # .04
    a.bit = abs(diff(xrange))*.15
    xrange = c(xrange[1]-a.bit,xrange[2]+a.bit)
    if (is.null(second.stat)) {
      yrange = xrange
    } else {
      if (xlim.equals.ylim) {
        yrange = quantile(unlist(all.data[,c(col.inds,col.inds.2)]),c(lim.probs[1],lim.probs[2])) 
        a.bit = abs(diff(yrange))*.15
        yrange = c(yrange[1]-a.bit,yrange[2]+a.bit)  
        xrange = yrange
        
      } else {
        yrange = quantile(unlist(all.data[,col.inds.2]),c(lim.probs[1],lim.probs[2])) 
        a.bit = abs(diff(yrange))*.15
        yrange = c(yrange[1]-a.bit,yrange[2]+a.bit)  
      }
      
    }
    if (stat.pos=="upperright") {
      stat.y = yrange[2]-abs(diff(yrange))/20
    } else if (stat.pos=="lowerright") {
      stat.y = yrange[1]+abs(diff(yrange))/20
    } else { stop("Invalid stat placement") }
    stat.y.placement = 
    smooth.col = add.alpha("firebrick2",.8)
    gg = ggplot(df, aes(x = x, y = y)) + geom_point(size=.1,alpha=.6) + theme_classic() + 
      geom_smooth(method="lm",se=F,colour=line.col) + # labs(title="Title") +
      geom_smooth(method="loess",se=F,span=span,colour=smooth.col,size=.8) + 
      scale_x_continuous(limits=xrange,labels=abbreviate) +
      scale_y_continuous(limits=yrange) + 
      xlab(bquote(.(stat.bquote.x) ~ .(title.text[1]))) +
      ylab(bquote(.(stat.bquote.y) ~ .(title.text[2]))) +
      # annotate("text",x=0,y=ymax,label=toupper(letters)[j],hjust = 0) +
      annotate("text",x=xrange[2],y=stat.y,label=paste("rho ==",sprintf("%0.2f",spear.cor$estimate)),
               hjust = 1,size=3,parse=T,colour=spear.col)
    if (fst.cols) {
      gg = gg + scale_color_continuous(name=bquote(F["ST"]),                 
                                       # breaks = c(0,0.5,1),     
                                       # labels = c(0,0.5,1),
                                       low = "blue",              
                                       high = "red")   +          
        theme(legend.position=c(0.01,.9),
              legend.justification=c(0,1),
              legend.direction = "vertical",
              legend.key.height = unit(.05,"npc"))
    }
    if (abline) {
      gg=gg+geom_abline()
    }
    plotList=c(plotList,list(gg))
  }
  
  # tiff(paste("figures/scratch/",nm,"_compare.tiff",sep=""),
  #      width=2*main.plot.layout[2],height=2*main.plot.layout[1],
  #      units = "in",res=400)
  cairo_pdf(paste("figures/scratch/",nm,"_compare.pdf",sep=""),
       width=2*main.plot.layout[2],height=2*main.plot.layout[1])
  toPlot = plotList
  # for (i in 1:length(plotList)) {
  #   toPlot[[i]] = toPlot[[i]] + annotate("text",x=xrange[1],y=yrange[2],label=toupper(letters)[i],hjust = 0)
  # }
  grobList = lapply(toPlot,ggplotGrob)
  do.call("grid.arrange",c(grobList,nrow=main.plot.layout[1],ncol=main.plot.layout[2]))
  dev.off()
  
  for (k in 1:length(sub.plot.list)) {
    cairo_pdf(paste("figures/scratch/",nm,"_compare_sub_",k,".pdf",sep=""),
         width=sub.plot.layout[2]*2,height=sub.plot.layout[1]*2)
    # tiff(paste("figures/scratch/",nm,"_compare_sub_",k,".tiff",sep=""),
    #      width=sub.plot.layout[2]*2,height=sub.plot.layout[1]*2,units = "in",res=400)
    subPlotList = plotList[which(plot.order%in%sub.plot.list[[k]])]
    for (i in 1:length(subPlotList)) {
      subPlotList[[i]] = subPlotList[[i]]+ annotate("text",x=xrange[1],y=yrange[2],label=toupper(letters)[i],hjust = 0)
    }
    grobList = lapply(subPlotList,ggplotGrob)
    do.call("grid.arrange",c(grobList,nrow=sub.plot.layout[1],ncol=sub.plot.layout[2]))
    dev.off()
  }
}


plot.order.2.title.text = function(plot.order,fst.type=F) {
  if (fst.type) {
    title.text.list=lapply(lapply(plot.order,function(x) { unlist(lapply(strsplit(x,split="\\."),function(y) tail(y,1))) }),
                           function(y) { unlist(lapply(strsplit(y,":"),function(X){ paste(substring(X,1,3),".",sep="",collapse=" & ")})) })
    title.text.list = lapply(title.text.list, function(x) gsub("Iri","Ire",x))
  } else {
    title.text.list = lapply(plot.order,function(x) { unlist(lapply(strsplit(x,split="\\."),function(y) tail(y,1))) })
    title.text.list = lapply(title.text.list, function(x) gsub("Irish","Ireland",x))
    title.text.list = lapply(title.text.list, function(x) gsub("Alb$","F. alb.",x))
    title.text.list = lapply(title.text.list, function(x) gsub("Hyp$","F. hyp.",x)) 
  }
  return(title.text.list)
}

