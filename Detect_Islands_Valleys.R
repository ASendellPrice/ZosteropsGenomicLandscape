#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------     
#
# Script used to detect genomic islands of divergence and genomic valleys of similarity as used in
# Sendell-Price et al. (2018) The genomic landscape of divergence across the speciation
# continuum in island-colonising birds.
#
# Based on script developed by Benjamin Van Doren (benjamin.vandoren@zoo.ox.ac.uk).
#
#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------   


library(ape)
library(dplyr)

setwd("YOUR FOLDER DIRECTORY WITH ALL FILES")

source("Island_Detecting_Functions.R")

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
                                 kernel.type="box",
                                 chr.col.name = "SCAFFOLD",
                                 Z.chr.names = c("Chr_Z", "Chr_Z_Un")) {
  # data.col is number or name of column
  
  # browser()
  
  y = dataset[, data.col]
  x = dataset$index
  cs = dataset[,chr.col.name]
  stopifnot(!NA %in% cs)
  if (length(treat.Z) > 1) {
    treat.Z == "same"
  }
  stopifnot(treat.Z %in% c("same", "exclude", "only", "separate"))
  
  onZ = cs %in% Z.chr.names
  
  # browser()
  
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
      # browser()
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
    plot(x, y, col = gray(.4), pch = ".",xlab="Position",ylab="Genomic statistic")
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
    title(peak.pop)
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
    abline(v = which(over), col = add.alpha("orange", .03))
    abline(v = which(under), col = add.alpha("blue", .01))
    
    abline(v = which(over.per.peak),
           col = add.alpha("orange", .1))
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


#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------

process.vcftools.fst = function(data,
                                hide.small.chrs = T,
                                na.rm = T,
                                remove.unmapped.scaffolds = T,
                                weighted.fst=T) {
  data$mean.pos = apply(data[,c("BIN_START","BIN_END")],1,function(X) { mean(c(X["BIN_START"]-1,X["BIN_END"])) })
  if (weighted.fst) { data = data[,c(1,7,4,5)] }
  else { data = data[,c(1,7,4,6)] }
  colnames(data) = c("chr","mean.pos","n.snps","mean.fst")
  if (na.rm) {
    data = data[complete.cases(data), ]
  }
  data$index = 1:nrow(data)
  # return(process.with.chrs(data, na.rm, hide.small.chrs,remove.unmapped.scaffolds))
  return(data)
}

#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------

ash.ppt = function(the.data,pop,win.size,bandwidth=20) {
  the.data = as.data.frame(the.data)
  set.seed(123)
  res = peak.permutation.test(the.data,
                              data.col="mean.fst",
                              bandwidth=bandwidth,
                              N=1e4,   #### <<<<-------- number of permutation iterations 
                              win.size=win.size,
                              peak.pop=pop,
                              plot=T,
                              treat.Z="separate",
                              chr.col.name = "chr",
                              Z.chr.names = c("Chr_Z", "Chr_Z_Un"))
  # browser()
  res$chr.start = the.data$chr[res$beg]
  res$pos.start = the.data$mean.pos[res$beg]-win.size/2
  res$chr.end = the.data$chr[res$end]
  res$pos.end = the.data$mean.pos[res$end]+win.size/2
  return(res %>% select(-c(peak,beg,end)))
  
}
  
#####Explanation of table columns:
# type = whether it's a `peak' (high Fst = `high') or `valley' (low Fst = `low')
# len = length (in number of windows)
# mean.val = mean value of Fst over the peak/valley
# pop = comparison
# chr.start = chromosome of the first window in the peak/valley
# pos.start = position at which the first window in the peak/valley starts
# chr.end = chromosome of the last window in the peak/valley (should be the same as chr.start, included just in case something goes wrong)
# pos.end = position at which the last window in the peak/valley ends


#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------
# Example island valley detection 
#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------

data = read.delim("FST_500K/SI_FP_FST_500k.windowed.weir.fst",header=T)
data$MEAN_FST[which(data$MEAN_FST<0)]=0 #Weir's FST can give values <0 so we reset to 0

data = process.vcftools.fst(data,weighted.fst = F,
                                  remove.unmapped.scaffolds = F,
                                  hide.small.chrs = F)

500kb = ash.ppt(data, pop="Pop 1 vs. Pop 2", win.size=500000)
500kb
write.csv(500kb,"output file name.csv")

# summary
500kb$peak = T
peak.data.summary(500kb, type.return = "both", win.size.multiplier = 500)

