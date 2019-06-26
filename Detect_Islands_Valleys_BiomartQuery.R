##################################################################################################################
# NOTES:
# Sendell-Price et al (2019) The genomic landscape of divergence across the speciation continuum in an 
# island-colonising bird

# The following script detects genomic islands and genomic valleys across the genomic landscape of divergence 
# populations via permutation. This script is based on that used by Van Doren et al. (2017) Correlated patterns of 
# genetic diversity and differentiation across an avian family
##################################################################################################################

#Load required packages and functions
library(ape)
library(dplyr)
library(reshape2)
source("Island_Detecting_Functions.R")

##################################################################################################################
# Define further functions
##################################################################################################################

peak.permutation.test = function(dataset,
                                 data.col,
                                 bandwidth,
                                 N,
                                 win.size,
                                 peak.pop,
                                 plot = T,
                                 treat.Z = c("same", "exclude", "only", "separate"),
                                 retain.peaks.greater.than = 0,
                                 fill.gaps.less.than = 10,
                                 kernel.type="box",
                                 chr.col.name = "ficed_chr",
                                 Z.chr.names = c("chr_Z", "chr_Z_un")) {
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

################################################################################################

process.stats = function(data,
                         hide.small.chrs = T,
                         na.rm = T,
                         remove.unmapped.scaffolds = T,
                         fst=T) {
  data$mean.pos = apply(data[,c("BIN_START","BIN_END")],1,function(X) { mean(c(X["BIN_START"]-1,X["BIN_END"])) })
  if (fst) { data = data[,c("CHR","mean.pos","N_VARIANTS",paste0("Fst_",Pop1,"_",Pop2), paste0("dxy_",Pop1,"_",Pop2), paste0("pi_",Pop1),paste0("pi_",Pop2))] }
  else { data = data[,c("CHR","mean.pos","N_VARIANTS",paste0("Fst_",Pop1,"_",Pop2), paste0("dxy_",Pop1,"_",Pop2), paste0("pi_",Pop1),paste0("pi_",Pop2))] }
  colnames(data) = c("chr","mean.pos","n.snps","mean.fst","mean.dxy","pi_Pop1","pi_Pop2")
  if (na.rm) {
    data = data[complete.cases(data), ]
  }
  data$index = 1:nrow(data)
  data$mean.fst[which(data$mean.fst<0)]=0 #Weir's FST can give values <0 so we reset to 0
  # return(process.with.chrs(data, na.rm, hide.small.chrs,remove.unmapped.scaffolds))
  return(data)
}

################################################################################################

ash.ppt = function(the.data,pop,win.size,bandwidth=20) {
  the.data = as.data.frame(the.data)
  set.seed(123)
  res = peak.permutation.test(the.data,
                              data.col="mean.fst",
                              bandwidth=bandwidth,
                              N=25000,   #### <<<<-------- number of permutation iterations 
                              win.size=win.size,
                              peak.pop=pop,
                              plot=T,
                              treat.Z="separate",
                              chr.col.name = "chr",
                              Z.chr.names = c("Z", "Z_un"))
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
# mean.val = mean value of Fst/dxy over the peak/valley
# pop = comparison
# chr.start = chromosome of the first window in the peak/valley
# pos.start = position at which the first window in the peak/valley starts
# chr.end = chromosome of the last window in the peak/valley (should be the same as chr.start, included just in case something goes wrong)
# pos.end = position at which the last window in the peak/valley ends


#################################################################################################
#This for loop will conduct island / valley detection for each comparison, extract stats from
#islands and valleys and write output
#################################################################################################

#First load windowed stats (FST, DXY, and Pi)
Windowed.Stats <- read.csv("FST_Pi_Dxy_All_Comps_50kb.csv")

#First specify window size (in KB) that we will use:
window_size <- 50

#Find centre of window:
Windowed.Stats$mean.pos = Windowed.Stats$BIN_START+(window_size*(1000/2)-1)
#When estimate FST negative values can sometimes be generated. We will reset these to 0. 
Windowed.Stats[Windowed.Stats <=0] = 0

#The loop
for (COMP in  c("CI_SI","SI_FP","HI_ML","GT_LF","VN_LH")){
  
  #Split comparison into individual populations
  Populations <- colsplit(COMP, "_", names = c("Pop1", "Pop2"))
  
  #Define each population in the comparison of interest
  Pop1 <- Populations$Pop1
  Pop2 <- Populations$Pop2
  
  #Create new dataframe containing info for comparison of interest
  data = Windowed.Stats[,c("CHR","mean.pos","N_VARIANTS",paste0("Fst_",Pop1,"_",Pop2), paste0("dxy_",Pop1,"_",Pop2), paste0("pi_",Pop1),paste0("pi_",Pop2),"BIN_START","BIN_END")]
  
  colnames(data) = c("chr","mean.pos","n.snps","mean.fst","mean.dxy","pi_Pop1","pi_Pop2","BIN_START","BIN_END")
  
  data <- na.omit(data)
  data$index = 1:nrow(data)
  
  #Run island / valley detection on dataframe
  Comp = ash.ppt(data, pop=paste0(Pop1,"_",Pop2), win.size=window_size*1000)
    Comp
  
  #Create new empty columns
  Comp$mean.dxy <- NA
  Comp$mean.Pi_Pop1 <- NA
  Comp$mean.Pi_Pop2 <- NA
  
  #Fill columns with values calculated for each island/valley
  for (row in 1:length(Comp$chr.start)) {
    data_within_islandvalley <- data[(data$chr == Comp$chr.start[row] & data$BIN_START >= Comp$pos.start[row] & data$BIN_END <= Comp$pos.end[row]), ] 
    Comp$mean.dxy[row] <- mean(data_within_islandvalley$mean.dxy)
    Comp$mean.Pi_Pop1[row] <- mean(data_within_islandvalley$pi_Pop1)
    Comp$mean.Pi_Pop2[row] <- mean(data_within_islandvalley$pi_Pop2)
  }
  
  #Save output to file
  write.table(Comp, paste0(Pop1,"_",Pop2,"_IslandValleyStats.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

  #Create indpendent datasets for islands and valleys
  islands <- Comp[(Comp$type == "high"), ]
  valleys <- Comp[(Comp$type == "low"), ]
  
  
  #Now for each 50kb window define if contained within island
  for (row in 1:nrow(islands)) {
    
    within_island <- data[(data$chr == islands$chr.start[row] & data$BIN_START >= islands$pos.start[row] & data$BIN_END <= islands$pos.end[row]), ] 
    within_island$OutlierStatus <- "Island"
    data <- full_join(data, within_island)
  }
  
  #Now for each 50kb window define if contained within valley
  for (row in 1:nrow(valleys)) {
    
    within_valley <- data[(data$chr == valleys$chr.start[row] & data$BIN_START >= valleys$pos.start[row] & data$BIN_END <= valleys$pos.end[row]), ] 
    within_valley$OutlierStatus <- "Valley"
    data <- full_join(data, within_valley)
  }

  #Fill NAs with "background" and write to file.
  data$OutlierStatus[is.na(data$OutlierStatus)] <- "Background"
  write.table(data, paste0(Pop1,"_",Pop2,"_windowstats_plus_IslandValley_status.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
}


#################################################################################################
#This for loop will query the biomart database to identify all genes contained within genomic
#islands and genomic valleys.
#################################################################################################

#Load biomaRt library
library(biomaRt)

#Get Zebra Finch genes from biomart:
ensembl = useEnsembl(biomart="ensembl", dataset="tguttata_gene_ensembl")

#Pull out info of interest. 
all.genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), 
                   mart = ensembl)

#For loop which will conduct query searches for each comparison. 
for (COMP in  c("CI_SI","SI_FP","HI_ML","GT_LF","VN_LH")){
  
  print(paste("Identifying genes for ",COMP))
  
  #re-load island/valley positions for given comparison
  islands_and_valleys <- read.delim(paste0("Island_valley_detection/50kb/June2019/",COMP,"_IslandValleyStats.txt"))
  
  #For speed we will copy all.genes to a new dataframe - this will be overrighten later so need to maintain all.genes
  all.genes2 <- all.genes
  
  #Create indpendent datasets for islands and valleys
  islands <- islands_and_valleys[(islands_and_valleys$type == "high"), ]
  valleys <- islands_and_valleys[(islands_and_valleys$type == "low"), ]
  
  #Search genes within each island
  for (row in 1:nrow(islands)) {
    
    genes.within.island <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), 
                                 mart = ensembl,
                                 filter = c('chromosome_name','start','end'),
                                 values = list(islands$chr.start[row], islands$pos.start[row], islands$pos.end[row]))
    
    genes.within.island$chromosome_name <- as.character(genes.within.island$chromosome_name)
    genes.within.island$hgnc_symbol <- as.character(genes.within.island$hgnc_symbol)
    
    
    #This if else statement prevents errors when there is no genes within island/valley. 
    if(nrow(genes.within.island) == 0){
      print("data.frame is empty")
    }else{
      genes.within.island$Within <- "Island"
      
      #join within column to all.genes2
      all.genes2 <- full_join(all.genes2, genes.within.island)
    }
  }
  
  
  #Do the same for valleys.
  for (row in 1:nrow(valleys)) {
    
    genes.within.valley <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), 
                                 mart = ensembl,
                                 filter = c('chromosome_name','start','end'),
                                 values = list(valleys$chr.start[row], valleys$pos.start[row], valleys$pos.end[row]))
    
    genes.within.valley$chromosome_name <- as.character(genes.within.valley$chromosome_name)
    genes.within.valley$hgnc_symbol <- as.character(genes.within.valley$hgnc_symbol)
    
    if(nrow(genes.within.valley) == 0){
      print("data.frame is empty")
    }else{
      genes.within.valley$Within <- "Valley"
      all.genes2 <- full_join(all.genes2, genes.within.valley)
    }
  }
  
  #Write to file (omitting genes outside islands/valleys using na.omit())
  write.table(na.omit(all.genes2), paste0("Island_valley_detection/50kb/June2019/",COMP,"_genes_within_islandsvalleys.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
}



