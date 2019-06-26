##################################################################################################################
# NOTES:
# Sendell-Price et al (2019) The genomic landscape of divergence across the speciation continuum in an 
# island-colonising bird

# The following script counts the number of PCAdapt outliers within genomic islands and genomic valleys.
# Requires output from Island_Valley_Detection.R and Detect_PCAdapt_Outliers.R
##################################################################################################################

for (COMP in c("CI_SI","SI_FP","HI_ML","GT_LF","VN_LH")){
  
  #Load outlier data for comparison
  Outliers <- read.delim(paste0("PCAdapt/MAF0.2_Alpha0.001/",COMP,"_outliers.txt"), sep = " ") 
  #Apply some filtering -- include only outlier SNPs and remove outliers on unmapped chromosomes "_Un"
  Outliers <- Outliers[!grepl("Neutral", Outliers$STATUS),]
  Outliers <- Outliers[!grepl("_Un", Outliers$ZFCHROM),]
  
  #Report number of outliers for comparison
  print(paste0(COMP," ##################################################"))
  print(paste0("Outliers = ",nrow(Outliers)))
  
  #Load island and valley positions
  Islands.Valleys <- read.delim(paste0("Island_valley_detection/50kb/June2019/",COMP,"_IslandValleyStats.txt")) 
  
  #Split into seperate dataframes for Islands and valleys
  Valleys <- Islands.Valleys[grepl("low", Islands.Valleys$type),]
  Islands <- Islands.Valleys[grepl("high", Islands.Valleys$type),]
  
  #Create empty dataframes which we will add outliers to if within islands or valleys
  In.Islands <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(In.Islands) <- c("ZFCHROM","ZFPOS","ZLCHROM","ZLPOS","PCAdapt.PVals","STATUS")
  In.Valleys <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(In.Valleys) <- c("ZFCHROM","ZFPOS","ZLCHROM","ZLPOS","PCAdapt.PVals","STATUS")
  
  #Find outliers (if any) contained within genomic islands
  for (row in 1:nrow(Islands)){
    temp <- filter(Outliers, ZFCHROM == as.character(Islands$chr.start[row]) & ZFPOS >= Islands$pos.start[row] & ZFPOS <= Islands$pos.end[row])
    In.Islands <- rbind(In.Islands, temp)
  }
  
  #Print number of outliers in islands
  print(paste0("of which ",nrow(In.Islands), " are in islands"))
  
  #Find outliers (if any) contained within genomic islands
  for (row in 1:nrow(Valleys)){
    temp <- filter(Outliers, ZFCHROM == as.character(Valleys$chr.start[row]) & ZFPOS >= Valleys$pos.start[row] & ZFPOS <= Valleys$pos.end[row])
    In.Valleys <- rbind(In.Valleys, temp)
  }
  
  #Print number of outliers in valleys
  print(paste0("and ", nrow(In.Valleys), " are in valleys"))
  
}
