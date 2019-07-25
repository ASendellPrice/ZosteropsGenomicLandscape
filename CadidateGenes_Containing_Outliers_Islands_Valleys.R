
setwd("~/Dropbox/DPhil/ZostGenomicLandscape/Post_Review/June2019/")

#Load candidate gene list - which contains chromosomal positions
Candidates <- read.delim("Candidate_Genes/ExistingCadidates.txt")


###########################################################################################################
# PART 1: DO CANDIDATE GENES CONTAIN OUTLIER SNPS?
###########################################################################################################


for (COMP in c("CI_SI","SI_FP","HI_ML","GT_LF","VN_LH")){
  
  print(paste0("For: ", COMP))
  
  #Load PCAdapt output for comparison
  Outliers <- read.delim(paste0("../PCAdapt/MAF0.2_Alpha0.001/",COMP,"_outliers.txt"), sep = " ") 
  #Apply some filtering -- include only outlier SNPs and remove outliers on unmapped chromosomes "_Un"
  Outliers <- Outliers[!grepl("Neutral", Outliers$STATUS),]
  Outliers <- Outliers[!grepl("_Un", Outliers$ZFCHROM),]
  
  for (row in 1:nrow(Candidates)){
    temp <- filter(Outliers, ZFCHROM == as.character(Candidates$CHROM[row]) & ZFPOS >= as.numeric(Candidates$START[row]) & ZFPOS <= as.numeric(Candidates$END[row]))
    
    if (nrow(temp) != 0) {
      print(paste0(Candidates$GENE[row]," contains ", nrow(temp),' SNPs'))
    } 
  }
  print("-------------------------------------------------")
  
}


###########################################################################################################
# PART 2: ARE CANDIDATE GENES WITHIN GENOMIC ISLANDS
###########################################################################################################

for (COMP in c("CI_SI","SI_FP","HI_ML","GT_LF","VN_LH")){
  
  print(paste0("For: ", COMP))
  
  #Load island and valley positions
  Islands.Valleys <- read.delim(paste0("Island_valley_detection/",COMP,"_IslandValleyStats.txt")) 
  #Split into seperate dataframes for Islands and valleys
  Valleys <- Islands.Valleys[grepl("low", Islands.Valleys$type),]
  Islands <- Islands.Valleys[grepl("high", Islands.Valleys$type),]
  
  #Create empty dataframes which we will add genes to if within islands or valleys
  In.Islands <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(In.Islands) <- c("ID","GENE","CHROM","START","END")
  In.Valleys <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(In.Valleys) <- c("ID","GENE","CHROM","START","END")
  
  
  for (row in 1:nrow(Islands)){
    Island.Range <- as.vector(c(Islands$pos.start[row], Islands$pos.end[row]))
    
    Candidates_on_chrom <- filter(Candidates, CHROM == as.character(Islands$chr.start[row]))
    
    if(nrow(Candidates_on_chrom) > 0){
      
      for (row.candidates in 1:nrow(Candidates_on_chrom)){
        Candidate_Range <- c(Candidates_on_chrom$START[row.candidates], Candidates_on_chrom$END[row.candidates])
        
        library(DescTools)
        
        Overlaps <- Island.Range %overlaps% Candidate_Range
        
        if (Overlaps == TRUE) {
          print(paste0("Island_",Islands$chr.start[row], ": ", Islands$pos.start[row], "-", Islands$pos.end[row]," contains: ", Candidates_on_chrom$GENE[row.candidates]))
        }
      }
      
    }
  }
  
  
  for (row in 1:nrow(Valleys)){
    Valley.Range <- as.vector(c(Valleys$pos.start[row], Valleys$pos.end[row]))
    
    Candidates_on_chrom <- filter(Candidates, CHROM == as.character(Valleys$chr.start[row]))
    
    if(nrow(Candidates_on_chrom) > 0){
      
      for (row.candidates in 1:nrow(Candidates_on_chrom)){
        Candidate_Range <- c(Candidates_on_chrom$START[row.candidates], Candidates_on_chrom$END[row.candidates])
        
        library(DescTools)
        Overlaps <- Valley.Range %overlaps% Candidate_Range
        
        if (Overlaps == TRUE) {
          print(paste0("Valley_",Valleys$chr.start[row], ": ", Valleys$pos.start[row], "-", Valleys$pos.end[row]," contains: ", Candidates_on_chrom$GENE[row.candidates]))
        }
      }
      
    }
  }
  print("-------------------------------------------------")
}
  

