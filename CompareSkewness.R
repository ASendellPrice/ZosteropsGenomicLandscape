##################################################################################################################
# NOTES:
# Sendell-Price et al (2019) The genomic landscape of divergence across the speciation continuum in an 
# island-colonising bird

# The following script tests for significant differences FST distributional skew between population
# comparisons via randomisation test.
##################################################################################################################

#Load required packages
library(e1071)
library(moments)
library(dplyr)

#Turns off scientific notation
options(scipen = 9999) 

#Load windowed stats:
Windowed.Stats.500kb <- na.omit(read.csv("Calculate_FST_dxy_Pi/500kb/FST_Pi_Dxy_All_Comps_500kb.csv"))

#Correct values <0 to 0:
Windowed.Stats.500kb[Windowed.Stats.500kb <=0] = 0

#Split into Autosome and Non-autosome datatsets:
Autosomes <- filter(Windowed.Stats.500kb, CHR != "Z")
ChrZ <- filter(Windowed.Stats.500kb, CHR == "Z")

#Define comparisons to compare:
Comp_1 = "GT_LF"
Comp_2 = "VN_LH"

#Define statistic (dxy or Fst)
Stat = "Fst"
#----------------------------------------------------------------------------------------
# FST Autosomes only
#----------------------------------------------------------------------------------------

#Create dataframe for Comp1 pulling out relevant column from Autosome d/f.
Comp1 <- Autosomes[paste0(Stat, "_",Comp_1)]
colnames(Comp1) <- "Divergence_Statistic"
Comp1$Pops <- Comp_1

#Do the same for Comp2
Comp2 <- Autosomes[paste0(Stat, "_",Comp_2)]
colnames(Comp2) <- "Divergence_Statistic"
Comp2$Pops <- Comp_2

#Join these together
Combined <- rbind(Comp1, Comp2) 

#Create two variable from combined dataframe
all.fsts <- Combined$Divergence_Statistic
dist.fact <- Combined$Pops

# now calculate skewness for each observed FST distribution
obs.skew <- tapply(all.fsts,dist.fact,skewness)

# next take the absolute difference between the two skews and the two kurtosises
obs.abs.diff.skew <- abs(obs.skew[1]-obs.skew[2])
# how many randomisations to conduct
n.samples <- 10000
# create arrays to store the randomised test statistics in
rand.abs.skew <- array(0,n.samples)
# Now I randomise the factor levels and repeat each the above steps to generate a random distribution
# of the absolute difference statistic to compare with the observed one above
for (i in 1:n.samples){
  # shuffle dist.fact
  rand.dist.fact <- sample(dist.fact)
  # calculate skews and kurtosises for the randomised distribution
  rand.skew <- tapply(all.fsts,rand.dist.fact,skewness)
  # next take the absolute difference between the two skews and the two kurtosises
  rand.abs.skew[i] <- abs(rand.skew[1]-rand.skew[2])
}

# do some plotting
par(mfrow=c(1,1))
hist(rand.abs.skew,xlab='Skew',main='Randomised test stat')
abline(v=obs.abs.diff.skew,col='red')


# now we want to find the percentile of the distribution the observed test statistic lies on.
percentile.skew <- ecdf(rand.abs.skew)
1-percentile.skew(obs.abs.diff.skew) # <--- returns Pvalue


#----------------------------------------------------------------------------------------
# FST Z Chromosome only
#----------------------------------------------------------------------------------------
#Create dataframe for Comp1 pulling out relevant column from ChrZ d/f.
Comp1 <- ChrZ[paste0(Stat, "_",Comp_1)]
colnames(Comp1) <- "Divergence_Statistic"
Comp1$Pops <- Comp_1

#Do the same for Comp2
Comp2 <- ChrZ[paste0(Stat, "_",Comp_2)]
colnames(Comp2) <- "Divergence_Statistic"
Comp2$Pops <- Comp_2

#Join these together
Combined <- rbind(Comp1, Comp2) 

all.fsts <- Combined$Divergence_Statistic
dist.fact <- Combined$Pops

# now calculate skewness for each observed FST distribution
obs.skew <- tapply(all.fsts,dist.fact,skewness)

# next take the absolute difference between the two skews and the two kurtosises
obs.abs.diff.skew <- abs(obs.skew[1]-obs.skew[2])
# how many randomisations to conduct
n.samples <- 10000
# create arrays to store the randomised test statistics in
rand.abs.skew <- array(0,n.samples)
# Now I randomise the factor levels and repeat each the above steps to generate a random distribution
# of the absolute difference statistic to compare with the observed one above
for (i in 1:n.samples){
  # shuffle dist.fact
  rand.dist.fact <- sample(dist.fact)
  # calculate skews and kurtosises for the randomised distribution
  rand.skew <- tapply(all.fsts,rand.dist.fact,skewness)
  # next take the absolute difference between the two skews and the two kurtosises
  rand.abs.skew[i] <- abs(rand.skew[1]-rand.skew[2])
}

# do some plotting
par(mfrow=c(1,1))
hist(rand.abs.skew,xlab='Skew',main='Randomised test stat')
abline(v=obs.abs.diff.skew,col='red')

# now we want to find the percentile of the distribution the observed test statistic lies on.
percentile.skew <- ecdf(rand.abs.skew)
1-percentile.skew(obs.abs.diff.skew) # <--- returns Pvalue

