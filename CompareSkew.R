#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------     
#
# Script used to conduct comparisons of distributional skew and kurtosis of FST distributions as used in
# Sendell-Price et al. (2018) The genomic landscape of divergence across the speciation
# continuum in island-colonising birds.
#
#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------   

setwd("YOUR FOLDER DIRECTORY WITH ALL FILES")


#STEP 1 - Comparing skew and kurtosis of randomly generated data

# number of loci
n.loci <- 10000
# now create a marker for locus name
locus.name <- rep(c(1:n.loci),times=2)
# I am just randomly generating some fst distributions here. So this might be FST between population 1 and 2
fst.dist.1 <- runif(n.loci)
# and this one between population 1 and 3
fst.dist.2 <- runif(n.loci)
# I now join these two distributions together
all.fsts <- c(fst.dist.1,fst.dist.2)
# next I create a factor that tells me whether the fst value is for the first or second comparision
dist.fact  <- as.factor(rep(c(1,2),each=n.loci))
# now create a dataframe to put all this in
data.frame <- data.frame(n.loci,locus.name,all.fsts,dist.fact)
# take a look at the top of it
head(data.frame)


# load library to calculate skew and kurtosis
# install.packages("e1071") -- i've commented this out as i have already downloaded the package
library(e1071)
# now calculate skewness for each observed FST distribution
obs.skew <- tapply(all.fsts,dist.fact,skewness)
obs.kurt <- tapply(all.fsts,dist.fact,kurtosis)



# Step 2 - doing the above for actual FST data:

Comp1 <- read.delim("FSTs/SI-CI.fst") # FST file for Comp1
Comp2 <- read.delim("FSTs/SI-FP.fst") # FST file for Comp2

Comp1$FST[which(Comp1$FST<0)]=0 #Weir's FST can give values <0 so we reset to 0
Comp1$Comp_Code <- "1"
Comp2$FST[which(Comp2$FST<0)]=0 #Weir's FST can give values <0 so we reset to 0
Comp2$Comp_Code <- "2"

All_FST <- Comp1
All_FST$Comp2_FST <- Comp2$FST
All_FST$Comp2_Code <- Comp2$Comp_Code
All_FST <- na.omit(All_FST)

FSTs <- c(All_FST$FST,All_FST$Comp2_FST)
dist.fact <- as.factor(c(All_FST$Comp_Code, All_FST$Comp2_Code))

library(e1071)
# now calculate skewness for each observed FST distribution
obs.skew <- tapply(FSTs,dist.fact,skewness)
obs.kurt <- tapply(FSTs,dist.fact,kurtosis)

# next take the absolute difference between the two skews and the two kurtosises
obs.abs.diff.skew <- abs(obs.skew[1]-obs.skew[2])
obs.abs.diff.kurt <- abs(obs.kurt[1]-obs.kurt[2])
# how many randomisations to conduct
n.samples <- 10000
# create arrays to store the randomised test statistics in
rand.abs.skew <- array(0,n.samples)
rand.abs.kurt <- array(0,n.samples)
# Now I randomise the factor levels and repeat each the above steps to generate a random distribution
# of the absolute difference statistic to compare with the observed one above
for (i in 1:n.samples){
  # shuffle dist.fact
  rand.dist.fact <- sample(dist.fact)
  # calculate skews and kurtosises for the randomised distribution
  rand.skew <- tapply(FSTs,rand.dist.fact,skewness)
  rand.kurt <- tapply(FSTs,rand.dist.fact,kurtosis)
  # next take the absolute difference between the two skews and the two kurtosises
  rand.abs.skew[i] <- abs(rand.skew[1]-rand.skew[2])
  rand.abs.kurt[i] <- abs(rand.kurt[1]-rand.kurt[2])
}
# we now have a load of randomised test statistics and one observed test statistic.  In this
# example we should expect the observed skew and kurtosis to sit somewhere near the left of 
# the distribution of randomised test statistics, simply because the two initial distributions 
# of fst distributions are expected to be similar, and consequently the absolutely value of 
# their tests statistics should be close to 0 -- that's the runif command higher up

quartz() # delete on a windoze machine
# do some plotting
par(mfrow=c(1,2))
hist(rand.abs.skew,xlab='Skew',main='Randomised test stat')
abline(v=obs.abs.diff.skew,col='red')
hist(rand.abs.kurt,xlab='Kurtosis',main='Randomised test stat')
abline(v=obs.abs.diff.kurt,col='red')

# now we want to find the percentile of the distribution the observed test statistic lies on.
# you are hoping for a value of >0.95
percentile.skew <- ecdf(rand.abs.skew)
percentile.kurt <- ecdf(rand.abs.kurt)
percentile.skew(obs.abs.diff.skew)
percentile.kurt(obs.abs.diff.kurt)

