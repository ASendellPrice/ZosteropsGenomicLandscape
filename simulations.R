#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------     
## (c) Claudio Quilodran and Eric C. Anderson 2018
# Department of Zoology, University of Oxford: claudio.quilodran@zoo.ox.ac.uk
# Fisheries Ecology Division, Southwest Fisheries Science Center, National Marine Fisheries Service: eric.anderson@noaa.gov      
#
# Individual-based simulation of two population divergence as used in 
# Sendell-Price et al. (2018) The genomic landscape of divergence across the speciation
# continuum in island-colonising birds.
#
# For further info contact: claudio.quilodran@zoo.ox.ac.uk
#
#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------   
setwd("YOUR FOLDER DIRECTORY WITH ALL FILES")
library(progress)
source("MainFunction.R")

#########################
### Orginizing datset ###
#########################
Geno<-read.table("165-129K_ZF_ified_17.geno", header=F, comment.char = "")

Geno.t<-t(Geno[,-1])
Geno.t[1,1] <- "ID"
colnames(Geno.t) <- Geno.t[1,]
Geno.t<-Geno.t[-1,]
rownames(Geno.t) <- Geno.t[,1]
Geno.t<-Geno.t[,-1]


n1<-apply(Geno.t,1:2,function(i) strsplit(i,split='')[[1]][1])
n2<-apply(Geno.t,1:2,function(i) strsplit(i,split='')[[1]][2]) 

n1[n1 %in% "A"] <- 1
n1[n1 %in% "C"] <- 2
n1[n1 %in% "G"] <- 3
n1[n1 %in% "T"] <- 4
n1[n1 %in% "N"] <- NA

n1<-apply(n1,2, as.numeric)
rownames(n1)<-rownames(n2)


n2[n2 %in% "A"] <- 1
n2[n2 %in% "C"] <- 2
n2[n2 %in% "G"] <- 3
n2[n2 %in% "T"] <- 4
n2[n2 %in% "N"] <- NA

n2<-apply(n2,2, as.numeric)
rownames(n2)<-rownames(n1)

struct.t0<-list(n1,n2)

L<-list(struct.t0[[1]], struct.t0[[2]])
struct.t00<-array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)), dimnames=list(rownames(L[[1]]),colnames(L[[1]]), 1:2 ))

#########################
###  Parameter values ###
#########################

loci.pos<-as.numeric(colnames(struct.t00) )
loci.pos<-loci.pos-loci.pos[1] + 1
chromo_mb<-max(loci.pos)

crossover <- 1.5/100000000.0 #recombination rate (1.5cM/MB)
sex.ratio = 0.5 			#ratio of males and females
Fledged = 1.9				#number of fledglings by clutch
Survival.summer = 0.96 		#Survival rate at summer
Survival.autoumn = 0.76 		#Survival rate at autoumn
Survival.winter = 0.63 		#Survival rate at winter

Generation.time = 3 			#generation time (years)		
mean.fitness = Fledged*Survival.summer*Survival.autoumn*Survival.winter*Generation.time 

n.gens = 100	#number of simulated generations
ve= 0.0013		#density-dependent demographic effect
ini.ind = 400 	#initial number of individuals at each population
disp=0.001 		#migration rate

############################################
##### Generating initial individuals #######
############################################

dat0<-as.data.frame(rbind(struct.t0[[1]],struct.t0[[2]]) )
rownames(dat0)<-1:nrow(dat0)
snp<-lapply(1:ncol(dat0), function(i){ table(dat0[[i]], exclude=NA) })

set.seed(1)
Pop<-sapply(1:ini.ind, function(i){ 
	print(i)
		sapply(1:ncol(struct.t0[[1]]), function(j){ 
			n1<-sample(names(snp[[j]]),1,T, prob=snp[[j]])			
			n2<-sample(names(snp[[j]]),1,T, prob=snp[[j]]) 			
	return(list(n1, n2))
			})	
	 })


Pop<-t(Pop)
even<-seq_len(ncol(Pop)) %% 2 ==0

n1<-Pop[, even]
n2<-Pop[, !even]

n1<-apply(n1,2, as.numeric)
n2<-apply(n2,2, as.numeric)

L<-list(n1,n2)

start<-array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)), dimnames=list(1:ini.ind,colnames(dat0), 1:2 ))

start.1=start #Initial population 1
start.2=start #Initial population 2

#########################
###    Simulations    ###
#########################

out<-main.function(start.1, start.2, sex.ratio, mean.fitness, n.gens, loci.pos, disp, chromo_mb, ve, crossover) 


###########################
###    Analyzing Fst    ###
###    Python Script    ###
###########################   
###Orginizing dataset
#First population
P1<-out[[1]]
P1<-sapply(1:dim(P1)[2], function(i){
	paste(P1[,i,1], P1[,i,2], sep="")
})
colnames(P1)<-colnames(start)
P1 <-gsub("NA", "N", P1)
P1 <-chartr("1234", "ACGT ", P1)
P1 <-t(P1)
P1<-as.data.frame( cbind(rep("Chr_17", nrow(P1)), colnames(start), P1) )
colnames(P1)<-c("#CHROM", "POS", paste("P1_", 1:(ncol(P1)-2), sep="") ) 
rownames(P1)<-1:dim(P1)[1]

#Second population
P2<-out[[2]]
P2<-sapply(1:dim(P2)[2], function(i){
	paste(P2[,i,1], P2[,i,2], sep="")
})
colnames(P2)<-colnames(start)
P2 <-gsub("NA", "N", P2)
P2 <-chartr("1234", "ACGT ", P2)
P2 <-t(P2)
colnames(P2)<-paste("P2_", 1:ncol(P2), sep="")

#both populations in a single file
Data_P1_P2<-as.data.frame(cbind(P1,  P2 ) )
filename=paste("Output.geno", sep="")
write.table(Data_P1_P2, filename, sep="\t", row.name=F, quote=F) 

#Calling the Python script
pop1<-dim(out[[1]])[1]
pop2<-dim(out[[2]])[1]
filename2<-paste(filename, ".gz", sep="")

	system(paste("gzip", filename))
	system( paste("python -W ignore calculate_FST.py", filename2, pop1+1, pop2 +1) )
	
	filename3=paste("FST_",  filename2, ".csv", sep="")
	Fst<-read.csv(filename3 , na.strings = "nan")[,6]
	
	system( paste("rm -r", filename2) )
	system( paste("rm -r", filename3) )

#########################
#####    Figure    ######
#########################

h <- hist(Fst, breaks=seq(0, 1, length=15), freq=T) 
h$counts=h$counts/sum(h$counts)
plot(h, col="grey", xlim=c(0,1),  ylim=c(0, 1), xlab="Fst", ylab="Proportion", main=NULL, bty="l", xaxs="i",yaxs="i", las=1, border=F)




