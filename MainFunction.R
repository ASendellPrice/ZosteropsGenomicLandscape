#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------     
## (c) Claudio Quilodran and Eric C. Anderson 2018
# Department of Zoology, University of Oxford: claudio.quilodran@zoo.ox.ac.uk
# Fisheries Ecology Division, Southwest Fisheries Science Center, National Marine Fisheries Service: eric.anderson@noaa.gov      
#
# Functions required to run individual-based simulation of two population divergence as used in 
# Sendell-Price et al. (2018) The genomic landscape of divergence across the speciation
# continuum in island-colonising birds.
#
# For further info contact: claudio.quilodran@zoo.ox.ac.uk
#
#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------   

library(Rcpp)
sourceCpp("Funcs_rcpp.cpp")

##################
### Functions  ###
##################

# a function to generate the initial population structure
initial.struct          <- function(n.N1,n.l,n.a.l){ # n.N1 - initial N, n.l - N loci, n.a.l - alleles / locus
  rand.ints             <- sample(1:n.a.l,n.N1*n.l*2,TRUE)
  struct                <- array(rand.ints,c(n.N1,n.l,2)) # an array
  return(struct)
}

# mating function
resample <- function(x, ...) x[sample.int(length(x), ...)]
mating <- function(fitn, struct){
  females      <- which(fitn[,"sex"] %in% 1)
  males        <- which(fitn[,"sex"] %in% 2)
  try(if(length(females)<1 || length(males)<1 ) stop("Population extinction. There is not enough females or males for mating.", call.=F))  
  female.fitn  <- fitn[females,]
  female.stru  <- struct[females,,, drop = F]
  male.fitn    <- fitn[males,]
  mum.stru     <- female.stru[c(rep(1:length(females),times=female.fitn[,"fit"])),,, drop = F]
  dads         <- resample(males,sum(female.fitn[,"fit"]),TRUE,prob=male.fitn[,"fit"])
  dad.stru     <- struct[c(dads),,, drop = F]
  return(list(mum.stru,dad.stru))
}

##########################################
#' disperal function that avoids slow abind
#' @param pop1 pop struct 1
#' @param pop2 pop struct 2
#' @param rate1to2  rate at which indivs in pop1 migrate to pop2
#' @param rate2tp1 rate at which indivds in pop2 migrate to pop1
#' @return Returns a list of pop structs
#' @export
fast_dispersal <- function(pop1,pop2,rate1to2,rate2to1){
  n1 <- dim(pop1)[1]
  n2 <- dim(pop2)[1]
  L <- dim(pop1)[2]
  g <- 2 #diploidia?
  p1 <- ifelse(runif(n1)<rate1to2,2,1)
  p2 <- ifelse(runif(n2)<rate2to1,1,2)

  ret <- rcpp_dispersal_placement(pop1, pop2, dim(pop1), dim(pop2), p1, p2);
  dim(ret[[1]]) <- c(sum(c(p1, p2)==1), L, g)
  dim(ret[[2]]) <- c(sum(c(p1, p2)==2), L, g)

  ret
}


######
fitness <- function(z,n,b0,b1,b2,b3,epsilon, n.loci){
  ze <- z[,1]    
  ng=n.loci
  a=b0
  b=b1
  c=b2 
  w <- round( a*exp(-((ze-b*ng)^2)/(2*(c*ng)^2) ) - b3*n + rnorm(n,0,epsilon) , 0)  
  w <- ifelse(w<0,0,w)
  return(w)
}

###
##########################################
#' disperal function
#' @param pop1 pop struct 1
#' @param pop2 pop struct 2
#' @param rate1to2  rate at which indivs in pop1 migrate to pop2
#' @param rate2tp1 rate at which indivds in pop2 migrate to pop1
#' @return Returns a list of pop structs
#' @export
fast_dispersal <- function(pop1,pop2,rate1to2,rate2to1){
  n1 <- dim(pop1)[1]
  n2 <- dim(pop2)[1]
  L <- dim(pop1)[2]
  g <- 2 #diploidia?
  p1 <- ifelse(runif(n1)<rate1to2,2,1)
  p2 <- ifelse(runif(n2)<rate2to1,1,2)

  ret <- rcpp_dispersal_placement(pop1, pop2, dim(pop1), dim(pop2), p1, p2);
  dim(ret[[1]]) <- c(sum(c(p1, p2)==1), L, g)
  dim(ret[[2]]) <- c(sum(c(p1, p2)==2), L, g)

  ret
}


######
#Example 1
doer.ex1 <- function(x){
  struct <- x[[1]]
  bvs <- x[[2]]
  nloc <- x[[3]]
  ve <- x[[4]]
  b0 <- x[[5]]
  b1 <- x[[6]]
  b2 <- x[[7]]
  b3 <- x[[8]]
  epsil <- x[[9]]
  theta <- x[[10]]
  add.loci<-x[[11]]
  add.pos<-x[[12]]
  
  struct.fit<-struct[, add.pos:(add.pos+(add.loci-1)),]
  
    
  zs             <- rcpp_g2p_map(struct.fit, dim(struct.fit), bvs, add.loci, ve)
  fits           <- cbind(zs,fit=fitness(zs,dim(zs)[1],b0,b1,b2,b3,epsil, add.loci))   
  pairs          <- mating(fits,struct)
  mums           <- pairs[[1]]
  dads           <- pairs[[2]]
  struct.rt      <- array(NA, dim(mums))

  struct.rt[,,1] <- rcpp_recombo_segregate(mums, dim(mums), theta)
  struct.rt[,,2] <- rcpp_recombo_segregate(dads, dim(dads), theta)  
    
  return(struct.rt)
}

##########
main.function.ex1<- function(Nsize,nloci,nalleles,Ve, Vd, time, fitness.param, migration.rate, recom.map, start.1, start.2, bv.for.alleles, add.loci, add.pos){

initial.population.size <- Nsize # NOTE: THIS WILL BECOME A VECTOR TWO NUMBERS WHEN WE EXTEND TO TWO POPULATIONS
n.loci                  <- nloci # how many linked genes we will deal with

loci.pos                <- 1:n.loci

n.alleles.per.locus     <- nalleles  
bv.for.alleles          <- bv.for.alleles

add.loci				    <- add.loci
add.pos					<- add.pos

V.e                     <- Ve # standard deviation in environmental component of the phenotype
V.d                     <- Vd # standard demographic variation
n.gens                  <- time # number of generations.

disp					<- migration.rate #rate of dispersal
theta					<- recom.map #recombination rate


struct.1				<- start.1
struct.2				<- start.2

param					<- fitness.param

param[1,2]

############################################
pb <- progress_bar$new(format = " Work in progress [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)

res <- list()
for (i in 1:n.gens){

	x1 <- list(struct.1,bv.for.alleles,n.loci,V.e,param[1,1],param[1,2],param[1,3],param[1,4],V.d, theta, add.loci, add.pos )
	x2 <- list(struct.2,bv.for.alleles,n.loci,V.e,param[2,1],param[2,2],param[2,3],param[2,4],V.d, theta, add.loci, add.pos )
  
  x <- list(x1,x2)

  out <- lapply(x, doer.ex1)  
  struct.1 <- out[[1]]
  struct.2 <- out[[2]]
  outd <- fast_dispersal(struct.1,struct.2,disp,disp)
  struct.1 <- outd[[1]]
  struct.2 <- outd[[2]]

  if (i %% n.gens==0) res<- list(struct.1,struct.2)
  pb$tick()
}

return(res)
}

#######
#Example 2
##########################
### Biallelic Mutation ###
##########################

mutate<-function(i){
  	um<-runif(1)
  	um <- ifelse(um>mutation.rate,0,1)
  	i <- if (um==0) i <- i else i <- ifelse(i==1,0,1)
	return(i)
	}

	
Mutation <- function(genotypes, mutation.rate){
  new.genotypes<-apply(genotypes, 2, mutate)
  return(new.genotypes)
}



######
doer.ex2 <- function(x){  

  struct <- x[[1]]
  sex.ratio <- x[[2]]
  mean.fitness <- x[[3]]  
  loci.pos <-  x[[4]]
  chromo_mb <- x[[5]]
  d.e <- x[[6]]
  crossover <- x[[7]]
  mutation.rate <- x[[8]]
  
  z <- as.data.frame(cbind(sex = rbinom(nrow(struct), 1, sex.ratio)+1)) #zs
  rownames(z)=rownames(struct)
  
  n = nrow(struct)
  z$fit <- round( pmax(rpois(n, lambda=mean.fitness) - d.e*n, 0) ) 
  
  pairs <- mating(z,struct)
  
  mums           <- pairs[[1]]
  dads           <- pairs[[2]]
  struct.rt      <- array(NA, dim(mums))
 
  
  struct.rt[,,1] <- rcpp_recombo_segregate_expo(mums, dim(mums), loci.pos, chromo_mb, crossover)
  struct.rt[,,2] <- rcpp_recombo_segregate_expo(dads, dim(dads), loci.pos, chromo_mb, crossover) 

  struct.rt[,,1] <- Mutation(struct.rt[,,1, drop = F], mutation.rate)
  struct.rt[,,2] <- Mutation(struct.rt[,,2, drop = F], mutation.rate)  
   
  return(struct.rt)
}



######
main.function.ex2<- function(start.1, start.2, sex.ratio, mean.fitness, time, migration.rate, loci.pos, chromo_mb, d.e, crossover, mutation.rate){

n.gens = time

struct.1 <- start.1
struct.2 <- start.2

sex.ratio <- sex.ratio
mean.fitness <- mean.fitness
  
loci.pos <-  loci.pos
chromo_mb <- chromo_mb

disp <- migration.rate
d.e <- d.e

crossover <- crossover
mutation.rate <- mutation.rate

############################################
pb <- progress_bar$new(format = " Work in progress [:bar] :percent eta: :eta",total =n.gens, clear = FALSE, width= 100)

res <- list()
for (i in 1:n.gens){

  	x1 <- list(struct.1, sex.ratio, mean.fitness, loci.pos, chromo_mb, d.e, crossover, mutation.rate)
	x2 <- list(struct.2, sex.ratio, mean.fitness, loci.pos, chromo_mb, d.e, crossover, mutation.rate)
  
  	x <- list(x1,x2)

	out <- lapply(x,doer.ex2)  
    struct.1 <- out[[1]]
    struct.2 <- out[[2]]

  	outd <- fast_dispersal(struct.1,struct.2,disp,disp)
  	struct.1 <- outd[[1]]
  	struct.2 <- outd[[2]]
	if (i %% n.gens==0) res<- list(struct.1,struct.2) 
 	 pb$tick()
	}

return(res)
}



##################
###    Fst     ###
### estimation ###
##################

########################################
library(pegas)
library(dplyr)
#' convert the two population "struct"s to a data frame of loci that the pegas package can use
#'
#' This should work with any number of alleles
#' @param P1 the pop struct, indexed by indiv, locus, gene copy, for pop 1
#' @param P2 the pop struct for pop 2
#' @export
struct2pegas <- function(P1, P2) {
  L <- dim(P1)[2]

  pop1 <- paste(P1[,,1], P1[,,2], sep = "/") %>%
    matrix(ncol = L) %>%
    as.data.frame(stringsAsFactors = FALSE)
  pop2 <- paste(P2[,,1], P2[,,2], sep = "/") %>%
    matrix(ncol = L) %>%
    as.data.frame(stringsAsFactors = FALSE)

  ret <- list("1" = pop1, "2" = pop2) %>%
    dplyr::bind_rows(.id = "population")

  names(ret)[-1] <- paste("locus", 1:(ncol(ret)-1), sep = "_")

  # write it out so we can read it into pegas
  temp <- tempfile()
  write.table(ret, sep = " ", file = temp, col.names = TRUE, quote = FALSE, row.names = FALSE)

  pegas::read.loci(temp)
}
#########
#' at each of the loci compute Fst
#'
#' uses pegas
#' @inheritParams struct2pegas
#' @export
fst_at_loci_with_pos <- function(P1, P2, pos) {
  dat <- struct2pegas(P1, P2)
  fst <- pegas::Fst(dat)
  dplyr::data_frame(locus = rownames(fst), pos = pos, Fst = fst[,"Fst"], Fit = fst[,"Fit"], Fis = fst[, "Fis"])

}
