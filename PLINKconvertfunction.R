#!/usr/bin/Rscript

convert <- function(dir, filename, genEffectRatio){
  #First, read in data from directory
  setwd(dir)
  data <- read.table(filename,sep=",",header=T)
  filename <- tools::file_path_sans_ext(filename)
  
  #Next, create PED file
  ped <- matrix(nrow=nrow(data),ncol=6)
  gen <- matrix(nrow=nrow(data),ncol=nrow(data)*2)
  ped[,1] <- c(rep(1,276),rep(2,276))
  ped[,2] <- data$TaxaOrder
  ped[,3] <- rep(0,nrow(ped))
  ped[,4] <- rep(0,nrow(ped))
  ped[,5] <- rep(3,nrow(ped))
  ped[,6] <- data$Phenotype
  gen[,1] <- ifelse(data$snpPresence==1, "A","B")
  gen[,2] <- ifelse(data$envE1Presence==1,"A","B")
  gen[,3] <- ifelse(data$BackgroundPresence==1,"A","B")
  for (i in 4:ncol(gen)){
    gen[,i] <- sample(LETTERS[1:2],nrow(gen),0.5)
  }
  ped<-cbind(ped,gen)
  
  #Write remaining two files for analysis
  map <- data.frame(rep(1,nrow(data)),data$Taxa,rep(0,nrow(data)),(1:nrow(data)))
  pheno <- data.frame(c(rep(1,276),rep(2,276)),data$TaxaOrder,data$Phenotype)
  
  #Finally, write the known-truth file with randomly generated, normally distributed effects
  SNPs <- as.character(sample(data$Taxa,15))
  KTeffects <- rnorm(15, mean=0, sd=0.5)
  knowntruth <- data.frame(SNPs, KTeffects)
  
  write.table(ped, paste(filename, "ped",sep="."),quote=F,row.names=F,col.names=F)
  write.table(map, paste(filename,"map",sep="."),quote=F,row.names=F,col.names=F)
  write.table(pheno, paste(filename, "pheno",sep="."),quote=F, row.names=F,col.names=F)
  write.table(knowntruth, paste(filename, "ote",sep="."),quote=F,row.names=F,col.names=F)
}