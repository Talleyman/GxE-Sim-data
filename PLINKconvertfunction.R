#!/usr/bin/Rscript

convert <- function(dir, filename, genEffectRatio, data_source=c("Maize","Pearl Millet"), data_type=c("QxSxE","SxE")){
  
  #First, read in data from directory
  setwd(dir)
  data <- read.table(filename,sep=",",header=T)
  filename <- tools::file_path_sans_ext(filename)
  
  #Next, create PED file
  ped <- matrix(nrow=nrow(data),ncol=6)
  gen <- matrix(nrow=nrow(data),ncol=nrow(data)*2)
  ped[,1] <- c(rep(1,nrow(data)/2),rep(2,nrow(data)/2))
  if (data_source=="Maize"){
    ped[,2] <- data$TaxaOrder
  } else {
    ped[,2] <- c(1:nrow(data))
  }
  ped[,3] <- rep(0,nrow(ped))
  ped[,4] <- rep(0,nrow(ped))
  ped[,5] <- rep(3,nrow(ped))
  ped[,6] <- data$Phenotype
  if (data_source=="Pearl Millet"){
    gen[,1] <- ifelse(data$snpPresence==1, "A","B")
  } else {
    gen[,1] <- ifelse(data$SNPPresence==1, "A","B")
  }
  
  if (data_source=="Maize"){
    if (data_type=="QxSxE"){
      gen[,2] <- ifelse(data$envE1Presence==1,"A","B")
    } else {
      gen[,2] <- ifelse(data$Env1Presence==1,"A","B")
    } 
    } else if (data_source=="Pearl Millet") {
      if (data_type=="QxSxE"){
        gen[,2] <- ifelse(data$envE1Presence==1,"A","B")
      } else {
        gen[,2] <- ifelse(data$EnvPresence==1,"A","B")
      } 
  }
  if (data_type=="QxSxE"){
    gen[,3] <- ifelse(data$BackgroundPresence==1,"A","B")
  } else {
    gen[,3] <- sample(LETTERS[1:2],nrow(gen),0.5)
  }
  for (i in 4:ncol(gen)){
    gen[,i] <- sample(LETTERS[1:2],nrow(gen),0.5)
  }
  ped<-cbind(ped,gen)
  
  #Write remaining two mandatory files for analysis
  if (data_source=="Maize"){
    map <- data.frame(rep(1,nrow(data)),data$Taxa,rep(0,nrow(data)),(1:nrow(data)))
    pheno <- data.frame(c(rep(1,nrow(data)/2),rep(2,nrow(data)/2)),data$TaxaOrder,data$Phenotype)
  } else {
    map <- data.frame(rep(1,nrow(data)),c(1:nrow(data)),rep(0,nrow(data)),(1:nrow(data)))
    pheno <- data.frame(c(rep(1,nrow(data)/2),rep(2,nrow(data)/2)),c(1:nrow(data)),data$Phenotype)
  }
  
  #Generate the covariate file

  multicov <- function(x, y, z){
    covxy <- sd(x)*sd(y)*cor(x,y)
    covyz <- sd(y)*sd(z)*cor(y,z)
    covxz <- sd(x)*sd(z)*cor(x,z)
    final <- sqrt((covxz^2+covyz^2-2*covxy*covxz*covyz)/(1-covxy^2))
    return(final)
  }

  famID <- c(rep(1,nrow(data)/2),rep(2,nrow(data)/2))
  indID <- data$Taxa
  if (data_type=="SxE"){
    if (data_source=="Maize"){
      covarvals <- cov(data$SNPPresence, data$Env1Presence)+rnorm(length(famID),mean=0,sd=0.5)
    } else {
      covarvals <- cov(data$SNPPresence, data$EnvPresence)+rnorm(length(famID),mean=0,sd=0.5)
    }
  } else if (data_type=="QxSxE"){
    covarvals <- multicov(data$snpPresence, data$envE1Presence, data$BackgroundPresence)+rnorm(length(famID),mean=0,sd=0.5)
  }
  covar <- data.frame(famID, indID, covarvals)
  
  #Finally, write the known-truth file with randomly generated, normally distributed effects
  SNPs <- as.character(sample(data$Taxa,15))
  KTeffects <- rnorm(15, mean=0, sd=0.5)
  knowntruth <- data.frame(SNPs, KTeffects)
  
  print("Writing output files...")
  write.table(ped, paste(filename, "ped",sep="."),quote=F,row.names=F,col.names=F)
  write.table(map, paste(filename,"map",sep="."),quote=F,row.names=F,col.names=F)
  write.table(pheno, paste(filename, "pheno",sep="."),quote=F, row.names=F,col.names=F)
  write.table(covar, paste(filename, "covar",sep="."),quote=F, row.names=F,col.names=F)
  write.table(knowntruth, paste(filename, "ote",sep="."),quote=F,row.names=F,col.names=F)
  print("Done!")
}