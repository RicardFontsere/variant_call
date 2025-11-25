# SCRIPT FOR INSPECTING AND FILTERING READS CONTAINING SEX-ASSOCIATED KMERS

# load packages

library(tidyverse)
library(patchwork)

# Read in data

kmer_counts <- read.delim("/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/most_assoc_kmers.table")

###check Kmers###

#make male and female sums
femaleInds <- read.delim("/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/female_list.txt", header=F)
maleInds <- read.delim("/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/male_list.txt", header=F)

#for each sex, count how many times a kmer is found in a read per individual, also note if that kmer is not found in other samples of same sex

#now females

female_files <- list.files(path = "/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/fbbduk", pattern = "_stats")
DF_Females <-  read.delim(paste0("/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/fbbduk/",female_files[1]), header=FALSE, comment.char="#")
colnames(DF_Females)<-c("KmerName",female_files[1],"perc")
DF_Females$perc<-NULL

for (i in female_files[-1]){
  df <- read.delim(paste0("/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/fbbduk/",i), header=FALSE, comment.char="#")      # read the file
  colnames(df)<-c("KmerName",i,"perc")
  print(head(df))
  df$perc<-NULL
  DF_Females <- merge(DF_Females, df, all=T, by="KmerName")    # append the current file
}

rownames(DF_Females)<-DF_Females$KmerName
DF_Females$KmerName<-NULL
##replace NA with 0 since it means no match
DF_Females[is.na(DF_Females)] <- 0

DF_Females$countszero<-rowSums(DF_Females == 0)
table(DF_Females$countszero)

#now males

male_files <- list.files(path = "/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/mbbduk", pattern = "_stats")
DF_Males <-  read.delim(paste0("/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/mbbduk/",male_files[1]), header=FALSE, comment.char="#")
colnames(DF_Males)<-c("KmerName",male_files[1],"perc")
DF_Males$perc<-NULL

for (i in male_files[-1]){
  df <- read.delim(paste0("/user/brussel/109/vsc10945/home/scratch/SNPs/results/PL/kmerGWAS/mbbduk/",i), header=FALSE, comment.char="#")      # read the file
  colnames(df)<-c("KmerName",i,"perc")
  print(head(df))
  df$perc<-NULL
  DF_Males <- merge(DF_Males, df, all=T, by="KmerName")    # append the current file
}

rownames(DF_Males)<-DF_Males$KmerName
DF_Males$KmerName<-NULL
##replace NA with 0 since it means no match
DF_Males[is.na(DF_Males)] <- 0

DF_Males$countszero<-rowSums(DF_Males == 0)
table(DF_Males$countszero)


####get kmer match statistics to see if we want to downfilter the kmers and the corresponding assembly

########## FILTERING #################

#### edit this for the numbers of samples

DF_Males$totalread<-rowSums(DF_Males[,1:3]) #Check number of males/females samples
DF_Females$totalread<-rowSums(DF_Females[,1:3]) #Check number of males/females samples

DF_combKmer<-data.frame(DF_Females$totalread, DF_Females$countszero)
DF_combKmer$kmer<-row.names(DF_Females)
DF_temp<-data.frame(DF_Males$totalread, DF_Males$countszero)
DF_temp$kmer<-row.names(DF_Males)
DF_combKmer<-merge(DF_combKmer,DF_temp, all.x = T, all.y=T)

DF_combKmer$DF_Females.totalread[is.na(DF_combKmer$DF_Females.totalread)] <- 0
DF_combKmer$DF_Males.totalread[is.na(DF_combKmer$DF_Males.totalread)] <- 0
DF_combKmer$DF_Females.countszero[is.na(DF_combKmer$DF_Females.countszero)] <- 3
DF_combKmer$DF_Males.countszero[is.na(DF_combKmer$DF_Males.countszero)] <- 3

##remove kmers which are present in all individuals

DF_combKmer<-subset(DF_combKmer, !(DF_combKmer$DF_Females.countszero==0 &  DF_combKmer$DF_Males.countszero==0))

### focus on each sex then, define female specific as present in all but maximum not in 3 females and in none or maximum 3 males

### CHANGE THESE VALUES DEPENDING ON YOUR NUMBER OF SAMPLES ###
DF_fem_specific<-subset(DF_combKmer, DF_combKmer$DF_Females.countszero==0 & DF_combKmer$DF_Males.countszero==3)
DF_male_specific<-subset(DF_combKmer, DF_combKmer$DF_Males.countszero==0 & DF_combKmer$DF_Females.countszero==3)

### write output files, to be used in read extraction and assembly

write.table(DF_fem_specific, file = paste0(SPECIES, "/FilteredFemaleKmers.txt"), sep="\t", quote=F, row.names = F)
write.table(DF_male_specific, file = paste0(SPECIES, "/FilteredMaleKmers.txt"), sep="\t", quote=F, row.names = F)

## go to cluster
## get the female ones from the female reads only and the male ones from the males
## assemble and reblast



