#!/usr/bin/Rscript

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##
## This script takes the output of SeedVicious and then SNPs that were previously 
## downloaded and  finds which SNPs are in the target regions:
##
## usage: Rscript ./Find_SNPs_in_Targets.R chr#
##
## Output: (these are put into the individual chromosome folders)
##	1)	DB_SNPinMiRtargs_chr#.txt 
##		-- A database of the SNPs that are in the target region of miRNA
##	2)	DB_SNPinMiRtargs_revSNPs_chr#.txt
##		-- The SNPs that are in reverse complement format (will need to be 
##			dealt with seperately- for now at least)
##	3)	DB_SNPinMiRtargs_errorSNPs_chr#.txt 
##		-- The SNPs which have a different major allele from my reference

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

args<-commandArgs(T)

# If not running as bash script, be sure to set the following manually:
chrNum <- args[1]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** RETRIEVING SeedVicious DATA ***", stdout())

# Import the output of the SeedVicious parser:
pathSeedV<- paste("./chr", chrNum ,"/Parser_seedV_miTall_chr", chrNum , 
						 "_biomartTargs_all.txt", sep = "")

seedVic<-read.table(pathSeedV, sep = "\t")					

# for test: seedVic <- read.table("~/Desktop/Parser_seedV_miRall_chr21_biomartTargs_p1.txt", sep = "\t")		
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** PREPARING SeedVicious DATA ***", stdout())
# Set the sequence up more usefully and name the columns:

# Extract sequence from surrounding descriptors in seedV output
seedVerboseParsed<-seedVic
seedVic[,6] <- unlist(lapply(as.vector(seedVerboseParsed[,6]), 
						function(x) strsplit(x, split = " ")[[1]][4]))
seedVic[,8] <- unlist(lapply(as.vector(seedVerboseParsed[,8]), 
						function(x) strsplit(x, split = " ")[[1]][5]))

# remove empty columns
seedVic[,5]<-NULL
seedVic[,6]<-NULL

# add column names
colnames(seedVic)<- c("tr", "pos", "miR", "type", "miRseq:3-5", "mRNAseq:5-3")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# There are some target UTRs that have multiple locations, these will be removed 
# for now - figure out what to do with these...

indexSingle <- grep(";", seedVic[,"tr"], fixed = T, invert = T)
single_loc <- seedVic[indexSingle,]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Now separate the dataset into 'near Target' and regular and add a column that
# indicates whether it is a target or a near target (and remove the * from the 
# targ type)


# first take the real targets and add a column indicating their group
indexT <- grep("*", single_loc[,"type"], fixed = T, invert = T)
targ_or_nt <- rep("targ", times = length(indexT))
targseedV <- cbind(single_loc[indexT,], targ_or_nt)

# now do the same with near targets (also removing the * from the type)
indexNt <- grep("*", single_loc[,"type"], fixed = T)
nearTseedVtemp <- single_loc[indexNt,] 
type <- gsub("*", replacement = "", nearTseedVtemp[,"type"], fixed = T)
targ_or_nt <- rep("*nt", times = length(indexNt))
nearTseedV <- cbind(nearTseedVtemp, targ_or_nt)
nearTseedV[,"type"]<-type

# Put them back together for the next part (and reset rownames)
seedV <- rbind(targseedV, nearTseedV)
row.names(seedV)<-NULL


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## Now find the target region (because seedV just says where in the transcript it 
## is, not overall) 

# As loops cannot be done with this size of data, I will seperate into groups and 
# do them invidually then recombine.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### SECTION 1: subset of only "R" (reverse complement) targets

rcSeedV<-seedV[grepl("R", as.vector(seedV[,"tr"])),]

# Get the genomic positions of the target parts
rcUTRstart <- unlist(lapply(as.vector(rcSeedV[,"tr"]), 
					function(x) as.numeric(strsplit(x , "|", fixed = T)[[1]][4])))
rcTARG_seed_start <- as.numeric(sub("\\).*", "", 
								sub(".*\\(", "", as.vector(rcSeedV[,"pos"]))))
# make seed region only the base 6nt inciated in bartel2009 for all targ sites
rcTARG_start <- rcUTRstart + rcTARG_seed_start
rcTARG_end <- rcUTRstart + rcTARG_seed_start + 5

# extract the rest of the relevant info direclty from the subset
rcmiR <- as.vector(rcSeedV[,"miR"])
rcmRNA <- as.vector(rcSeedV[,"tr"])
rctype <- as.vector(rcSeedV[,"type"])
rctarg_or_nt <- as.vector(rcSeedV[,"targ_or_nt"])
rcmiRseq.3.5 <- as.vector(rcSeedV[,"miRseq:3-5"])
rcmRNAseq.5.3 <- as.vector(rcSeedV[,"mRNAseq:5-3"])

# now put them all together and rename for combining below
rcTargInter <- data.frame(rcmiR, rcmRNA, rcTARG_start, rcTARG_end, rctype, 
						  rctarg_or_nt, rcmiRseq.3.5, rcmRNAseq.5.3)
colnames(rcTargInter) <- c("miR", "mRNA", "TARG_start", "TARG_end", "type", 
						   "targ_or_nt", "miRseq.3.5", "mRNAseq.5.3")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### SECTION 2: subset of all the rest of the targets (not reverse complement)

nSeedV<-seedV[!grepl("R", as.vector(seedV[,"tr"])),]

# Get the genomic positions of the target parts
nUTRstart <- unlist(lapply(as.vector(nSeedV[,"tr"]), 
					function(x) as.numeric(strsplit(x , "|", fixed = T)[[1]][4])))
nTARG_seed_start <- as.numeric(as.vector(nSeedV[,"pos"]))
	# -1 in next two to account for first position index (subtraction stuff...)
# make seed region only the base 6nt inciated in bartel2009 for all targ sites
nTARG_end <- nUTRstart - 2 + nTARG_seed_start
nTARG_start <- nUTRstart - 2 + nTARG_seed_start - 5

# extract the rest of the relevant info direclty from the subset
nmiR <- as.vector(nSeedV[,"miR"])
nmRNA <- as.vector(nSeedV[,"tr"])
ntype <- as.vector(nSeedV[,"type"])
ntarg_or_nt <- as.vector(nSeedV[,"targ_or_nt"])
nmiRseq.3.5 <- as.vector(nSeedV[,"miRseq:3-5"])
nmRNAseq.5.3 <- as.vector(nSeedV[,"mRNAseq:5-3"])

# now put them all together and rename for combining below
nTargInter <- data.frame(nmiR, nmRNA, nTARG_start, nTARG_end, ntype, ntarg_or_nt, 
						 nmiRseq.3.5, nmRNAseq.5.3)
colnames(nTargInter) <- c("miR", "mRNA", "TARG_start", "TARG_end", "type", 
						  "targ_or_nt", "miRseq.3.5", "mRNAseq.5.3")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# now combine the two subsets
targInter <- rbind(rcTargInter, nTargInter)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Now load the biomaRt SNPs file for this chromosome:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** LOADING SNP DATA ***", stdout())

# create the file name based on chromosome number
SNPs_output <- paste("./chr", chrNum ,"/SNPdb_BiomaRt_chr", chrNum , ".txt", 
					 sep = "")

# for test: SNPs_output <- "~/Desktop/SNPdb_BiomaRt_chr21.txt"

snp_db <- read.table(SNPs_output, header = T, sep = "\t")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Remove all the indels:
snp_db_indel <- snp_db[!grepl("-", snp_db[,"allele"]),]

# Now take only those with just one position change
snp_db_single <- snp_db_indel[!grepl("[ATCG][ATCG]", snp_db_indel[,"allele"]),]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# See if there are any reverse SNPs, justmake a list of them so can check back 
# later (this seems to be exceedingly uncommon). Leave them in the regular db for 
# the analysis though.
revSNPs <- snp_db_single[snp_db_single[,"chrom_strand"] == "-1",]
onlySNPs <- snp_db_single

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Take the SNPs that are located within one of the target regions 
# Thanks to this page: 
# http://stackoverflow.com/questions/3916195/finding-overlap-in-ranges-with-r
# source("http://bioconductor.org/biocLite.R")  # biocLite("IRanges")

library(IRanges)

write("*** FINDING SNPS THAT ARE IN TARGET REGIONS ***", stdout())

# Set up the ranges of the SNPs:
SNPs_ranges <- IRanges(names = c(onlySNPs$refsnp_id, "rsTEST"), 
					  c(onlySNPs$chrom_start, 10), c(onlySNPs$chrom_end, 10))

# Set up the ranges of the target regions:
targ_ranges <- IRanges(c(targInter$TARG_start, 9), c(targInter$TARG_end, 15))

# Find the 'interactions' between SNPs and target regions: order = query, subject
over <- findOverlaps(SNPs_ranges, targ_ranges, type = "within")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Now get the info which corresponds to the SNP/targ interactions in 'over'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# set all variables to empty, just in case used above
miRNA <- c()
utr3 <- c()
targ_start <- c()
targ_end <- c()
chrom <- c()
snpName <- c()
snpPos <- c()
snpDB_alleles <- c()
utr_range <- c()
snp_minorAllele <- c()
snp_minorAlleleFreq <- c()
targ_or_nt <- c()
miRseq.3.5 <- c()
mRNAseq.5.3 <- c()
utrENSG <- c()
type <- c()
ancestAllele <- c()

### First extract info from the 'subjectHits' aka 'targInter' ###
miRNA <- as.vector(targInter[subjectHits(over)[],"miR"])
targ_start <- as.vector(targInter[subjectHits(over)[],"TARG_start"])
targ_end <- as.vector(targInter[subjectHits(over)[],"TARG_end"])
targ_or_nt <- as.vector(targInter[subjectHits(over)[],"targ_or_nt"])
type <- as.vector(targInter[subjectHits(over)[],"type"])
miRseq.3.5 <- as.vector(targInter[subjectHits(over)[],"miRseq.3.5"])
mRNAseq.5.3 <- as.vector(targInter[subjectHits(over)[],"mRNAseq.5.3"])


### Then extract info from the 'queryHits' aka 'onlySNPs' ###
snpName <- as.vector(onlySNPs[queryHits(over)[],"refsnp_id"])
snpPos <- as.vector(onlySNPs[queryHits(over)[],"chrom_start"])
snpDB_alleles <- as.vector(onlySNPs[queryHits(over)[],"allele"])
snp_minorAllele <- as.vector(onlySNPs[queryHits(over)[],"minor_allele"])
snp_minorAlleleFreq <- as.vector(onlySNPs[queryHits(over)[],"minor_allele_freq"])
ancestAllele <- as.vector(onlySNPs[queryHits(over)[],"allele_1"])


### Now take info about/from the UTR ###
# extract the name so I dont have to do it in each part
utr3_name <- as.vector(targInter[subjectHits(over)[], "mRNA"])
	
# to only take parts of the utr name:
utr3 <- unlist(lapply(utr3_name, function(x) 
				paste(strsplit(x, "|", fixed = T)[[1]][1], "_", 
				strsplit(x, "|", fixed = T)[[1]][2], sep = "")))
				
# get the utr range out of the name as well (will need this later)	
utr_range <- unlist(lapply(utr3_name, function(x) 
					paste(strsplit(x, "|", fixed = T)[[1]][4], ":", 
					strsplit(x, "|", fixed = T)[[1]][5], sep = "")))
	
# get the chromosome number:
chrom <- unlist(lapply(utr3_name, function(x) 
				strsplit(x, "|", fixed = T)[[1]][3]))
	
# to get just the ENSG - for removing duplicates (below)
utrENSG <- unlist(lapply(utr3_name, function(x) 
				  strsplit(x, "|", fixed = T)[[1]][1]))
				  
### Now put it all together ###
DB_interactions_temp <- data.frame(miRNA, utr3, targ_start, targ_end, chrom, 
								   snpName, snpPos, snpDB_alleles, utr_range, 
								   snp_minorAllele, snp_minorAlleleFreq, 
								   ancestAllele, targ_or_nt, type, miRseq.3.5, 
								   mRNAseq.5.3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Remove duplicated SNP/ranges (because of multiple transcripts of one gene)
# the function duplicated() can be used to find the colums that are duplicated,
# and lucky for me, it only marks as 'TRUE' for the second one (so first 
# duplicate). Therefore I can just take all the 'FALSE' ones :D

# I will consider duplicates to be the same in miR, ENSG (NOT ENST), targ start, 
# targ end, SNP name, and SNP position

repTEST <- data.frame(miRNA, utrENSG, targ_start, targ_end, chrom, snpName, 
						snpPos)

DB_interactions <- DB_interactions_temp[!duplicated(repTEST),] 

# To fix the rownames to the adjusted rows:
rownames(DB_interactions)<-NULL

# if last is NA, remove it (it causes probs later on): tail(DB_interactions)
if(is.na(as.vector(DB_interactions[length(DB_interactions[,1]),1]))){
	DB_interactions <- DB_interactions[-length(DB_interactions[,1]),]
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Because the cases are different and I cannot use loops, I will seperate my 
# database into four groups -- First by spliting into 'real target' and 'near 
# target' then each of those subsequently into reverse and forward)

# Just the 'real targets' (PART A)					
DB_interactions_targ <- DB_interactions[!grepl("*nt", 
									DB_interactions[,"targ_or_nt"], fixed = T),]


# Just the near target interactions (PART B)
DB_interactions_nt <- DB_interactions[grepl("*nt", DB_interactions[,"targ_or_nt"], 
						fixed = T),]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# library(Biostrings)
# Because DNAstring/RNAstring from Biostrings takes so long and dont really need 
# it, I will make my own revComp function. And, so I don't have to turn the Us to 
# Ts first, this accounts for that (also I dont want any Us in reverse complements)

revComp <- function(x){
	x2<-x
	x2[x=="A"] <- "T"
	x2[x=="T"] <- "A"
	x2[x=="U"] <- "A"
	x2[x=="C"] <- "G"
	x2[x=="G"] <- "C"
	return(x2)
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Find nt at the SNP position (for the real targets)


#### PART A: SECTION 1: Real target, reverse complement #####

# First extract the reverse complement accessions
rcDB_interactions_targ <- DB_interactions_targ[grepl("R", 
								 as.vector(DB_interactions_targ[,"utr_range"])),]
								 								 				 
# Get the position of the SNP in the target sequence (from the seedVparser)
rcmRseqPos <- rcDB_interactions_targ[,"snpPos"] - 
				rcDB_interactions_targ[,"targ_start"] + 1

# Get the position in the seq (starting from the beginning of the seq -- in the 
# seedVparser output it is in the last potion)
rcmRseqPos2 <- nchar(as.vector(rcDB_interactions_targ[,"mRNAseq.5.3"])) - 
					 rcmRseqPos

# Take the nucleotide at that position:
rcrefNT <- revComp(substr(as.vector(rcDB_interactions_targ[,"mRNAseq.5.3"]), 
					rcmRseqPos2, rcmRseqPos2))

# As they are the same (since this is a real target) just make them the same 	
rctargNT <- rcrefNT

rcDB_interactions_targ_ref <- cbind(rcDB_interactions_targ, rctargNT, rcrefNT)

# Give it the general column names for combination below								
colnames(rcDB_interactions_targ_ref)<-c(colnames(DB_interactions_nt), 
									    "targNT", "refNT")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

##### PART A: SECTION 2: Real target, 'normal' #####

# Extract the normal (not rev. comp) accessions
nDB_interactions_targ <- DB_interactions_targ[!grepl("R", 
								as.vector(DB_interactions_targ[,"utr_range"])),]

# Get the position of the SNP in the target sequence (from the seedVparser)
nmRseqPos <- nDB_interactions_targ[,"targ_end"] - nDB_interactions_targ[,"snpPos"] + 1

# Get the position in the seq (starting from the beginning of the seq -- in the 
# seedVparser output it is in the last potion)
nmRseqPos2 <- nchar(as.vector(nDB_interactions_targ[,"mRNAseq.5.3"])) - nmRseqPos

# Take the nucleotide at that position:
nrefNT <- substr(as.vector(nDB_interactions_targ[,"mRNAseq.5.3"]), nmRseqPos2, 
				  nmRseqPos2)
nrefNT[nrefNT == "U"] <- "T"	

# As they are the same (since this is a real target) just make them the same 
ntargNT <- nrefNT

nDB_interactions_targ_ref <- cbind(nDB_interactions_targ, ntargNT, nrefNT)

# Give it the general column names for combination below									  
colnames(nDB_interactions_targ_ref)<-c(colnames(DB_interactions_nt), 
									    "targNT", "refNT")
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Now put PART A back together:
DB_interactions_targ_partA <- rbind(rcDB_interactions_targ_ref, 
									nDB_interactions_targ_ref)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 
# Now work with the near targets to see which are actually interesting - meaning
# which have a SNP in a near targ change postion that makes them real targets

# First I need to take the position of the 'near target' location, and what it is 
# changing from and to

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** FINDING 'NEAR TARGETS' WITH SNPS MAKING THEM REAL TARGETS ***", 
		stdout())

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#### PART B: SECTION 1: Reverse complement #####

# First extract the reverse complement accessions
rcDB_interactions_nt <- DB_interactions_nt[grepl("R", 
								 as.vector(DB_interactions_nt[,"utr_range"])),]

# Then find the near target mRNA and miRNA seq from seedV output
rcmRseq <- as.vector(rcDB_interactions_nt[, "mRNAseq.5.3"])
rcmiRseq <- as.vector(rcDB_interactions_nt[, "miRseq.3.5"])	

# Find the potion of the nt change in the seedV seq output
rcntPos_inSeq <- unlist(lapply(rcmRseq, function(x) 
						gregexpr(pattern = "[a,u,c,g]", x)[[1]][1]))
						
# Find the genomic position of the nt change
rcntPos <- rcDB_interactions_nt[,"targ_start"] + (nchar(rcmRseq) - rcntPos_inSeq - 1)

# Take the nt of the original seq
rcmR_NTchange <- substr(rcmRseq, rcntPos_inSeq, rcntPos_inSeq)
rcntFrom <- revComp(toupper(rcmR_NTchange))			

# Take the nt that it needs to change to to be a target
rcmiR_NTchange <- substr(rcmiRseq, rcntPos_inSeq, rcntPos_inSeq)
rcntTo <- revComp(toupper(rcmiR_NTchange))
			
# Add the new info to the database	
rcDB_interactions_ntdiff <- data.frame(rcDB_interactions_nt, rcntPos, rcntFrom, 
									   rcntTo)
									   
# Give it the general column names for combination below		   
colnames(rcDB_interactions_ntdiff)<-c(colnames(DB_interactions_nt), "ntPos", 
									  "ntFrom", "ntTo")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#### PART B: SECTION 2: Normal (not reverse complement) #####

# Extract the normal (not rev. comp) accessions
nDB_interactions_nt <- DB_interactions_nt[!grepl("R", 
								as.vector(DB_interactions_nt[,"utr_range"])),]

# Then find the near target mRNA and miRNA seq from seedV output
nmRseq <- as.vector(nDB_interactions_nt[, "mRNAseq.5.3"])
nmiRseq <- as.vector(nDB_interactions_nt[, "miRseq.3.5"])

# Find the potion of the nt change in the seedV seq output
nntPos_inSeq <- unlist(lapply(nmRseq, function(x) 
					   gregexpr(pattern = "[a,u,c,g]", x)[[1]][1]))	

# Find the genomic position of the nt change
nntPos <- nDB_interactions_nt[, "targ_end"] - (nchar(nmRseq) - nntPos_inSeq - 1)

# Take the nt of the original seq (and change Us to Ts)
nmR_NTchange <- substr(nmRseq, nntPos_inSeq, nntPos_inSeq)
nntFrom <- toupper(nmR_NTchange)   
nntFrom[nntFrom == "U"] <- "T"   
				   
# Take the nt that it needs to change to to be a target
nmiR_NTchange <- substr(nmiRseq, nntPos_inSeq, nntPos_inSeq)	
nntTo <- revComp(toupper(nmiR_NTchange))
	
# Add the new info to the database	
nDB_interactions_ntdiff <- data.frame(nDB_interactions_nt, nntPos, nntFrom, 
									  nntTo)
 
# Give it the general column names for combination below  
colnames(nDB_interactions_ntdiff)<-c(colnames(DB_interactions_nt), "ntPos", 
									 "ntFrom", "ntTo")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Put PART B back together (but not with PART A until the next part is done)
DB_interactions_ntdiff <- rbind(rcDB_interactions_ntdiff, nDB_interactions_ntdiff)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Now find which SNPs are in the near target locations
possibles <- DB_interactions_ntdiff[DB_interactions_ntdiff[, "snpPos"] == 
									DB_interactions_ntdiff[, "ntPos"],]
rownames(possibles) <- NULL

# First confirm that the major allele is the one from my UTR seq.
# These may be because of databases with different reference sequences or errors 
# or reverse SNPs which I will have to deal with (although I have not yet found 
# any, they are saved here: revSNPs)

ErrorSNPS <- possibles[unlist(lapply(as.vector(possibles[,"snpDB_alleles"]), 
						function(x) ((strsplit(x, split = "/")[[1]][1]))) != 
						as.vector(possibles[,"ntFrom"])),]

# To avoid a loop, I am making three logical vectors to check the three possible 
# positions of the alleles (ex:C/T/A/G) Then combining them and turning NAs to Fs
# and the extracting the Ts from possibles. 

pos1 <- unlist(lapply(as.vector(possibles[,"snpDB_alleles"]), function(x) 
		((strsplit(x, split = "/")[[1]][2]))) == as.vector(possibles[,"ntTo"]))
pos2 <- unlist(lapply(as.vector(possibles[,"snpDB_alleles"]), function(x) 
		((strsplit(x, split = "/")[[1]][3]))) == as.vector(possibles[,"ntTo"]))
pos3 <- unlist(lapply(as.vector(possibles[,"snpDB_alleles"]), function(x) 
		((strsplit(x, split = "/")[[1]][4]))) == as.vector(possibles[,"ntTo"]))
allPos <- pos1 | pos2 | pos3
allPos[is.na(allPos)]<- FALSE

nt_SNP_to_targ <- possibles[allPos,]
# nt_rejects <- possibles[!allPos,]

rownames(nt_SNP_to_targ) <- NULL
# rownames(nt_rejects) <- NULL

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Reorganize the columns to match with PART A

DB_interactions_neartarg_partB <- cbind(nt_SNP_to_targ[,1:16], 
										nt_SNP_to_targ[,"ntTo"], 
										nt_SNP_to_targ[,"ntFrom"])
colnames(DB_interactions_neartarg_partB) <- c(colnames(nt_SNP_to_targ[,1:16]), 
												"targNT", "refNT")

# Connect the two parts:

DB_interactions_withNT <- rbind(DB_interactions_targ_partA, 
								DB_interactions_neartarg_partB)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


### FOR BASH SCRIPT OUTPUT: ###

# The output files will be:
# 1) The SNP/target interactions: DB_interactions_withNT
# 2) The 'mistakes' file, with any reverse SNPs (that can be dealt with if they 
#		exist): revSNPs
# 3) The SNPs which have a different major allele from my reference: ErrorSNPS		

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Create the file names based on chromosome number
SNPtargs_output <- paste("./chr", chrNum ,"/DB_SNPinMiRtargs_chr", chrNum , 
						 ".txt", sep = "")

SNPtargs_output_revSNPs <- paste("./chr", chrNum ,
					"/DB_SNPinMiRtargs_revSNPs_chr", chrNum , ".txt", sep = "")
					
SNPtargs_output_errorSNPs <- paste("./chr", chrNum ,
					"/DB_SNPinMiRtargs_errorSNPs_chr", chrNum , ".txt", sep = "")
					
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# SAVE the database and error files:
write.table(DB_interactions_withNT, file = SNPtargs_output, 
			sep = "\t", row.names = F, quote = F)
			
write.table(revSNPs, file = SNPtargs_output_revSNPs, 
			sep = "\t", row.names = F, quote = F)
			
write.table(ErrorSNPS, file = SNPtargs_output_errorSNPs, 
			sep = "\t", row.names = F, quote = F)

write("*** DONE: CHECK OUTPUT FILES ***", stdout())

# for test: write(SNPs, file = "~/Desktop/data_SNPs_in_miR_targets/chr21/DB_SNPinMiRtargs_chr21.txt")


















	
