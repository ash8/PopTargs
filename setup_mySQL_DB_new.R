#!/usr/bin/Rscript
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# This script Takes the big database of all the interactions (per chromosome) and 
# separates it into the tables that will go into the mySQL database.

# usage: Rscript ./setup_mySQL_DB.R chr#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Take the chromosome number as an argument

args<-commandArgs(T)
chrNum <- args[1]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** RETRIEVING INTERACTIONS DATABASE ***", stdout())

# Import the big database of all interactions (for the chromosome)

pathInterDB <- paste("./chr", chrNum ,"/DB_SNPinMiRtargs_chr", chrNum , ".txt", 
						sep = "")
InterDB <- read.table(pathInterDB, header = T, sep = "\t")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** EXTRACTING INTERACTION INFO ***", stdout())

# Take portions of each of these components
miRname <- unlist(lapply(as.vector(InterDB[,"miRNA"]), function(x) 
					strsplit(x, " ")[[1]][1]))
ensg <- unlist(lapply(as.vector(InterDB[,"utr3"]), function(x) 
                    strsplit(x, "_")[[1]][1]))
enst <- unlist(lapply(as.vector(InterDB[,"utr3"]), function(x) 
					strsplit(x, "_")[[1]][2]))

# Record where target nt is the same as the reference nt
targNT <- as.vector(InterDB[,"targNT"])
refNT <- as.vector(InterDB[,"refNT"])
targ_is_ref <- as.numeric(targNT == refNT)

# Take the sequences as well (changing name)
miRseq_3_5 <- as.vector(InterDB[,"miRseq.3.5"])
mRNAseq_5_3 <- as.vector(InterDB[,"mRNAseq.5.3"])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Make empty columns for the allele frequencies (put all the possible frequencies 
# as factors so that the needed ones are able to be added later)
freqA <- factor(rep("NA", times = length(InterDB[,1])), 
                levels = unique(c(unique(InterDB[,"snp_minorAlleleFreq"]), 
                (1 - unique(InterDB[,"snp_minorAlleleFreq"])), "NA")))
freqT <- factor(rep("NA", times = length(InterDB[,1])), 
                levels = unique(c(unique(InterDB[,"snp_minorAlleleFreq"]), 
                (1 - unique(InterDB[,"snp_minorAlleleFreq"])), "NA")))
freqC <- factor(rep("NA", times = length(InterDB[,1])), 
                levels = unique(c(unique(InterDB[,"snp_minorAlleleFreq"]), 
                (1 - unique(InterDB[,"snp_minorAlleleFreq"])), "NA")))
freqG <- factor(rep("NA", times = length(InterDB[,1])), 
                levels = unique(c(unique(InterDB[,"snp_minorAlleleFreq"]), 
                (1 - unique(InterDB[,"snp_minorAlleleFreq"])), "NA")))
                
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Create the database:
# Take the important columns and add on the extras to fill later:
              
InterDB_snp <- cbind(miRname, ensg, enst, chrNum, 
                  InterDB[ ,c("utr_range", "snpName", "snpPos", "snpDB_alleles", 
                  "type", "targ_start", "targ_end", "targNT", "refNT", 
                  "ancestAllele", "snp_minorAllele", "snp_minorAlleleFreq")], 
                  targ_is_ref, miRseq_3_5, mRNAseq_5_3, freqA, freqT, freqC, freqG)
          
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** GETTING ALLELE FREQUENCIES ***", stdout())

# Take all the accessions that have SNP allele frequencies and ancestral alleles:

boolFreqs <- (as.vector(InterDB_snp[,"snp_minorAllele"]) != "") & (as.vector(InterDB_snp[,"ancestAllele"]) != "")
InterDB_onlyFreqs_all <- InterDB_snp[boolFreqs,]


# Do the single allele alternative SNPs first 
# (Not sure yet how I will do the others...)

InterDB_onlyFreqs <- InterDB_onlyFreqs_all[!grepl("/[ATCG]/", 
                                  InterDB_onlyFreqs_all[,"snpDB_alleles"]), ]

# Take only the MINOR ALLELES and add its frequency to the table
# - sadly the "snpDB_alleles" is not always correct... work-around below

minorAllele <- as.vector(InterDB_onlyFreqs[,"snp_minorAllele"])

# For each allele, replace the "NA" in the frequency column with the minor allele
# The rownames in these limited sets are still the same from InterDB_snp, so can 
# change them there in the overall matrix

InterDB_onlyFreqs[rownames(InterDB_onlyFreqs[minorAllele == "A", ]),"freqA"] <- 
            InterDB_onlyFreqs[rownames(InterDB_onlyFreqs[minorAllele == "A", ]), 
            "snp_minorAlleleFreq"]
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs[minorAllele == "T",]),"freqT"] <- 
            InterDB_onlyFreqs[rownames(InterDB_onlyFreqs[minorAllele == "T",]),
            "snp_minorAlleleFreq"]
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs[minorAllele == "C",]),"freqC"] <- 
            InterDB_onlyFreqs[rownames(InterDB_onlyFreqs[minorAllele == "C",]),
            "snp_minorAlleleFreq"]
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs[minorAllele == "G",]),"freqG"] <- 
            InterDB_onlyFreqs[rownames(InterDB_onlyFreqs[minorAllele == "G",]),
            "snp_minorAlleleFreq"]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Back to the problem with the minor allele in "snp_minorAllele" not always being 
# the second listed in "snpDB_alleles". (btw it seems like this is cases where 
# the freq is close to 50/50 - therefore may be because with more data the minor 
# has changed...) To solve this I will take them seperately - the ones that do 
# match the pattern and those that do not (this is just to avoid loops)

minorAllele <- as.vector(InterDB_onlyFreqs[,"snp_minorAllele"])
minorAllele2 <- unlist(lapply(as.vector(InterDB_onlyFreqs[,"snpDB_alleles"]), 
						function(x) strsplit(x, split = "/")[[1]][2]))

# The two sets 'regular' and 'irregular'
InterDB_onlyFreqs_irreg <- InterDB_onlyFreqs[!minorAllele == minorAllele2,]
InterDB_onlyFreqs_reg <- InterDB_onlyFreqs[minorAllele == minorAllele2,]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Now find the MAJOR ALLELES and do the same:

# Start with the ones in 'regular' format:

majorAllele_reg <- unlist(lapply(as.vector(InterDB_onlyFreqs_reg[ ,
                           "snpDB_alleles"]), function(x) strsplit(x,split = "/")
                           [[1]][1]))

InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_reg[majorAllele_reg == "A", ]),"freqA"] <- 
            (1 - InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_reg[majorAllele_reg == 
            "A", ]), "snp_minorAlleleFreq"])
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_reg[majorAllele_reg == "T", ]),"freqT"] <- 
            (1 - InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_reg[majorAllele_reg == 
            "T", ]), "snp_minorAlleleFreq"])
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_reg[majorAllele_reg == "C", ]),"freqC"] <- 
            (1 - InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_reg[majorAllele_reg == 
            "C", ]), "snp_minorAlleleFreq"])
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_reg[majorAllele_reg == "G", ]),"freqG"] <- 
            (1 - InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_reg[majorAllele_reg == 
            "G", ]), "snp_minorAlleleFreq"])


# Then those in 'irregular' format:

majorAllele_irreg <- unlist(lapply(as.vector(InterDB_onlyFreqs_irreg
                            [,"snpDB_alleles"]), function(x) 
                            strsplit(x, split = "/")[[1]][2]))

InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_irreg[majorAllele_irreg == "A", ]), 
            "freqA"] <- (1 - InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_irreg
            [majorAllele_irreg == "A", ]), "snp_minorAlleleFreq"])
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_irreg[majorAllele_irreg == "T", ]), 
            "freqT"] <- (1 - InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_irreg
            [majorAllele_irreg == "T", ]), "snp_minorAlleleFreq"])
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_irreg[majorAllele_irreg == "C", ]), 
            "freqC"] <- (1 - InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_irreg
            [majorAllele_irreg == "C", ]), "snp_minorAlleleFreq"])
InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_irreg[majorAllele_irreg == "G", ]), 
            "freqG"] <- (1 - InterDB_onlyFreqs[rownames(InterDB_onlyFreqs_irreg
            [majorAllele_irreg == "G", ]), "snp_minorAlleleFreq"])

# re-set the rownames
rownames(InterDB_onlyFreqs) <- NULL

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Now get the derived target allele frequencies:

# Find the allele freqs of the 'target' and 'derived' alleles
targ_allele <- unlist(lapply(InterDB_onlyFreqs[,"targNT"], function(x) 
                      paste("freq", x, sep = "")))
                      
targ_allele_freq <- unlist(lapply(1:length(InterDB_onlyFreqs[,1]), function(x) as.numeric(as.vector(InterDB_onlyFreqs[ x , targ_allele[x]]))))
                           
derived_allele <- unlist(lapply(InterDB_onlyFreqs[,"ancestAllele"], 
                         function(x) paste("freq", x, sep = "")))
                         
derived_allele_freq <- unlist(lapply(1:length(InterDB_onlyFreqs[,1]), function(x)
                    1 - as.numeric(as.vector(InterDB_onlyFreqs[ x , derived_allele[x]]))))


InterDB_final <- cbind(InterDB_onlyFreqs, targ_allele_freq, derived_allele_freq)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Create the file names:

interactions_table <- paste("./chr", chrNum ,"/mysqlTable_complete_chr", 
                            chrNum , ".txt", sep = "")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Output the databases:

# I want rownames only for DB_interactions because that will be the table key 
write.table(InterDB_final, file = interactions_table, sep = "\t", row.names = F, 
            col.names = F, quote = F )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** DONE. CHECK OUTPUT FILES ***", stdout())

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



