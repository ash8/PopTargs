#!/usr/bin/Rscript

# Script for taking in all the 1000 genomes snp info 

# usage: Rscript ./PopData_2ndRound.R 3

# - - - - - - - - - - - 


# Take in the part number:
args<-commandArgs(T)
chrNum <- args[1] 

write("*** Starting import of the SNP and population data for the chromosome  ***", stdout())

# Load the data by chr
fileNameChr <- paste("PopData_Chr_", chrNum, ".txt", sep = "")

# ### for rev comp!
# fileNameChr <- paste("revComp_PopData_Chr_", chrNum, ".txt", sep = "")

all_SNPpop_freqs_T <- read.table(fileNameChr, header = T)


# Load the full interactions DB
 allDB <- read.table("CompleteSetALLINFO_jan18.txt", header = T)   

# ### for rev comp
# allDB <- read.table("CompleteSetALLINFO_plusRevComp.txt", header = T) 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


### NOTE: The 'AF' in the 1000 genomes (all_SNPpop_freqs) is the MINOR ALLELE FREQUENCY!!! NOT THE REF!!!


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# this section is all about speeding up the big loop by removing all extraneous data

# make a simple small version of the allDB: Take only those in this chr and only some columns
allDB_chr_T <- unique(allDB[allDB[ ,"chrNum"] == chrNum, c("chrNum", "miRname", "ensg", "enst", "snpName", "targ_is_ref", "ancestralAllele", "targNT", "refNT")])


# and only those SNPS in the populations doc
allDB_chr_T2 <- allDB_chr_T[allDB_chr_T[ ,"snpName"] %in% all_SNPpop_freqs_T[,"ID"],]

# see where the ref not the same
probList <- unlist(lapply(1:length(all_SNPpop_freqs_T[,1]), function(x) sum(as.vector(allDB_chr_T2[allDB_chr_T2[,"snpName"] == as.vector(all_SNPpop_freqs_T[x,"ID"]), "refNT"]) != as.vector(all_SNPpop_freqs_T[x,"REF"]))))

# write to output how many:
write(paste("*** There were ", length(all_SNPpop_freqs_T[which(probList != 0),1]), " SNPs removed because 1000Genomes ref != DB ref ***", sep=""), stdout())

# remove the problem ones:
all_SNPpop_freqs <- all_SNPpop_freqs_T[which(probList == 0),]
rownames(all_SNPpop_freqs)<- NULL

# and now remove those from the allDB as well
allDB_chr <- allDB_chr_T2[allDB_chr_T2[ ,"snpName"] %in% all_SNPpop_freqs[,"ID"],]
rownames(allDB_chr)<- NULL


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write("*** Computing DAF and TAFs of each interaction ***", stdout())

# Split it up to run: 
n <- 10000
DBofDBs <- split(allDB_chr, rep(1:ceiling(nrow(allDB_chr)/n), each=n, length.out=nrow(allDB_chr)))

partNum <- 0
for(db in DBofDBs){
    partNum <- partNum + 1
    allDB3 <- db

    ### Make a db with my info of interest:
    # Start with the cases where the target is the ref:
    allDB3_Tref <- allDB3[allDB3[ , "targ_is_ref"] == '1', ]

    # Only do this if there are this type:
    if (length(allDB3_Tref[,1]) > 0){
        # add columns for TAF and DAF of each population
        for (group in c("EAS", "AMR", "AFR", "EUR", "SAS")){
            allDB3_Tref[, paste(group, "_TAF", sep = "")] <- NA
            allDB3_Tref[, paste(group, "_DAF", sep = "")] <- NA
        }   
        
        # For each population, fill in the TAF for each SNP  ** The 1000genomes AFs are of the NON-REF ALLELE **
        for (snp in unique(allDB3_Tref[ ,"snpName"])){
            allDB3_Tref[allDB3_Tref[ ,"snpName"] == as.vector(snp), "EAS_TAF"] <- 1 - as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] ==  as.vector(snp), "EAS_AF"]))
            allDB3_Tref[allDB3_Tref[ ,"snpName"] == as.vector(snp), "AMR_TAF"] <- 1 - as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "AMR_AF"]))
            allDB3_Tref[allDB3_Tref[ ,"snpName"] == as.vector(snp), "AFR_TAF"] <- 1 - as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "AFR_AF"]))
            allDB3_Tref[allDB3_Tref[ ,"snpName"] == as.vector(snp), "EUR_TAF"] <- 1 - as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "EUR_AF"]))
            allDB3_Tref[allDB3_Tref[ ,"snpName"] == as.vector(snp), "SAS_TAF"] <- 1 - as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "SAS_AF"]))
        }
       # For each population, fill in the DAF line by line 
        for( i in 1:length(allDB3_Tref[,1])){
             if(as.vector(allDB3_Tref[i,"targNT"]) != as.vector(allDB3_Tref[i,"ancestralAllele"])){
                allDB3_Tref[i, "EAS_DAF"] <- as.numeric(as.vector(allDB3_Tref[i, "EAS_TAF"]))
                allDB3_Tref[i, "AMR_DAF"] <- as.numeric(as.vector(allDB3_Tref[i, "AMR_TAF"]))
                allDB3_Tref[i, "AFR_DAF"] <- as.numeric(as.vector(allDB3_Tref[i, "AFR_TAF"])) 
                allDB3_Tref[i, "EUR_DAF"] <- as.numeric(as.vector(allDB3_Tref[i, "EUR_TAF"]))
                allDB3_Tref[i, "SAS_DAF"] <- as.numeric(as.vector(allDB3_Tref[i, "SAS_TAF"]))       
             }else{
                allDB3_Tref[i, "EAS_DAF"] <- 1 - as.numeric(as.vector(allDB3_Tref[i, "EAS_TAF"]))
                allDB3_Tref[i, "AMR_DAF"] <- 1 - as.numeric(as.vector(allDB3_Tref[i, "AMR_TAF"]))
                allDB3_Tref[i, "AFR_DAF"] <- 1 - as.numeric(as.vector(allDB3_Tref[i, "AFR_TAF"])) 
                allDB3_Tref[i, "EUR_DAF"] <- 1 - as.numeric(as.vector(allDB3_Tref[i, "EUR_TAF"]))
                allDB3_Tref[i, "SAS_DAF"] <- 1 - as.numeric(as.vector(allDB3_Tref[i, "SAS_TAF"]))  
                }
        }
    }

    # and now where the target is derived allele
    allDB3_NTref <- allDB3[allDB3[ , "targ_is_ref"] == '0', ]

    # Only do this if there are this type:
    if (length(allDB3_NTref[,1]) > 0){
        # add columns for TAF and DAF of each population
        for (group in c("EAS", "AMR", "AFR", "EUR", "SAS")){
            allDB3_NTref[, paste(group, "_TAF", sep = "")] <- NA
            allDB3_NTref[, paste(group, "_DAF", sep = "")] <- NA
        }
        # For each population, fill in the TAF for each SNP   ** The 1000genomes AFs are of the NON-REF ALLELE **
        for (snp in unique(allDB3_NTref[ ,"snpName"])){
            allDB3_NTref[allDB3_NTref[ ,"snpName"] == as.vector(snp), "EAS_TAF"] <- as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "EAS_AF"]))
            allDB3_NTref[allDB3_NTref[ ,"snpName"] == as.vector(snp), "AMR_TAF"] <- as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "AMR_AF"]))
            allDB3_NTref[allDB3_NTref[ ,"snpName"] == as.vector(snp), "AFR_TAF"] <- as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "AFR_AF"]))
            allDB3_NTref[allDB3_NTref[ ,"snpName"] == as.vector(snp), "EUR_TAF"] <- as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "EUR_AF"]))
            allDB3_NTref[allDB3_NTref[ ,"snpName"] == as.vector(snp), "SAS_TAF"] <- as.numeric(as.vector(all_SNPpop_freqs[all_SNPpop_freqs[ ,"ID"] == as.vector(snp), "SAS_AF"]))
        }
       # For each population, fill in the DAF line by line  
        for( i in 1:length(allDB3_NTref[,1])){ 
             if(as.vector(allDB3_NTref[i,"targNT"]) != as.vector(allDB3_NTref[i,"ancestralAllele"])){
                allDB3_NTref[i, "EAS_DAF"] <- as.numeric(as.vector(allDB3_NTref[i, "EAS_TAF"]))
                allDB3_NTref[i, "AMR_DAF"] <- as.numeric(as.vector(allDB3_NTref[i, "AMR_TAF"]))
                allDB3_NTref[i, "AFR_DAF"] <- as.numeric(as.vector(allDB3_NTref[i, "AFR_TAF"]))
                allDB3_NTref[i, "EUR_DAF"] <- as.numeric(as.vector(allDB3_NTref[i, "EUR_TAF"]))
                allDB3_NTref[i, "SAS_DAF"] <- as.numeric(as.vector(allDB3_NTref[i, "SAS_TAF"]))
             }else{
                allDB3_NTref[i, "EAS_DAF"] <- 1 - as.numeric(as.vector(allDB3_NTref[i, "EAS_TAF"]))
                allDB3_NTref[i, "AMR_DAF"] <- 1 - as.numeric(as.vector(allDB3_NTref[i, "AMR_TAF"]))
                allDB3_NTref[i, "AFR_DAF"] <- 1 - as.numeric(as.vector(allDB3_NTref[i, "AFR_TAF"]))
                allDB3_NTref[i, "EUR_DAF"] <- 1 - as.numeric(as.vector(allDB3_NTref[i, "EUR_TAF"]))
                allDB3_NTref[i, "SAS_DAF"] <- 1 - as.numeric(as.vector(allDB3_NTref[i, "SAS_TAF"]))        
                }
        }
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    
    # add them together:
    allDB_fin <- rbind(allDB3_Tref, allDB3_NTref)
    
    # Name the file part to save
    fileSave <- paste("PopData_2ndRound/chr_", chrNum, "_p_", partNum, ".txt", sep = "")
    fileSaveH <- paste("PopData_2ndRound/chr_", chrNum, "_Header_", partNum, ".txt", sep = "")
#      ### FOR REV.COMP
#     fileSave <- paste("PopData_2ndRound/RevComp_chr_", chrNum, "_p_", partNum, ".txt", sep = "")
#     fileSaveH <- paste("PopData_2ndRound/RevComp_chr_", chrNum, "_Header_", partNum, ".txt", sep = "")
#    

    # For the output, keep the columnnames only for the first one (for the overll header when cat together - and for x and y)
       # added format, because for some strange reason it was adding a bunch of digdgits to some, eventhough they were not in allDB_fin
    if ((partNum == "1") & (chrNum == "1")){
        write.table(format(allDB_fin, digits = 4), file = fileSaveH, 
                    sep = "\t", row.names = F, quote = F, col.names = T)
        }else{
            write.table(format(allDB_fin, digits = 4), file = fileSave, 
                    sep = "\t", row.names = FALSE, quote = F, col.names = F)
        }
      
}


write("*** DONE, check output files (each chr is split into sepctions of 10000 for speed) ***", stdout())

