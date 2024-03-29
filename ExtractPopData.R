#!/usr/bin/Rscript


# Script for taking in all the 1000 genomes snp info (one chromosome at a time) 
# (chr Y has a special part because only the male individuals are included).

# usage: R script ./ExtractPopData.R 21

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Take in the chromosome number:
args<-commandArgs(T)
chrNum <- args[1] 

# If doing normal or rev. comp change the snp db that is loaded below! AND CHANGE OUTPUT!
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

write("*** Starting import of the SNP data for the chromosome  ***", stdout())

# Import the SNP info for that chromosome:
if (chrNum == "X"){
        pathChrSNPs <- "RefComparisonWork/SNPsIndividualsALL/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf"
    }else if(chrNum == "Y"){
        pathChrSNPs <- "RefComparisonWork/SNPsIndividualsALL/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf"
    }else{
        pathChrSNPs <- paste("RefComparisonWork/SNPsIndividualsALL/ALL.chr", chrNum ,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf", sep = "")
}

chrSNPs <- read.table(pathChrSNPs, sep = "\t")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Take only the first 8 columns (takes pop info without genotype of each individual)
# Also select only those that have a snp ID, so remove those that are "."
chrSNPsTEMP <- chrSNPs[which(chrSNPs[,3]!= "."),1:8]
colnames(chrSNPsTEMP) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

# so not so much stored:
chrSNPs <- 0

# To speed up the process:
################################
# # Load the list of all DB snps (called DBsnps)
load("snpNames_from_CompleteDB.Rdata")

# FOR THE REV. COMP VERSION!! (change back to other if doing normal)
#load("snpNames_from_CompleteDB_RevComp.Rdata")
################################

# Take only those snps in my overall DB:
snpsChrPop <- chrSNPsTEMP[which(chrSNPsTEMP[,"ID"] %in% DBsnps), ]
rownames(snpsChrPop) <- NULL

# Split the 8th column into seperate columns
newCols <- lapply(as.vector(snpsChrPop[,8]), function(x) strsplit(x, split = ";")[[1]])

# to be sure the name is correct, assign each to what is is the list
for(i in 1:(length(newCols[[1]]) - 1)){
   assign(strsplit(newCols[[1]], split = "=")[[i]][1], as.numeric(unlist(lapply(newCols, function(x) strsplit(x, split = "=")[[i]][2]))))
}

# Re-combine into one db:
all_SNPpop_freqs <- cbind(snpsChrPop[ ,c("CHROM", "POS", "ID", "REF", "ALT")], AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF, EAS_AF)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Create the output file:
SNP_output <- paste("PopData_byChr/PopData_Chr_", chrNum , ".txt", sep = "")

write.table(all_SNPpop_freqs, file = SNP_output, sep = "\t", row.names = F, quote = F, col.names = T)

# if (chrNum == "X" | chrNum == "Y" | chrNum == "1"){
#     write.table(all_SNPpop_freqs, file = SNP_output , 
#             sep = "\t", row.names = F, quote = F, col.names = T)
#     }else{
#     write.table(all_SNPpop_freqs, file = SNP_output , 
#             sep = "\t", row.names = F, quote = F, col.names = F)
#         }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

write("*** DONE, check optupt ***", stdout())




