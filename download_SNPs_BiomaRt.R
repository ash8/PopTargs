# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# # Script for building a biomaRt query the SNPs:

# usage (ex chr 21): download_SNPs_BiomaRt.R 21

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# To run this in a bash script, include the path to the seedV output and 
# chromosome number to use (ex chr 1)
# usage: Rscript ./Find_SNPs_in_Targets.R seedV_path.txt 1
args<-commandArgs(T)

# If not running as bash script, be sure to set the following manually:
chrNum <- args[1] 

## Because most of the chromosomes are too big to import all of their SNPs at 
## once, I will do them in smaller portions. And as creating empty DBs doesnt seem
## to cause a problem, I will just do it for all:

library("biomaRt")

write("*** DOWNLOADING SNPS (this may take a while) ***", stdout())

# # choose the dataset and BioMart database:  
# # (to see the options for each use: listMarts() and listDatasets(snpmart))
snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# to see the filters available: listFilters(snpmart)
# to see attributes available: listAttributes(snpmart)

snp_attributes <- c("refsnp_id", "refsnp_source", "chr_name", "chrom_start", "chrom_end", "chrom_strand", "allele", "minor_allele", "minor_allele_freq", "allele_1")
					
#### DESCRIPTION OF 'snp_attributes' ####
# # refsnp_id			Variant Name
# # refsnp_source		Variant source
# # chr_name			Chromosome name
# # chrom_start			Chromosome position start (bp)
# # chrom_end			Chromosome position end (bp)
# # chrom_strand		Strand
# # allele				Variant Alleles
# # minor_allele		Minor allele (ALL)
# # minor_allele_freq 	1000 Genomes global Minor Allele Frequency (all individuals)
# # allele_1 			Ancestral allele
		
snp_filters <- c("chromosomal_region")

regionPart1 <- paste(chrNum, ":1:50000000", sep = "")
regionPart2 <- paste(chrNum, ":50000000:100000000", sep = "")
regionPart3 <- paste(chrNum, ":100000000:150000000", sep = "")
regionPart4 <- paste(chrNum, ":150000000:200000000", sep = "")
regionPart5 <- paste(chrNum, ":200000000:", sep = "")

snp_db1 <- getBM(attributes = snp_attributes, filters= snp_filters, values = regionPart1, 
				mart = snpmart)
snp_db2 <- getBM(attributes = snp_attributes, filters= snp_filters, values = regionPart2, 
				mart = snpmart)
snp_db3 <- getBM(attributes = snp_attributes, filters= snp_filters, values = regionPart3, 
				mart = snpmart)
snp_db4 <- getBM(attributes = snp_attributes, filters= snp_filters, values = regionPart4, 
				mart = snpmart)
snp_db5 <- getBM(attributes = snp_attributes, filters= snp_filters, values = regionPart5, 
				mart = snpmart)
								
snp_db <- rbind(snp_db1, snp_db2, snp_db3, snp_db4, snp_db5)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FOR BASH SCRIPT OUTPUT:

# create the file name based on chromosome number
SNPs_output <- paste("./chr", chrNum ,"/SNPdb_BiomaRt_chr", chrNum , ".txt", sep = "")

# SAVE the database:
write.table(snp_db, file = SNPs_output, 
			sep = "\t", row.names = F)



