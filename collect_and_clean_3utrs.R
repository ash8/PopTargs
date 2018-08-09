#!/usr/bin/Rscript

##
## This script takes 3'UTRs from BiomaRt and filters them. 
## 
## Filters: only utrs with sequence avaialble, and for protein coding, and if it 
## is negative strand, takes the reverse complement, so all are positive strand.
## Output: fasta format sequence with sequence information in the header


#############################
# To run this in a bash script, include the chromosome number to use (ex chr 1)
# usage: Rscript ./collect_and_clean_3utrs.R 1
chrNum<-commandArgs(T)

# IF not running as bash script, be sure to set chrNum to the chromosome you want 

#############################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Build a biomaRt query:
# Just like in web version- there are filters and attributes.

# from the manual: "The getBM function is the main query function in biomaRt. 
# It has four main arguments:
	# 1. attributes: is a vector of attributes that one wants to retrieve 
	#    (= the output of the query).
	# 2. filters: is a vector of lters that one wil use as input to the query.
	# 3. values: a vector of values for the lters. In case multple lters are in 
	#    use, the values argument requires a list of values where each position 
	#    in the list corresponds to the position of the lters in the lters 
	#    argument (see examples below).
	# 4. mart: is an object of class Mart, which is created by the useMart function.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get the 3'UTR sequences:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library("biomaRt")

# choose the dataset and BioMart database:  
# (to see the options for each use: listMarts() and listDatasets(ensembl))
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

# to see the filters available: 
# listFilters(ensembl)

# to see attributes available:
# listAttributes(ensembl)


# set the attributes and filters:
utr_attributes = c("3utr", "ensembl_gene_id", "ensembl_transcript_id", 
					"chromosome_name", "3_utr_start", "3_utr_end", 
					"transcript_biotype", "strand")
					
utr_filters = c("chromosome_name")

write("*** RETRIEVING 3UTR DATA ***", stdout())

# Retrieve the data:  ### here is where you choose the chromosome ###
utrs_db = getBM(attributes = utr_attributes, filters = utr_filters, 
				values = chrNum, mart = ensembl)
			
write("*** REMOVING UNWANTED RECORDS AND CREATING FASTA FILE ***", stdout())

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Put the database into fasta format for use with Seed Vicious
# also remove irrelevant accessions 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Remove seqs with no sequence available
remove_utrs<-c()
keep_utrs<-c()
for (i in 1:length(utrs_db[,1])) {
	
	if(utrs_db[i,1] == "Sequence unavailable") {
		remove_utrs<-rbind(remove_utrs, utrs_db[i,])
		
	}else{
		keep_utrs<-rbind(keep_utrs, utrs_db[i,])
	}
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Remove seqs that are not protein coding
keep_utrs2<-c()
for (i in 1:length(keep_utrs[,1])) {
	
	if(keep_utrs[i,"transcript_biotype"] == "protein_coding") {
		keep_utrs2<-rbind(keep_utrs2, keep_utrs[i,])
		
	}else{
		remove_utrs<-rbind(remove_utrs, keep_utrs[i,])
	}
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Put it into fasta format and take the reverse complement of '- strand' versions

library("Biostrings")

utrs_set<-DNAStringSet()
utrs_set_names<-c()
for (i in 1:length(keep_utrs2[,1])) {
	
	newSeqName<-paste(keep_utrs2[i,c("ensembl_gene_id","ensembl_transcript_id",
	                  "chromosome_name", "3_utr_start", "3_utr_end")],
	                  collapse = "|")
	# or keep_utrs2[1,2:6] if I dont want to call them specifically.
	
	if(keep_utrs2[i,"strand"] == "1") {
		utrs_set<-c(utrs_set, DNAStringSet(keep_utrs2[i,"3utr"]))
		utrs_set_names<-c(utrs_set_names, newSeqName)
		
	}else{
		utrs_set<-c(utrs_set, 
					reverseComplement(DNAStringSet(keep_utrs2[i,"3utr"])))
		utrs_set_names<-c(utrs_set_names, newSeqName)
	}
} 
names(utrs_set)<-utrs_set_names


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Output the file for use in SeedVicious


outFilename <- paste("./chr", chrNum ,"/biomart_3utrs_chr", chrNum , ".txt", sep = "")
writeXStringSet(utrs_set, file = outFilename)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 





