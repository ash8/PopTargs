# PopTargs
# These are the scripts used to create the MySQL database that is used by PopTargs.essex.ac.uk.
# The pipeline is a long process and scripts may need to be adjusted to fit your file names. 

# The general pipeline is as follows:

1)	Get the 3’ UTRs with: “collect_and_clean_3utrs.R”
2)  Download the mature human miRNA from miRBase.org. Split the fasta file into many small files for parallel SeedVicious runs.
3)  Download teh SNPs with: download_SNPs_BiomaRt.R 
4)	Run SeedVicious 
5)	Find the SNP/Target interactions: Find_SNPs_in_Targets.R 
6)	Put each chromosome into the MySQL database ./setup_mySQL_DB_new.R 
7)	Combine the chromosomes into one big table for the database 
8)	Download all the populations data
9)	Extract the populations data. (Using ExtractPopData_array.sh to run ExtractPopData.R
10)	Take the relevant populations data and set it up for the database. (Use PopData2ndRound_array.sh to run PopData2ndRound.R) 	
11)	Cat populations together
12)	Upload to MySQL
