# PopTargs #
These are the scripts used to create the MySQL database that is used by PopTargs.essex.ac.uk.
The pipeline can be altered to create similar databases with different species, it may need to be adjusted to fit your file names. 

## The general pipeline is as follows: ##

1)	Get the 3’ UTRs with: "collect_and_clean_3utrs.R chr#"
2)  Download the mature human miRNA from miRBase.org
3)  Download the SNPs with: "download_SNPs_BiomaRt.R chr#"
4)	Run SeedVicious 
5)	Find the SNP/Target interactions: "Find_SNPs_in_Targets.R chr#"
6)	Put each chromosome into the MySQL database "setup_mySQL_DB_new.R chr#"
7)	Combine the chromosomes into one big table for the database 
8)	Download all the populations data from 1,000 Genomes: ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/
9)	Extract the populations data: "ExtractPopData.R chr#"
10)	Take the relevant populations data and set it up for the database "PopData2ndRound.R chr#"
11)	Combine the chromosomes into one big table for the database
12)	Upload to MySQL
