
LORCAN "LOng Read Consenus ANalysis" 	 
====================================
LORCAN v1.8

# DESCRIPTION:
Using nanopore reads in fastq format, the pipeline produces consensus sequences     
for each barcoded sample, based on best match to a reference database.    

# PARAMETERS:
*Input dir 		
   The path to the directory containing the fastq files (not demultiplexed)    
		
*Output dir path
   Path to the output directory. The "output" directory is automatically created and appended to the path.    
		
*Number of cpus
The max value is the value available on the current system.		
		
*Minimum number of reads per barcoded sample (threshold)
The number of reads per taxonomic group under which the consensus per taxonomic group 		
will not be calculated.		
		
*Databases
A database flat file (fasta) needs to be obtained and formatted to be compatible with BLASTN     

This is an example of formating for SSU-rDNA-mk37/16S_TS-stringent_dedup.fasta (leBiBi) database:   
`module load Blast/blast/2.6.0` #system specific     
`makeblastdb -in  /path2dir/SSU-rDNA-mk37/16S_TS-stringent_dedup.fasta -dbtype nucl`     
The database path should be indicated in the lib/config.pm file.     

A corresponding .taxdict.csv file must be prepared using the same database file (fasta).  

For each entry of the database file,the .taxdict.csv file contains the following fields:     
[Accession number],[taxonomic classication],[full fasta header information]    
e.g. URS00000B1AF5,Abiotrophia_defectiva,Abiotrophia_defectiva~v~TT~URS00000B1AF5=Bacteria-Fir...
			
A utility script can be found in lib/CreateTaxoDict_LeBiBI16S.pl to prepare such file using the LeBiBiref database.
		
Other databases can be added by modifying the config.pm file appropriately.   
		
*Email address
Under normal mode, an email is sent to the provided email address when the job is complete.    
When verbose mode is on, no email is sent.   

*Verbose mode
When on, more detailed output is provided. When off, several intermediary files are removed.    
		
		
# CONFIGURATION:
**config.pm**   
This file contains the absolute paths to the dependencies (executable, databases).    
	Each field may need to be changed depending on your own configuration and system.    

# Example:
`
FastqDir=/storage/20190213_x3_Amplikon/fastq_pass/   
cwd=/Analyses/Pipelines/dev/LORCAN/v1.8   
/pathto/main.pl -V -i $FastqDir -o $cwd/20190213_x3_Amplikon_2 -L log_main.txt -I $cwd/20190213_x3_Amplikon_2/sample_id.txt -n 20 -m 10 -M 3000 -P 100 -D 5  -d BiBi16SLong
`     

		
