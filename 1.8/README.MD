
LORCAN "LOng Read Consenus ANalysis" 	 
====================================
LORCAN v1.8

# DESCRIPTION   
Using nanopore reads in fastq format, the pipeline produces consensus sequences..     
for each barcoded sample, based on best match to a reference database.    

# PARAMETERS  

Arguments:
	-h help information....
	-V verbose mode (for debugging)
	-v version number
	-d choice of the database for alignment (e.g. BiBi16S, ADV...) Must match declaration in config.pm. Must be primer customized (no primer) 
	-D choice of delta value around the modal sequence length (e.g. 5 means modal sequence length +-5 bp). 
	-E Email address: e.g.  your.name\@your.institute.com 
	-i directory containing the fastq files (full path)
	-I the full path to the sample id file (what is after # is a comment not read; RunName:/Date:/BCnber: no space)
	-L log file name e.g. log.txt (not the path)  
	-m minimum number of reads per barcode to further proceed  
	-M maximum number of reads to retain per sample (e.g. 3000)
	-n Nber of threads e.g. 20 
	-o Full path to the output directory e.g. /path_to/output  
	-P minimum  number of read aligned for a reference for the latter to be further considered (e.g 100 for 100 reads)  
Usage example:
	perl lib/main.pl -V -i \$testDataDir -n 20 -m 10 -o $cwd/output_test_main_16S_delta5 -I $cwd/sample_id.txt -L log_main.txt -P 30 -D 5 -d BiBi16S
	
**Input dir**     		
The path to the directory containing the fastq files (not demultiplexed).          
		
**Output dir path**   
Path to the output directory. The "output" directory is automatically created and appended to the path.    
		
**Number of cpus**   
The max value is the value available on the current system.		
		
**Minimum number of reads per barcoded sample (threshold)**   
The number of reads per taxonomic group under which the consensus per taxonomic group 		
will not be calculated.		
		
**Databases**   
A database flat file (fasta) needs to be obtained and formatted to be compatible with BLASTN, as follows.     

`module load Blast/blast/2.6.0` #system specific     

`makeblastdb -in  /path2dir/SSU-rDNA-mk37/16S_TS-stringent_dedup.fasta -dbtype nucl`     

The database path should be indicated in the lib/config.pm file.     

A corresponding .taxdict.csv file must be prepared using the same database file (fasta). See below.      

Other databases can be added by modifying the config.pm file appropriately.    

**taxdict.csv file**   
For each entry of the database file, the taxdict.csv file contains the following fields:     
[Accession number],[taxonomic classication],[full fasta header information]      
e.g. URS00000B1AF5,Abiotrophia_defectiva,Abiotrophia_defectiva~v~TT~URS00000B1AF5=Bacteria-Fir      
			
A utility script can be found in lib/CreateTaxoDict_LeBiBI16S.pl to prepare such file using the LeBiBiref database.     
		
      
		
**Email address**   
Under normal mode, an email is sent to the provided email address when the job is complete.    
When verbose mode is on, no email is sent.   

**Verbose mode**
When selected, more detailed output is provided. When off, several intermediary files are removed.    
		
		
# CONFIGURATION   
**config.pm**   
This file contains the absolute paths to the dependencies (executable, databases).    
Each field may need to be changed depending on your own configuration and system.    

# Example

`FastqDir=/storage/20190213_x3_Amplikon/fastq_pass/`   

`cwd=/Analyses/Pipelines/dev/LORCAN/v1.8`   

`/pathto/main.pl -V -i $FastqDir -o $cwd/20190213_x3_Amplikon_2 -L log_main.txt -I $cwd/20190213_x3_Amplikon_2/sample_id.txt -n 20 -m 10 -M 3000 -P 100 -D 5  -d BiBi16SLong`
     

		
