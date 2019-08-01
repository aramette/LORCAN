
LORCAN "LOng Read Consenus ANalysis" 	 
====================================
LORCAN version: 1.8

## SHORT DESCRIPTION   
From basecalled nanopore reads in FASTQ format, the pipeline produces consensus sequences     
for each barcoded sample, based on read alignment to best matches from a customized reference database.  
All analysis steps are automated and a final report (PDF) is produced for each barcoded sample.      

## INSTALLATION INSTRUCTIONS      
**Disclaimer**: The pipeline was developed and tested with Scientific Linux 7.4 (x86-64 server architecture),    
and may only work in a GNU/Linux environment. Other OS have not been tested and are not likely to be compatible with the software.    

The recommended way is to use a conda (anaconda, miniconda) environment to install the software and its dependencies.     

Here is a step-by-step installation guide using a conda environment:     
1.	If not yet done, install conda     
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html     
To test your installation, in your terminal window or Anaconda Prompt, run `conda list`     
If you read “bash: conda: command not found...”, then the conda installation was not successful and LORCAN cannot be installed with conda.     
Otherwise, you should see an output starting with "#packages in environment at...".       

2.	Create environment     

module load conda/3 #this depends on your system     
conda create -n LORCAN-env      
source activate LORCAN-env      
3.	Install the dependencies (**NOTETOS: cf LORCAN main script for Perl scripts, and config.pm file too      
**NOTETOS: (Give precise version numbers too for the installation)       

[OK] conda install seqkit=0.8.0; seqkit version #OK       
[OK] conda install porechop; porechop --version #OK       
[OK] conda install mafft ; mafft --version #OK       			
[OK] conda install Gblocks ;   # OK needed ?       		
[OK] conda install iqtree=1.6.11	      #OK       	
       [OK] conda install minimap2=2.17	       
       [OK] conda install samtools       
       [OK] conda install BLAST       
       [OK] cpanm Parallel::ForkManager #OK	       
       [OK] cpanm FindBin       
       [OK] cpanm PDF::API2       



source deactivate LORCAN-env        
conda clean --tarballs; conda clean --package #to save some space       

4.	Pull the git       
5.	Test -h       
6.	Close env       

What to modify before starting:       

Taxdict file        
Db (custom if needed, blast)       

Config.pm       
(create s script to find the dependencies to be put in the config.file “which blast”)       
our $PERL			="/software/KM/perl5/perlbrew/perls/perl-5.28.0/bin/perl";       
our $PORECHOP		="/software_conda/3/bin/porechop";       
our $samtools		="module load UHTS/samtools/1.4;samtools";       		
our $BLASTN			="module load Blast/blast/2.6.0;blastn";       
our $Rscript 		="/software/KM/R-3.5.0_KM_WGS/bin/Rscript";       
our $mafft 			="/software_conda/3/bin/mafft";       
our $Gblocks		="/software_conda/3/bin/Gblocks";       
our $iqtree			="/software_conda/3/bin/iqtree";       
also path to db and blast       


# PARAMETERS        

**Arguments:**      
**-h** help information            
**-V** verbose mode (for debugging)      
**-v** version number      
**-d** choice of the database for alignment (e.g. BiBi16S, ADV...) Must match declaration in config.pm. Must be primer customized (no primer)       
**-D** choice of delta value around the modal sequence length (e.g. 5 means modal sequence length +-5 bp).       
**-E** email address: e.g.  your.name\@your.institute.com       
**-i** directory containing the basecalled FASTQ files (full path)      
**-I** the full path to the sample id file. The sample id file is made of the first two lines with RunName and RunDate, followed by      
	the sample names associated with each barcode, separated by colon (:). The whole BC line is copied to the final report. 
	This is an example of file with 3 barcoded samples:     
	`
	RunName:20190220_x4_Amplikon      
	RunDate:21.02.2019      	
	#      
	BC17:2019-067      	
	BC18:2019-068      	
	BC19:2019-069      
	`
**-L** log file name e.g. log.txt (not the path)        
**-m** minimum number of reads per barcode to further proceed        
**-M** maximum number of reads to retain per sample (e.g. 3000)      
**-n** nber of threads e.g. 20       
**-o** full path to the output directory e.g. /path_to/output        
**-P** minimum  number of read aligned for a reference for the latter to be further considered (e.g 100 for 100 reads)        

# Usage example   
`perl lib/main.pl -V -i \$testDataDir -n 20 -m 10 -o $cwd/output_test_main_16S_delta5 -I $cwd/sample_id.txt -L log_main.txt -P 30 -D 5 -d BiBi16S`      
	
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
     

# How to cite?       		
xxxx.    


