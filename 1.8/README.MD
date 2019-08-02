
LORCAN "LOng Read Consensus ANalysis" 	 
====================================
LORCAN version: **1.8**

## SHORT DESCRIPTION   
From basecalled nanopore reads in FASTQ format, the pipeline produces consensus sequences     
for each barcoded sample, based on read alignment to best matches from a customized reference database.  
All analysis steps are automated and a final report (PDF) is produced for each barcoded sample.      

## INSTALLATION       
***Disclaimer***: The pipeline was developed and tested with Scientific Linux 7.4 (x86-64 server architecture),    
and may only work in a GNU/Linux environment. Other OS have not been tested and are not likely to be compatible with the software.    

The recommended way is to use a *conda* (*anaconda*, *miniconda*) environment to install the software and its dependencies.     

Here is a step-by-step installation guide using a *conda* environment:     
1.	If not yet done, install *conda*:     
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html     
To test your installation, in your terminal, type `conda list`     
If you read “bash: conda: command not found...”, then the conda installation was not successful and LORCAN cannot be installed with conda.     
Otherwise, you should see some text in the console, starting with "#packages in environment at...".       

2.	Create a dedicated environment (to do once):     
```
conda create -n LORCAN-env         
```

3.	Install the dependencies in the activated environment:      
```
source activate LORCAN-env         
conda install seqkit=0.8.0; seqkit version #OK       
conda install porechop=0.2.3_seqan2.1.1
conda install mafft=v7.313
conda install Gblocks           		
conda install iqtree=1.6.11	             	
conda install minimap2=2.17	       
conda install samtools       
conda install blast       
cpanm Parallel::ForkManager 	       
cpanm FindBin
conda install -c anaconda gcc_linux-64
cpanm PDF::API2       
```

To save some space, optionally remove unneeded files:    
```
conda clean --tarballs; conda clean --package        
```


4.	Clone the *LORCAN* binaries from GitHub to your own machine:          
(move first to the directory where you want to install LORCAN)      
```
git clone https://github.com/aramette/LORCAN.git    
```

5.	Test whether the pipeline is correctly installed:    
```
cd LORCAN/1.8/lib
perl lorcan.pl -h
```

6. Optional, but recommended:      
To be able to run the command everywhere on your system with the current username, edit your `~./bashrc` file by running the following command that adds an alias to the bashrc file:
```
printf "alias LORCAN=\'perl /path_to_your_directory/LORCAN/1.8/lib/lorcan.pl\'">> ~/.bashrc

```
Then run:
```
source ~./bashrc
```
to take the modification into account (next time your start a session, the changes will already be written in your .bashrc file).      
now you can call the program directly using the alias *LORCAN* on the command line:         
```
LORCAN -h
```

## CONFIGURATIONS





```
```
5.	Test *LORCAN*            
```
```

What to modify before starting:       



To exit from the LORCAN-env, type: 
```
source deactivate LORCAN-env        
```


# PARAMETERS        

**Arguments:**      
**-h** help information            
**-V** verbose mode (for debugging). When selected, more detailed output is provided. When off, several intermediary files are removed.          
**-v** version number      
**-d** choice of the database for alignment (e.g. BiBi16S, ADV...) Must match declaration in config.pm. Must be primer customized (no primer)       
**-D** choice of delta value around the modal sequence length (e.g. 5 means modal sequence length +-5 bp).       
**-E** email address: e.g.  your.name\@your.institute.com. Under normal mode, an email is sent to the provided email address when the job is complete. When verbose mode is on, no email is sent.   
**-i** directory containing the basecalled FASTQ files (full path) not yet demultiplexed.          
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
**-o** full path to the output directory e.g. /path_to/output. The "output" directory folder is automatically created and appended to the indicated path.            
**-P** minimum  number of read aligned for a reference for the latter to be further considered (e.g 100 for 100 reads). In other words, this indicates the number of reads per taxonomic group under which the consensus per taxonomic group 		
will not be calculated.		
        

## Example   
```
FastqDir=/storage/20190213_x3_Amplikon/fastq_pass/
cwd=/Analyses/Results   
perl lorcan.pl -V -i $FastqDir -o $cwd/myOutput -L log_main.txt -I $cwd/20190213_x3_Amplikon_2/sample_id.txt -n 20 -m 10 -M 3000 -P 100 -D 5  -d BiBi16SLong`
```     
This command will produce an output structure as follows:     
The results will be found in $cwd/myOutput/...        
XXXXXXXXXX             

# How to cite?       		
xxxx.    


