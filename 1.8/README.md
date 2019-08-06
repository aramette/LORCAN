*LORCAN* "LOng Read Consensus ANalysis" 	 
====================================
version: **1.8**

## SHORT DESCRIPTION   
From basecalled nanopore reads in FASTQ format, the pipeline produces consensus sequences for each barcoded sample, based on alignment of the reads to a customized reference sequence database. All analysis steps are automated and a final report (PDF) is produced for each barcoded sample.      


## Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [USAGE](#USAGE)
- [ARGUMENTS](#ARGUMENTS)
- [INSTALLATION](#INSTALLATION)
- [CONFIGURATION](#CONFIGURATION)
- [TEST](#TEST)
- [FEEDBACK](#FEEDBACK)
- [CITATION](#CITATION)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## USAGE
```
FastqDir=/storage/20190213_x3_Amplikon/fastq_pass/
cwd=/Analyses/Results   

perl lorcan.pl -V -i $FastqDir -o $cwd/myOutput -L log_main.txt -I $cwd/20190213_x3_Amplikon_2/sample_id.txt -n 20 -m 10 -M 3000 -P 100 -D 5  -d MyDataBase
```     

## ARGUMENTS                 
<P ALIGN=RIGHT><a href="#top">↥ back to top</a></P>          

| Flag | Description | 
| :--------------: |:-----------| 
|**-h**| help information|
|**-V** |verbose mode (for debugging). When selected, more detailed output is provided. When off, several intermediary files are removed.|          
|**-v**| version number |     
|**-d**| choice of the database for alignment (e.g. BiBi16S, ADV...) Must match declaration in config.pm. Must be primer customized (no primer) |      
|**-D**| choice of delta value around the modal sequence length (e.g. 5 means modal sequence length +-5 bp).|       
|**-E**| email address: e.g.  your.name\@your.institute.com. Under normal mode, an email is sent to the provided email address when the job is complete. When verbose mode is on, no email is sent.|   
|**-i**| directory containing the basecalled FASTQ files (full path) not yet demultiplexed.|          
|**-I**| the full path to the sample id file. The sample id file is made of the first two lines with RunName and RunDate, followed by the sample names associated with each barcode, separated by colon (:). The whole BC line is copied to the final report.|
|**-L**| log file name e.g. log.txt (not the path) |       
|**-m**| minimum number of reads per barcode to further proceed |       
|**-M**| maximum number of reads to retain per sample (e.g. 3000)|      
|**-n**| nber of threads e.g. 20|       
|**-o**| full path to the output directory e.g. /path_to/output. The "output" directory folder is automatically created and appended to the indicated path.|            
|**-P**| minimum  number of read aligned for a reference for the latter to be further considered (e.g 100 for 100 reads). In other words, this indicates the number of reads per taxonomic group under which the consensus per taxonomic group will not be calculated.|		


This is an example of *sample_id.txt* file with 3 barcoded samples:     
```
RunName:20190220_x4_Amplikon      
RunDate:21.02.2019      	
#      
BC17:2019-067      	
BC18:2019-068      	
BC19:2019-069      
```

## INSTALLATION             
<P ALIGN=RIGHT><a href="#top">↥ back to top</a></P>          

:warning: ***Disclaimer***: The pipeline was developed and tested with Scientific Linux 7.4 (x86-64 server architecture), and may only work in a GNU/Linux environment. Other OS have not been tested and are not likely to be compatible with the software.    

The recommended way is to use a *conda* (*anaconda*, *miniconda*) environment to install the software and its dependencies.     

Here is a step-by-step installation guide using a *conda* environment:     
1.	If not yet done, install *conda*:     
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html     
To test your installation, in your terminal, type `conda list`     
If you read “bash: conda: command not found...”, then the conda installation was not successful and LORCAN cannot be installed with conda.     
Otherwise, you should see some text in the console, starting with "#packages in environment at...".       


There are two options to install the dependencies (3A: facilitated, 3B: manual).    

3A.	Install all dependencies at once using the *yaml* file provided:     
```
conda env create -f lorcan.yml
```
Note: the name of the environment can be modified in the yaml file.

3B.	Create a dedicated environment (to do only once) and install the dependencies one by one in the environment:     
```
conda create -n LORCAN-env         
```

Install the dependencies in the activated environment (the pipeline was only tested with the versions listed below):      
```
source activate LORCAN-env         
conda install seqkit=0.8.0       
conda install porechop=0.2.3_seqan2.1.1
conda install mafft=7.407
conda install Gblocks=0.91b           		
conda install iqtree=1.6.11	             	
conda install minimap2=2.17	       
conda install samtools=1.9       
conda install blast=2.9.0       
cpanm Parallel::ForkManager@2.02 	       
cpanm FindBin@1.51
conda install -c anaconda gcc_linux-64=7.3.0
cpanm PDF::API2@2.034       
```

To save some space, optionally remove unneeded files:    
```
conda clean --tarballs; conda clean --package        
```


4.	Clone the *LORCAN* binaries from GitHub to your own machine:          
(create and move first to the directory where you want to install LORCAN)      
```
git clone https://github.com/aramette/LORCAN.git    
```

5.	Test whether the pipeline is correctly installed:    
```
cd LORCAN/1.8/lib
perl lorcan.pl -h
```
At this stage, you should see the *help* menu being listed in your console. Otherwise, the download was not successful.        


6. Optional, but recommended:      
To be able to run the command everywhere on your system with the current username, edit your `~/.bashrc` file by running the following command that adds an alias to the .bashrc file:           
```
PWD=`pwd`    # /pathto/LORCAN/1.8/lib
printf "alias LORCAN=\'perl $PWD/lorcan.pl\'">> ~/.bashrc

```
Then run:    
```
source ~/.bashrc
```
to take the modification into account (next time your start a session, the changes will already be written in your .bashrc file).      
now you can call the program directly using the alias *LORCAN* on the command line:         
```
LORCAN -h
```

## CONFIGURATION         
<P ALIGN=RIGHT><a href="#top">↥ back to top</a></P>          

**1) General pipeline configuration**     
This step needs to be done only once when installing LORCAN.       
Edit the lib/config.pm file in the LORCAN/lib folder with your favorite text editor.      
This file contains the absolute paths to the dependencies (executable, databases).    
Each field starting with "our" may need to be changed depending on your own configuration and system. Make sure you exactly write the paths as indicated in the config.pm template.       
e.g.        
our $LIBDIRECTORY	="/path_to/LORCAN/1.8/lib"; 

If not yet done, activate the environment:         
```
source activate LORCAN-env 
```
Then edit the path description for each executable. 
To check for the exact paths, type in the console:       
```
which perl
which seqkit
which porechop
which minimap2
which samtools
which cutadapt
which blastn 
which Rscript
which mafft 	
which Gblocks
which iqtree
```
Then report the paths to the *lib/config.pm* file. An example file is provided in *lib/*.     



**2) Assay-specific configuration**    
<P ALIGN=RIGHT><a href="#top">↥ back to top</a></P>          

For each assay (i.e. 16S amplicon based on a specific set of primers), the following configuration needs to be done once.     
Here we will use the LeBiBi database as an example for 16S rRNA gene sequences.        
Example files are provided in the folder *example_files/*.      

Copy the Custom16S.tar.gz file to the LORCAN/DB/16S/ folder, and uncompress the tar.gz file:   
```
mkdir -p LORCAN/DB/16S/
cd LORCAN/DB/16S/
cp ../../1.8/example_files/Custom16S.tar.gz .    
tar xvzf Custom16S.tar.gz 
```
In the subfolder, you should see two example files:   
- 16S_stringent_custom.fasta  
- 16S_stringent_custom_simplified.names   

**Prepare the corresponding BLAST database**      
<P ALIGN=RIGHT><a href="#top">↥ back to top</a></P>          

Please refer to the document `Preparation of custom 16S database`[https://github.com/aramette/LORCAN/blob/master/1.8/Create_custom_db/custom_16S.md]  for further indications.    

```
makeblastdb -in  16S_stringent_custom.fasta -dbtype nucl  
```

Now indicate the path to the fasta file in the lib/config.pm file:    

```
%RefDB=(
	My16SDB => "/home/alban/test_software/LORCAN_install_test/LORCAN/DB/16S/BiBi/16S_stringent_custom.fasta", #<-------------
	MyOtherGeneDB => "/storage/databases/adenovirus/2018_10_10/Human_Adenovirus_genomes_20181010.fasta",
);
```
**Note**: Additional databases can be also listed in the config.pm file (separated by a comma; see example above). The name of the database (e.g. *"My16SDB"*) is what needs to be specified after the *-d* flag  when using *lorcan.pl* script).          


Finally, prepare the **TaxDictFile** (Taxonomy Dictionary file) as follows:      
The TaxDictFile is a comma-separated-values file with 3 fields: **`Field_1,Field_2,Field_3`**      
with:       

|Field| Description|
|:-----:|:----------|
|**Field_1**| Accession number (or unique id) of the sequence |
|**Field_2**| Taxonomic_grouping with no space (e.g. Pseudomonas_aeruginosa). This information is used to regroup reference sequences by level during the LORCAN analysis |
|**Field_3**| The full fasta header from the corresponding reference fasta file |

e.g.
```
URS00000B1AF5,Abiotrophia_defectiva,Abiotrophia_defectiva\~v\~TT\~URS00000B1AF5=Bacteria-Firmicutes-Bacilli-Lactobacillales-Aerococcaceae-Abiotrophia-Abiotrophia_defectiva       
```

Here is a small command-line script to extract and format the necessary information to creeate a taxdict.csv file:       
```
grep ">"  16S_stringent_custom.fasta | sed 's/>//' | awk -F\~ '{print $4","$1","$1"~"$2"~"$3"~"$4}' > My16S_taxdict.csv
```

*Field_2* of taxdict.csv can be further edited to keep the subspecies information or not.          
To print Field_2 with *subsp* information:    
```
cat My16S_taxdict.csv | cut -f2 --delim=","  | grep "subsp" | sort | uniq

## output:
# Acidovorax_avenae_subsp._avenae
# Acinetobacter_calcoaceticus_subsp._anitratus
# Actinobacillus_equuli_subsp._equuli
# Actinobacillus_equuli_subsp._haemolyticus
# …
```

Now indicate the path to the taxdict.csv file in the lib/config.pm file:           
```
pwd
our $TaxDictFile="/home/alban/test_software/LORCAN_install_test/LORCAN/DB/16S/BiBi/My16S_taxdict.csv";
```
The *config.pm* file is now properly edited.      



## TEST                  
<P ALIGN=RIGHT><a href="#top">↥ back to top</a></P>          

The test data (in example_files/input_test/) consists of a FASTQ file containing 8000 sequences from several barcoded 16S amplicon samples, with read length of 635.8 bp on average, ranging from 14-1,456 bp, and 5,086,542 bases in total.      

Uncompress the file:
```
gunzip fastq_0.fastq.gz   
```   

A sample_id text file should be prepared. One example is provided in LORCAN/1.8/example_files.          

Run LORCAN with the following parameters:      
```
WD=/path_of_your working_directory
EF=/pathto/LORCAN/1.8/example_files
FastqDir=$EF/input_test/

LORCAN -V -i $FastqDir -o $WD/myOutput -L log_main.txt -I $EF/sample_id.txt -n 20 -m 10 -M 3000 -P 100 -D 5 -d My16SDB
```
The results will be available in *$WD/myOutput* directory after 2-3 min or so.       

The PDF reports for each sample can be found in *myOutput/3_PDF_reports/*. They should match the information in the table below:         

| Sample ID       | Input reads | Consensus sequence |
| :-------------- |:-----------:| :---------------------------------------------------------|
| BC16:mysample1  | 148 | L=462 N=1 G=1 Taxo=Mycobacterium_chimaera
| BC17:mysample2  |  54 | Not enough reads matching. No consensus was built for this sample.
| BC18:mysample3  |  96 | Not enough reads matching. No consensus was built for this sample.
| BC19:mysample4  | 239	| L=462 N=0 G=1 Taxo=Mycobacterium_intracellulare_subsp._intracellulare
| BC22:mysample5  |	127	| No taxonomic group obtained in total >100 reads. No consensus was built for this sample.
| BC23:mysample6  |	168	| No taxonomic group obtained in total >100 reads. No consensus was built for this sample.
| BC24:mysample7  |	232	| L=462 N=1 G=1 Taxo=Mycobacterium_intracellulare_subsp._intracellulare	
| BC25:mysample8  |	191	| L=462 N=0 G=1 Taxo=Mycobacterium_intracellulare_subsp._intracellulare	
| BC26:mysample9  |	206	| No taxonomic group obtained in total >100 reads. No consensus was built for this sample.
| BC27:mysample10 | 142	| No taxonomic group obtained in total >100 reads. No consensus was built for this sample.
| BC28:mysample11 | 246	| L=483 N=0 G=0 Taxo=Streptococcus_thermophilus
| BC29:mysample12 | 272	| L=448 N=0 G=4 Taxo=Corynebacterium_tuberculostearicum
| BC30:mysample13 | 392	| L=448 N=0 G=4 Taxo=Corynebacterium_tuberculostearicum

- Input read = Total number of reads in the input fasta file after demultiplexing
- L = length (bp) of the consensus sequence produced by LORCAN 
- N = number of unknown bases, i.e. not  A/T/G/C in the consensus sequence ; 
- G = number of bases present in the reference sequence, but absent in the consensus sequence (deletions)
- Taxo = Taxonomic assignment based on read mapping, BLAST similarity analysis, and phylogenetic tree positioning


<P ALIGN=RIGHT><a href="#top">↥ back to top</a></P>          

The /LORCAN/myOutput folder should contain the following structure:    
```
 myOutput/
 |_ 0_logs/
 | 	 |_ log_main.txt		 (overall run log)       
 |	 |_ 1_initial_fastq_stats.txt 	 (fastq stats before demultiplexing)       
 |	 |_ 2_stats_fastq_postporechop.txt   (fastq stats after demultiplexing)       
 |	 |_ 3_Barcodes_with_modal_sequences.txt (BC number with modal sequences)       
 |_ 1_fasta/ 		  (for each barcode, 3 files are produced)       
 |	 |_ BC30.fasta		           (FASTA sequences from basecalled FASTQ)       	
 |	 |_ BC30.fasta_mode_closest.fasta  (FASTA sequences the closest to the modal length)
 |	 |_ BC30_stats.txt		   (Sequence statistics: total input, mode closest, etc.)     
 |_ 2_individual_barcodes/					 
 |	|_ BC30/    (for each BC, one folder is created)      
 |	    |_ 01_mapping_reads_to_refDB/  (all files used for mapping all mode closest reads to the chosen reference database)      
 |	    |	   |_ aln.bam		(only if -V flag is used)
 |	    |	   |_ aln.prim.bam	(only if -V flag is used)
 |	    |	   |_ aln.prim.bam.bai	(only if -V flag is used)
 |	    |	   |_ aln.prim.bam.stats  (only if -V flag is used)
 |	    |	   |_ aln.prim.sam	(only if -V flag is used)
 |	    |	   |_ aln.sam		(only if -V flag is used)
 |	    |	   |_ aln.sorted.bam	(only if -V flag is used)
 |	    |	   |_ list.mapped	(only if -V flag is used)
 |	    |	   |_ list.mappedS	(only if -V flag is used)
 |	    |	   |_ list.mapped.taxo	
 |	    |	   |_ new.prim.bam	(only if -V flag is used)
 |	    |	   |_ new.prim.bam.bai	(only if -V flag is used)
 |	    |	   |_ new.sam		(only if -V flag is used)
 |	    |_ 02_individual_consensus/	  (for each taxonomic level with enough reads, a consensus sequence is derived)
 |	    |	   |_ Corynebacterium_tubercuostearicum_consensus/  (one folder per taxonomic level with alignment files)      
 |	    |	   |        |_ aln1.bam		(only if -V flag is used)
 |	    |	   |        |_ aln1.mpileup	(only if -V flag is used)
 |	    |	   |        |_ aln1.mpileup.parsed	(only if -V flag is used)
 |	    |	   |        |_ aln1.sam			(only if -V flag is used)
 |	    |	   |        |_ aln1.sorted.bam		(only if -V flag is used)
 |	    |	   |        |_ BC30_ref_Corynebacterium_tubercuostearicum_consensus.fasta 
 |	    |	   |        |_ ConsensusPlot.html		
 |	    |	   |        |_ list_Corynebacterium_tuberculostearicum 
 |	    |	   |        |_ new1.bam	(only if -V flag is used)
 |	    |	   |        |_ new1.fasta	(only if -V flag is used)
 |	    |	   |        |_ R7g.Rmd		(only if -V flag is used)
 |	    |	   |        |_ REF1.fasta				
 |	    |	   |        |_ REF1.fasta.fai			
 |	    |	   |_ BC30_consensus_sequences.fasta (this file lists all consensus sequences that could be derived for the barcoded sample)
 |	    |_ 03_BLAST_analysis/  (BLASTN analysis of the consensus sequences against the reference database)     
 |	    |		    |_ blast_alignment.txt
 |	    |			    |_ blast_nt_output_3.txt
 |	    |			    |_ blast_nt_output_7_new.txt
 |	    |			    |_ blast_nt_output_7.txt
 |	    |_ 04_phylogenetic_analysis/  (all consensus sequences are used with their closest BLASTN hits)	
 |	    |	   |_ 01_alignment/
 |	    |	   |       |_ mafft.fasta
 |	    |	   |       |_ mafft.fasta-gb
 |	    |	   |       |_ mafft.fasta-gb.htm
 |	    |	   |       |_ mafft.fasta-gb.uniqueseq.phy
 |	    |	   |       |_ mafftWithHeaders.clw
 |	    |	   |_ 02_phylogeny/	(Phylogenetic analysis results produced by iqtree)
 |	    |	   |       |_ mafft.fasta-gb.bionj   (NJ tree)    
 |	    |	   |       |_ mafft.fasta-gb.contree  (Ultrafast bootstrap approximation, Consensus tree in Newick format)    
 |	    |	   |       |_ mafft.fasta-gb.iqtree  (IQ-TREE report)    
 |	    |	   |       |_ mafft.fasta-gb.log	(Screen log file)    
 |	    |	   |       |_ mafft.fasta-gb.mldist  (Likelihood distances)    
 |	    |	   |       |_ mafft.fasta-gb.treefile (Maximum-likelihood tree)    
 |	    |	   |       |_ README_file_info.txt
 |	    |	   |_ AllBlastREFs.fasta	(only if -V flag is used)
 |	    |	   |_ AllBlastREFsUniq.fasta	(only if -V flag is used)
 |	    |	   |_ CONSREF1.fasta
 |	    |	   |_ CONSREF1.header
 |	    |	   |_ CONSREF1_simple.fasta (only if -V flag is used)
 |	    |	   |_ OUTREFs.fasta	 (only if -V flag is used)
 |	    |_ 05_report/
 |	    |	   |_ BC30.report		(text file report)			
 |	    |	   |_ BC30_report.pdf	(PDF file report)
 |	    |_ 16S_stringent_custom.fasta	(reference database that was used ; can be deleted)
 |	    |_ 16S_stringent_custom.fasta.fai 	(reference database index that was used ; can be deleted)
 |	    |_ log.txt	(log for the specific barcode; here BC30)
 |_ 3_PDF_reports/	(copy of PDF reports from all barcode folders)
	    |_ BC30_report.pdf
```


To exit from the *LORCAN-env*, type: 
```
source deactivate LORCAN-env        
```

## FEEDBACK
Please file questions, bugs or ideas to the [Issue Tracker](https://github.com/aramette/LORCAN/issues).

## CITATION     
<P ALIGN=RIGHT><a href="#top">↥ back to top</a></P>          

xxxx.    


