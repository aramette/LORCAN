LORCAN "LOng Read Consenus ANalysis" 	 
====================================
LORCAN v1.8

# DESCRIPTION:
Using nanopore reads in fastq format, the pipeline produces consensus sequences 
for each barcoded sample, based on best match to a reference database 

# PARAMETERS:
Input dir 		
		The path to the directory containing the fastq files (not demultiplexed)
		
Output dir path
		The path to the output directory. The output directory will be named then "output".
		If the "output"	directory already exists, the START button will be disabled.
		
Number of cpus
		The max value is the value available on the current system.
		
Minimum number of reads per barcoded sample (threshold)
		The number of reads per taxonomic group under which the consensus per taxonomic group 
		will not be calculated.
		
Databases
		A database flat file (fasta) needs to be obtained and formatted to be compatible with BLASTN
			This is an example of formating for SSU-rDNA-mk37/16S_TS-stringent_dedup.fasta (leBiBi) database
				module load Blast/blast/2.6.0 #system specific
				makeblastdb -in  /path2dir/SSU-rDNA-mk37/16S_TS-stringent_dedup.fasta -dbtype nucl
			The database path should be indicated in the lib/config.pm file. 
		A corresponding .taxdict.csv file must be prepared using the same database file (fasta). 
			For each entry of the database file,the .taxdict.csv file contains the following fields:
				[Accession number],[taxonomic classication],[full fasta header information]
			e.g. URS00000B1AF5,Abiotrophia_defectiva,Abiotrophia_defectiva~v~TT~URS00000B1AF5=Bacteria-Fir...
			
			A utility script can be found in lib/CreateTaxoDict_LeBiBI16S.pl to prepare such file 
			using the LeBiBiref database.
			
		Other databases in the GUI pull-down menus can be added by modifying the config.pm file appropriately.
		
Email address
		Under normal mode, an email is sent to the provided email address when the job is complete.
		When verbose mode is on, no email is sent.

Verbose mode
		When on, more detailed output is provided. When off, several intermediary files are removed.
		
START
		This launches the script main.pl with all passed parameters. In case the output directory already 
		exists, the start button will be disabled. The user may want to choose another output directory path.		
		
		
CONFIGURATION:
==============
config.pm
		This file contains the absolute paths to the dependencies (executable, databases).
		Each field may need to be changed depending on your own configuration and system.

Perl Tk		
		To start LORCAN with the GUI, the Tk module of Perl must be installed. Otherwise, 
		the following error message may be seen:
			>perl lib/LORCAN_GUI.pl
				Can't locate Tk.pm in @INC ...
		The utility script "start_LORCAN_GUI.sh" may be modified according to the local Perl exe path.
		
		Perl Tk is not required if the command line is used instead of the GUI 
			>perl lib/main.pl
							Some input parameters are missing:
				 ->     NberThreads:    ....?
				 ->     readDIR:        ....?
				 ->     MinRead_threshold: .?
				 ->     Cdir:           ....?
				 ->     LOG:            ....?
				 ->     ReadThreshold:  ....?
				 ->     DBchoice:       ....?
				(use -h flag to get help with usage)

GUI	
		to obtain the right list of available databases on the GUI,
			the script LORCAN_GUI.pl under "=== choice db" should be modified according to the input in the config.pm file
				e.g. in config.pm we can read
				our %RefDB=(
						BiBi16S => "/storage/databases/leBIBI/SSU-rDNA/002/SSU-rDNA-mk37_stringent/16S_stringent_dedup.fasta",
						ADV => "/storage/databases/adenovirus/2018_10_10/Human_Adenovirus_genomes_20181010.fasta",
				);
				in LORCAN_GUI.pl we should also have:
				-options => [[my16Sdb=>"BiBi16S"], [Adenovirusdb=>"ADV"], ],
				with the name on the left side of the arrows being freely customable (should start by a letter character), 
				while on the right side, the names should match exactly what has been entered in the %RefDB hash definition.
		
