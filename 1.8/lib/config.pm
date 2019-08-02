package config;

use 5.013;

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

use strict; 
#---------------------------------------------
$VERSION     = 1.08;
@ISA         = qw(Exporter);
@EXPORT      = qw(			
					$LIBDIRECTORY
					$PERL
					$SEQKIT
					$PORECHOP
					$minimap2
					$samtools
					$BLASTN
					$Rscript
					$mafft
					$Gblocks
					$iqtree
					$TaxDictFile
					%RefDB
			
);   		
#########################################################################
#########################################################################
#  Edit as needed the paths below as needed (see README.MD for instructions):
############################################# 
our $LIBDIRECTORY ="/home/alban/test_software/LORCAN_install_test/LORCAN/1.8/lib"; #<-------
############
our $PERL	="/home/alban/.conda/envs/LORCAN-env/bin/perl";
our $SEQKIT	="/home/alban/.conda/envs/LORCAN-env/bin/seqkit";
our $PORECHOP	="/home/alban/.conda/envs/LORCAN-env/bin/porechop";		
our $minimap2	="/home/alban/.conda/envs/LORCAN-env/bin/minimap2";		
our $samtools	="/home/alban/.conda/envs/LORCAN-env/bin/samtools";		
our $BLASTN	="/home/alban/.conda/envs/LORCAN-env/bin/blastn";
our $Rscript 	="/software_conda/3/bin/Rscript";
our $mafft 	="/home/alban/.conda/envs/LORCAN-env/bin/mafft";
our $Gblocks	="/home/alban/.conda/envs/LORCAN-env/bin/Gblocks";
our $iqtree	="/home/alban/.conda/envs/LORCAN-env/bin/iqtree";
##########
our $TaxDictFile="/home/alban/test_software/LORCAN_install_test/LORCAN/DB/16S/BiBi/My16S_taxdict.csv";
##########
our %RefDB=(
	My16SDB => "/home/alban/test_software/LORCAN_install_test/LORCAN/DB/16S/BiBi/16S_stringent_custom.fasta", 
	MyOtherGeneDB => "/home/alban/databases/adenovirus/2018_10_10/Human_Adenovirus_genomes_20181010.fasta",
);
####################################################
1;
