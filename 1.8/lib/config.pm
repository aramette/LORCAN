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
					$cutadapt
					$BLASTN
					$Rscript
					$mafft
					$Gblocks
					$iqtree
					$TaxDictFile
					%RefDB

					
);   		

#--------------------------------------------- edit as needed the following paths
# 
our $LIBDIRECTORY	="/software/KM/Transfer/cron_lorcan/v1.8.1/lib"; #<------------------ 

our $PERL			="/software/KM/perl5/perlbrew/perls/perl-5.28.0/bin/perl";
our $SEQKIT			="/software_conda/3/pkgs/seqkit-0.8.0-0/bin/seqkit";
our $PORECHOP		="/software_conda/3/bin/porechop";		
our $minimap2		="module load MultipleSequenceAlignment/minimap2/2.5-r601-dirty;minimap2";		
our $samtools		="module load UHTS/samtools/1.4;samtools";		
our $cutadapt		="/software_conda/2/bin/cutadapt";		
our $BLASTN			="module load Blast/blast/2.6.0;blastn";
our $Rscript 		="/software/KM/R-3.5.0_KM_WGS/bin/Rscript";
our $mafft 			="/software_conda/3/bin/mafft";
our $Gblocks		="/software_conda/3/bin/Gblocks";
our $iqtree			="/software_conda/3/bin/iqtree";
our $TaxDictFile	="/software/KM/Transfer/cron_lorcan/v1.8/LeBiBiRef.taxdict.csv"; #created using CreateTaxoDict_LeBiBI16S.pl
######
our %RefDB=(
	BiBi16S => "/storage/Analyses/Pipelines/dev/LORCAN/v1.8/leBiBicustomDB/16S_stringent_custom.fasta", #<-------------
	BiBi16SLong => "/storage/databases/leBIBI/SSU-rDNA/002_long16S_Amplicon/16S_stringent_custom_long.fasta", #<-------------
	ADV => "/storage/databases/adenovirus/2018_10_10/Human_Adenovirus_genomes_20181010.fasta",
);
1;
