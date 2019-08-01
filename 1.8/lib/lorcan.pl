#!/usr/bin/env perl
#===============================================================================
#
#         FILE: main.pl
#
#
#  DESCRIPTION: main script file of LORCAN
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#       AUTHOR: Alban Ramette, alban.ramette@ifik.unibe.ch
# ORGANIZATION: IFIK
my      $version="1.8";
#      CREATED: 12/24/2018 
#      REVISED: 01/30/2019 
#         BUGS: ---
#         TODO: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use feature qw(say);
use Parallel::ForkManager;

#-------------
use FindBin;  # locate this script
use lib "$FindBin::Bin/../lib";
use config;
use LORCAN_subs;
#------------- 
use Cwd;
my $cwdir = getcwd;
#print "cwd: $cwdir\n"; # where the script was run
#-------------- flags
use Getopt::Std;

my %args;
getopts('n:i:I:m:E:o:L:P:vVhd:r:D:M:', \%args);

our ($NberThreads,$readDIR,$BCReadThreshold,$EMAIL,$Cdir,$LOG,$ReadThreshold);
our ($DBchoice,$MainLog,$Verbose);
our (@BCnamesOK,$BCfasta,$delta,$MaxRead,$sampleID);
our ($Time);

my $help= <<HELP;
Version: $version

What it does:
	demultiplex ONT fastq files,
	applies porechop,
	extract the modal sequences, 
	map the sequences to db, 
	and calculate specific consensus, 
	BLAST back the consensus to db,
	create report

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
	perl lib/main.pl -V -i \$testDataDir -n 20 -m 10 -o \$cwd/output_test_main_16S_delta5 -I \$cwd/sample_id.txt -L log_main.txt -P 30 -D 5 -d BiBi16S
HELP

if ($args{h}) {
	say $help; exit;

}else{ 
	if ($args{V}) {$Verbose=$args{V}; chomp $Verbose;}
	
	if ($args{i}) { #directory containing the fastq files;
       #if ($args{V}) {print "-i has a value of: $args{i} (directory containing the fastq files)\n";}
		$readDIR=$args{i}; chomp $readDIR;
	}
	if ($args{I}) { #the full path to the sample id file
       #if ($args{V}) {print "-I has a value of: $args{I} (directory containing the fastq files)\n";}
		$sampleID=$args{I}; chomp $sampleID;
	}
	if ($args{n}) { #NberThreads;
		#if ($args{V}) {print "-n has a value of: $args{n} (Nber threads)\n";}
		$NberThreads=$args{n}; chomp $NberThreads;
	}
	if ($args{m}) { # minimum number of reads per barcode to further proceed 
			#if ($args{V}) {print "-m has a value of: $args{m} (minimum number of reads per barcode to further proceed )\n";}
			$BCReadThreshold=$args{m}; chomp $BCReadThreshold;
	}
	if ($args{M}) { # Maximum number of reads to retain per sample
			$MaxRead=$args{M}; chomp $MaxRead;
	}	
	if ($args{P}) { # minimum  read threshold for a reference to be further considered
			#if ($args{V}) {print "-P has a value of: $args{P} ($ReadThreshold is the minimal read nber per taxonomic level to keep a specific ref sequence )\n";}
			$ReadThreshold=$args{P}; chomp $ReadThreshold;
	}
	if ($args{D}) { # # +- delta bp from the mode
			$delta=$args{D}; chomp $delta;
	}

	if ($args{E}) { # 
			#if ($args{V}) {print "-E has a value of: $args{E} (email address)\n";}
			$EMAIL=$args{E}; chomp $EMAIL;
	}
	if ($args{o}) { #e.g. "/storage/Analyses/Pipelines/dev/NABA_GUI/v2.5/SUBSET/output"
			#if ($args{V}) {print "-o has a value of: $args{o} (output directory)\n";}
			$Cdir=$args{o}; chomp $Cdir;
	}
	if ($args{L}) { #e.g. "log.txt"
			#if ($args{V}) {print "-L has a value of: $args{L}\n";}
			$LOG=$args{L}; chomp $LOG;
	}	
	if ($args{v}) { #e.g. 
			say $version;exit;
	}
	if ($args{d}) { #e.g. BiBi16S, ADV
			#if ($args{V}) {print "-d has a value of: $args{d}\n";}
			$DBchoice=$args{d}; chomp $DBchoice;
	}	
}

######################################################################
#### MAIN ############################################################
if ($NberThreads && $readDIR &&  $BCReadThreshold && $Cdir && $LOG && $ReadThreshold && $DBchoice && $delta && $MaxRead){
	if ($args{V}) {say "OK: all arguments are present";print "current working directory: $cwdir\n";}
} else {say "Some input parameters are missing:";
	if($NberThreads){say "\tNberThreads: \t$NberThreads";}			else{say " ->\tNberThreads: \t....?";}
	if($readDIR){say "\treadDIR: \t$readDIR";}  					else{say " ->\treadDIR: \t....?";}
	if($BCReadThreshold){say "\tBCReadThreshold:  $BCReadThreshold"; }  else{say " ->\tBCReadThreshold: .?";}
	if($Cdir){say "\tCdir: \t\t$Cdir";}   								else{say " ->\tCdir: \t\t....?";}
	if($LOG){say "\tLOG: \t\t$LOG"; }  									else{say " ->\tLOG: \t\t....?";}
	if($ReadThreshold){say "\tReadThreshold: \t$ReadThreshold"; }  	else{say " ->\tReadThreshold: \t....?";}
	if($DBchoice){say "\tDBchoice: \t$DBchoice"; } 					else{say " ->\tDBchoice: \t....?";}
	if($delta){say "\tdelta: \t$delta"; } 					else{say " ->\tdelta: \t....?";}
	if($MaxRead){say "\tMaxRead: \t$MaxRead"; } 					else{say " ->\tMaxRead: \t....?";}
	if($sampleID){say "\tsampleID: \t$sampleID"; } 					else{say " ->\tsampleID: \t....?";}
	say "(use -h flag to get help with usage)";
	exit;
}

unless(-e "$Cdir/0_logs") {system("mkdir -p $Cdir/0_logs");	}	
$MainLog = "$Cdir/0_logs/$LOG";
open MAINLOG,">>$MainLog"; # generate and open report file; do not overwrite
print MAINLOG "=========================================================\n";
print MAINLOG `date`;


chdir $Cdir;	
if($Verbose){
	say "\t(i) readDIR:           $readDIR";  
	say "\t(I) sampleID:          $sampleID";  
	say "\t(o) Cdir:              $Cdir";  
	say "\t(L) Main log name:     $LOG";  
	say "\t(n) NberThreads:       $NberThreads"; 
	say "\t(m) Min nber reads per BC: $BCReadThreshold";  
	say "\t(M) Consensus MaxRead: $MaxRead";  
	say "\t(P) Consensus MinRead: $ReadThreshold";  
	say "\t(D) delta:             $delta";  
	say "\t(d) DBchoice:          $DBchoice";  
}

	say MAINLOG "   (i) readDIR:         		$readDIR";  
	say MAINLOG "   (I) sampleID:         		$sampleID";  
	say MAINLOG "   (o) Cdir:           		$Cdir";  
	say MAINLOG "   (L) Main log name:   		$LOG";  
	say MAINLOG "   (n) NberThreads:     		$NberThreads"; 
	say MAINLOG "   (m) Min nber reads per BC: 	$BCReadThreshold (minimum number of reads per barcode to further proceed)";  
	say MAINLOG "   (M) Consensus MaxRead:      $MaxRead (CONS: maximum number of reads to retain per sample)";  
	say MAINLOG "   (P) Consensus MinRead:   	$ReadThreshold (CONS: minimum  number of read aligned to a reference for the latter to be further considered)";  
	say MAINLOG "   (D) delta:           		$delta";  
	say MAINLOG "   (d) DBchoice:        		$DBchoice";  
	say MAINLOG "";
	
my ($RunName,$RunDate,%BCinfo); #info from the sample_id.txt file
if(-e $sampleID) {
	open ID,"<$sampleID";
	foreach (<ID>){
		chomp;
		if(/^RunName\:(.+)$/){$RunName=$1; $RunName=~s/\s+//g;}
		if(/^RunDate\:(.+)$/){$RunDate=$1; $RunDate=~s/\s+//g;}
		if(/^(BC\d+)\:(.+)$/){$BCinfo{$1} = $2 if defined $2;}
	}
	
}	
	#checking	
	 #say "RunName=$RunName";
	 #say "RunDate=$RunDate";
	 my $Hashsize = scalar(keys %BCinfo);
	 foreach my $R (keys %BCinfo) {
		 say MAINLOG "$R => $BCinfo{$R}";   
	 }	
	 if($Hashsize == 0){say "sample id file was not parsed correctly. Use BCii:sampleinfo:db for the descriptor. Exiting";exit;}
 

#R0--------------------
say MAINLOG "\n=>0) Checking script dependencies";
unless (-e "$Cdir/0_logs/$LOG") {
	if($Verbose){say "=> Checking script dependencies..."; }
	$Time=localtime; 
	say MAINLOG " started: ...\t $Time";												
	LORCAN_subs::R0_test_dependencies({	
			Verbose 	 => $Verbose,
			LIBDIRECTORY => $LIBDIRECTORY,
			PERL 		 => $PERL,	
			SEQKIT 		 => $SEQKIT,
			PORECHOP 	 => $PORECHOP,
			minimap2 	 => $minimap2,
			samtools 	 => $samtools,
			cutadapt 	 => $cutadapt,
			BLASTN 		 => $BLASTN,
			Rscript 	 => $Rscript,	
			REF_FP 		 => $RefDB{$DBchoice}, # %RefDB defined in config.pm file
	});
 	$Time=localtime; say MAINLOG " finished: ...\t $Time";
} else {if($Verbose){say "=> (R0_Checking script dependencies) already completed"; }}

#R1--------------------
say MAINLOG "\n=>1) Running R1_Seqkit";
unless (-e "$Cdir/0_logs/1_initial_fastq_stats.txt") {
	if($Verbose){say "=> running R1_Seqkit..."; }
	$Time=localtime; say MAINLOG " started: ...\t $Time";															
	LORCAN_subs::R1_Seqkit({
		Verbose 		=> $Verbose,
		IN 				=> $readDIR,
		OUT 			=>"$Cdir/0_logs/1_initial_fastq_stats.txt",
		NberThreads 	=> $NberThreads,
		SEQKIT 			=> $SEQKIT,
	});
	$Time=localtime; say MAINLOG " finished: ...\t $Time";	
	unless (-e "$Cdir/0_logs/1_initial_fastq_stats.txt") {say "!!! Missing file (exit!): \n\t$Cdir/0_logs/1_initial_fastq_stats.txt";exit}	
} else {if($Verbose){say "=> (R1_Seqkit) already completed"; }}

#R2--------------------
say  MAINLOG "\n=>2) Running R2_porechop\n";
$BCfasta=$Cdir."/1_fasta"; mkdir $BCfasta;			
my $R2_test=qx(ls $Cdir/1_fasta/ | grep -c "fasta"); chomp $R2_test;
if($R2_test == 0){
	if($Verbose){say "=> running R2_porechop..."; }
	$Time=localtime; say MAINLOG " started: ...\t $Time";
	LORCAN_subs::R2_porechop({
		Verbose => $Verbose,
		IN => $readDIR,
		OUT =>"$Cdir/1_fasta/",
		PORECHOP =>$PORECHOP,
		NberThreads => $NberThreads,
	});
	$Time=localtime; say MAINLOG " finished: ...\t $Time";
} else {if($Verbose){say "=> (R2_porechop) already completed"; }}

my $R2_test1=qx(ls $Cdir/1_fasta/ | grep -c "fasta"); chomp $R2_test1;
if($R2_test1 == 0) {say "!!! No fasta file generated in $BCfasta (exit!)";exit}


#R3--------------------
print MAINLOG "\n=>3) Running R3_Seqkit_after_porechop\n";

unless (-e "$Cdir/0_logs/2_stats_fastq_postporechop.txt") {
	if($Verbose){say "=> running R3_Seqkit_after_porechop..."; }
	$Time=localtime; say MAINLOG " started: ...\t $Time";
		LORCAN_subs::R1_Seqkit({
		Verbose => $Verbose,
		IN => "$Cdir/1_fasta/",
		OUT =>"$Cdir/0_logs/2_stats_fastq_postporechop.txt",
		NberThreads => $NberThreads,
		SEQKIT => $SEQKIT,
	});
	$Time=localtime; say MAINLOG " finished: ...\t $Time";
	unless (-e "$Cdir/0_logs/2_stats_fastq_postporechop.txt") {say "!!! Missing file - R3 test0 (exit!): \n\t$Cdir/0_logs/2_stats_fastq_postporechop.txt";exit}	
} else {if($Verbose){say "=> (R3_Seqkit_after_porechop) already completed"; }}
my $R3_test1=qx(cat $Cdir/0_logs/2_stats_fastq_postporechop.txt| wc -l ); chomp $R3_test1;
if($R3_test1 == 0) {say "!!! Empty file generated at R3_test: \n$Cdir/0_logs/2_stats_fastq_postporechop.txt (exit!)";exit}

#R4--------------------
print MAINLOG "\n=>4) Running R4_Parsing stats for min number of reads per barcode\n";
print MAINLOG "chosen minimum threshold = $BCReadThreshold reads\n";
my $R4_test= (-d "$Cdir/2_individual_barcodes"); 

unless($R4_test){
	if($Verbose){say "=> running R4_parseMinNberReads..."; }
	$Time=localtime; say MAINLOG " started: ...\t $Time";
	LORCAN_subs::R4_parseMinNberReads({
		Cdir => $Cdir,
		BCReadThreshold => $BCReadThreshold,
		IN => "$Cdir/0_logs/2_stats_fastq_postporechop.txt",
		OUT  => "$Cdir/0_logs/$LOG"
	});
	$Time=localtime; say MAINLOG " finished: ...\t $Time";
	$R4_test= qx(grep -c "Created barcode directories accordingly" $MainLog); chomp $R4_test;
	if ($R4_test == 0) {say "!!! No barcoded sample remaining (exit!): R4_parseMinNberReads";exit}	
} else {if($Verbose){say "=> (R4_parseMinNberReads) already completed"; }}

our @BCDIR=qx(ls $Cdir/2_individual_barcodes/);chomp @BCDIR;
my $R4_test2= scalar @BCDIR; 
unless ($R4_test2>0) {say "R4_test2 failed. No directory created (exit)"; exit;}

#R5--------------------
print MAINLOG"\n=>5) running R5_Extract_modal_sequences\n";
say MAINLOG "(produce the mode +-delta bp sequences)";

our $BCDIRref = \@BCDIR;
unless(-e "$Cdir/0_logs/3_Barcodes_with_modal_sequences.txt"){
	if($Verbose){say "=> running R5_Extract_modal_sequences..."; }
	$Time=localtime; say MAINLOG " started: ...\t $Time";
	LORCAN_subs::R5_Extract_modal_sequences({
		Verbose 	=> $Verbose,
		IN 			=> $BCDIRref,#array reference to @BCDIR
		OUT  		=> "$Cdir/0_logs/$LOG",  #adding to the MAINLOG!!
		Cdir 		=> $Cdir,		
		NberThreads => $NberThreads,
		Rscript 	=> $Rscript,
		delta 		=> 	$delta,# +- delta bp from the mode
		MaxRead 	=> 	$MaxRead,# Maximum number of reads to retain per sample
		LIBDIRECTORY => $LIBDIRECTORY,
	});
	$Time=localtime; say MAINLOG " finished: ...\t $Time";
	my $filetest1=$BCDIR[0].".fasta_mode_closest.fasta";
	unless(-e "$Cdir/1_fasta/$filetest1"){ say "!!! No modal sequences produced for barcode sample $BCDIR[0] - test 1";}	
	my $R5_test=qx(cat $Cdir/0_logs/3_Barcodes_with_modal_sequences.txt | wc -l); chomp $R5_test;
	unless($R5_test >0){ say "!!! No modal sequences produced for barcode sample $BCDIR[0] - test 2"}	
} else {if($Verbose){say "=> (R5_Extract_modal_sequences) already completed"; }}

foreach(@BCDIR){ # checking production of files again
	my $filetest2=$_.".fasta_mode_closest.fasta";
	unless(-e "$Cdir/1_fasta/$filetest2"){ say "!!! No modal sequences produced for barcoded sample: $_ (exit!) - test 3";}	
}

#-------------------- analysis for each barcoded sample
print MAINLOG"\n=> Running Steps for each barcoded sample in parallel\n";
if ($args{V}) {	say "=> running Parallel steps for each barcode...";}

#-------------------------#only for debugging and development 
	#@BCDIR=("BC28","BC25");	#OK until PDF
	#@BCDIR=("BC16","BC17");	# do not produce PDF, stopping at list.mapped (empty file created)
	#@BCDIR=("BC28","BC16");	
	#@BCDIR=("BC28");	
#-------------------------
	my $NBc=scalar(@BCDIR);
	my $pm = new Parallel::ForkManager($NberThreads);
	if($NberThreads gt $NBc ){$pm = new Parallel::ForkManager($NBc);}
	my $pid;
	
########## parallel processing ######################################################
foreach my $BC (@BCDIR){ 
	$pid = $pm->start and next; 	
	say MAINLOG "=> Analysis of barcoded sample: $BC ========";
	my $INPUTFASTA="$Cdir/1_fasta/$BC.fasta_mode_closest.fasta";
	my $TotalNberReadsinFasta=qx( grep -c ">" $INPUTFASTA); chomp $TotalNberReadsinFasta;
	if($Verbose) {say "  =>$BC: Total Nber of modal fasta sequences: $TotalNberReadsinFasta";}
	say MAINLOG "  =>$BC: Total Nber reads in Fasta: $TotalNberReadsinFasta";
	unless(-e "$Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/new.prim.bam" && -e "$Cdir/2_individual_barcodes/$BC/list.mapped"){
		if($Verbose){say "  =>$BC: running R6_MapWithMinimap2..."; }
		say MAINLOG "  =>$BC: running R6_MapWithMinimap2...";
		$Time=localtime; say MAINLOG "$BC: started: ...\t $Time";
		LORCAN_subs::R6_MapWithMinimap2({
				Verbose 	=> $Verbose,
				IN 			=> $INPUTFASTA,
				OUT1 		=> "new.prim.bam",
				OUT2 		=> "list.mapped",
				Wdir 		=> "$Cdir/2_individual_barcodes/$BC",		
				REF_FP 		=> $RefDB{$DBchoice},
				PERL 		=> $PERL,	
				minimap2 	=> $minimap2,
				NberThreads => $NberThreads,
				samtools 	=>  $samtools,
		});		
		$Time=localtime; say MAINLOG "$BC: finished: ...\t $Time";
	}else {if($Verbose){say "  =>$BC: (R6_MapWithMinimap2) already completed"; }}
	$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;

say MAINLOG "  [",scalar(localtime),"] R6 complete\n";

###############################
#-------------------------R7_calc_stats_perContig
foreach my $BC (@BCDIR){#make a test of -z list.mapped and continue only with non empty ones... 
	$pid = $pm->start and next;	
	my $TOTALNberReads=qx($samtools view  $Cdir/2_individual_barcodes/$BC/new.prim.bam | wc -l);  chomp $TOTALNberReads;
	if ($Verbose){say "  =>$BC: Total number of reads that mapped after step R6: $TOTALNberReads"; }
			
	if (-z "$Cdir/2_individual_barcodes/$BC/list.mapped.taxo"){
			if($Verbose) {say "  =>$BC: !!! list.mapped.taxo is empty. No further processing of the sample.";} #*************** need to write a report on this
			say MAINLOG "  =>$BC: list.mapped.taxo is empty after step R6 (no enough reads mapped to reference). No further processing of the sample. ";
	} else{ #continue the pipeline
			unless(-e "$Cdir/2_individual_barcodes/$BC/log.txt"){
				unless (-e "$Cdir/2_individual_barcodes/$BC/list.mapped.taxo") {
					if($Verbose){say "  =>$BC: running MatchTaxo..."; }
					say MAINLOG "  =>$BC: running MatchTaxo...";
					LORCAN_subs::Matchtaxo({
							Verbose 	=> $Verbose,
							IN1 		=> "$Cdir/2_individual_barcodes/$BC/list.mapped",
							IN2 		=> $TaxDictFile,
							OUT 		=> "$Cdir/2_individual_barcodes/$BC/list.mapped.taxo",		
					});		
					unless (-e "$Cdir/2_individual_barcodes/$BC/list.mapped.taxo") {say "!!! missing file: $Cdir/2_individual_barcodes/$BC/list.mapped.taxo";}
				}else {if($Verbose){say "  =>$BC: (MatchTaxo) already completed"; }}	

			#-------------------------	
			chdir "$Cdir/2_individual_barcodes/$BC"; 
				if($Verbose){say "  =>$BC: running R7_calc_stats_perContig..."; }
					say MAINLOG "  =>$BC: R7_calc_stats_perContig...";
					$Time=localtime; say MAINLOG "$BC: started: ...\t $Time";
					LORCAN_subs::R7_calc_stats_perContig({
									Verbose 			=> $Verbose,
									IN0					=> "$Cdir/2_individual_barcodes/$BC",
									IN1 				=> "$Cdir/2_individual_barcodes/$BC/list.mapped.taxo",
									IN2 				=> "$Cdir/2_individual_barcodes/$BC/new.prim.bam",
									OUT1 				=> "$Cdir/2_individual_barcodes/$BC/log.txt",
									OUT2 				=> "aln1.mpileup.parsed",
									ReadThreshold 		=> $ReadThreshold,# $ReadThreshold is the minimal read nber per taxonomic level to keep a specific ref sequence
									BC 					=> $BC,
									SEQKIT 				=> $SEQKIT,
									minimap2 			=> $minimap2,
									samtools 			=> $samtools,
									REF_FP 				=> $RefDB{$DBchoice},
					});		
					unless (-e "$Cdir/2_individual_barcodes/$BC/list.mapped.taxo") {say "!!! missing file: $Cdir/2_individual_barcodes/$BC/list.mapped.taxo";}

				$Time=localtime; say MAINLOG "$BC: finished: ...\t $Time";
			}else {if($Verbose){say "  =>$BC: (R7_calc_stats_perContig) already completed"; }}	
	}
	$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;
say MAINLOG "  [",scalar(localtime),"] R7 complete\n";

########### move files to appropriate directory
foreach my $BC (@BCDIR){#make a test of -z list.mapped and continue only with non empty ones... 
	unless(-e"$Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB"){system("mkdir $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB");}
	system("mv $Cdir/2_individual_barcodes/$BC/aln* $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/");
	system("mv $Cdir/2_individual_barcodes/$BC/new* $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/");
	system("mv $Cdir/2_individual_barcodes/$BC/list* $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/");
}
#after this point, if folder 02_individual_consensus/ does not exist => not enough read for the BC

###############################R8_add consensus quality
# also produces for each $BC, 
# $Cdir/2_individual_barcodes/$BC/02_individual_consensus/$BC_consensus_sequences.fasta

if($Verbose){say "  =>All:  running R8_Consensus_format_extract..."; }
		say MAINLOG "  =>All:  R8_consensus_format_extract...";
		$Time=localtime; say MAINLOG "  started: ...\t $Time";
		LORCAN_subs::R8_Consensus_format_extract({
										Verbose 	=> $Verbose,
										IN 			=> $BCDIRref,
										OUT1 		=> "log.txt",
										OUT2 		=> "consensus_sequences.fasta",
										Cdir 		=> $Cdir,
		});
		$Time=localtime; say MAINLOG "  finished: ...\t $Time";
#-------------------------------
# Loop to identify those BC_conseensus files that are missing	
my @BCnoCons;
foreach my $BC (@BCDIR){
	my $Allconseqfilename="$BC\_consensus_sequences.fasta";
	if(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename")
	{ #say "$BC: 02_individual_consensus/BC_consensus_sequences.fasta exists";
	} else {
	  #say "$BC: 02_individual_consensus/BC_consensus_sequences.fasta not found";
	  push(@BCnoCons,$BC);
	}
}
foreach my $NBC(@BCnoCons){#remove not OK from BCDIR
	@BCDIR =grep(!/^$NBC/,@BCDIR);
}
#----------------------------R9_check_seq_blast_taxo
foreach my $BC (@BCDIR){
	$pid = $pm->start and next;	
	my $Allconseqfilename="$BC\_consensus_sequences.fasta";
	unless(-e"$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis"){system("mkdir $Cdir/2_individual_barcodes/$BC/03_BLAST_analysis");}
	if(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename"){ # 
			chdir "$Cdir/2_individual_barcodes/$BC"; 
				unless(-e "$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis/blast_nt_output_7.txt"){
				if($Verbose){say "  =>$BC: running R9_check_seq_blast_taxo..."; }
						say MAINLOG "  =>$BC: R9_check_seq_blast_taxo...";
						$Time=localtime; say MAINLOG "$BC: started: ...\t $Time";
						LORCAN_subs::R9_check_seq_blast_taxo({
							Verbose 		=> $Verbose,
							IN 				=> "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename",
							OUT 			=> "$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis",
							REF_FP 			=> $RefDB{$DBchoice},				
							NberThreads 	=> $NberThreads,
							PERL 			=> $PERL,	
							BLASTN 			=> $BLASTN,	
						});
						$Time=localtime; say MAINLOG "$BC: finished: ...\t $Time";
						unless(-e "$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis/blast_nt_output_7.txt"){say "$BC: blast_summary.txt (R9) not produced. ";}
				}else {if($Verbose){say "  =>$BC: (R9_check_seq_blast_taxo) already completed"; }}
	}
$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;				

#----------------------------R9_2_check_seq_phylogeny
foreach my $BC (@BCDIR){
	$pid = $pm->start and next;	
	my $Allconseqfilename="$BC\_consensus_sequences.fasta";
		unless(-e"$Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis"){system("mkdir $Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis");}
	if(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename"){ # 
			chdir "$Cdir/2_individual_barcodes/$BC"; 
	unless(-e "$Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/02_phylogeny/mafft.fasta-gb.treefile"){
				if($Verbose){say "  =>$BC: running R9_2_check_seq_phylogeny..."; }
						say MAINLOG "  =>$BC: R9_2_check_seq_phylogeny...";
						$Time=localtime; say MAINLOG "$BC: started: ...\t $Time";
						LORCAN_subs::R9_2_check_seq_phylogeny({
							Verbose 	=> $Verbose,
							IN1 		=> "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename",
							IN2			=> "$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis/blast_nt_output_7.txt",
							OUT 		=> "$Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis",
							NberThreads => $NberThreads,
							PERL 		=> $PERL,	
							SEQKIT 		=> $SEQKIT,
							mafft 		=> $mafft,
							Gblocks 	=> $Gblocks,
							iqtree 		=> $iqtree,
							REF_FP 		=> $RefDB{$DBchoice},								
						});
						$Time=localtime; say MAINLOG "$BC: finished: ...\t $Time";

						unless(-e "$Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/02_phylogeny/mafft.fasta-gb.treefile"){say "$BC: 02_phylogeny/mafft.fasta-gb.treefile (R9_2) not produced. ";}
				}else {if($Verbose){say "  =>$BC: (R9_2_check_seq_phylogeny) already completed"; }}
	}
$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;				


#----------------------------R10_createTxtReport $BC.report
foreach my $BC (@BCDIR){
	$pid = $pm->start and next;	
	#chdir "$Cdir/2_individual_barcodes/$BC";  	
	my $Allconseqfilename="$BC\_consensus_sequences.fasta";
	unless(-e"$Cdir/2_individual_barcodes/$BC/05_report"){system("mkdir $Cdir/2_individual_barcodes/$BC/05_report");}
	if(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename"){
				unless(-e "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report"){
					#reformat the Blast7 output
						open(BLAST7,"$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis/blast_nt_output_7.txt");
						open (NEWB7,">$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis/blast_nt_output_7_new.txt");
						foreach(<BLAST7>){
							if (/^#/){
								if (/^# Query/){print NEWB7 $_;}
								if (/^# Fields/){
									say NEWB7 "\# Fields\:\n\% id\/Nber identical\/alignment length\/Nber mismatches\/Nber gaps  subject_id";
								}
							} else{
							my $line=$_;
							$line =~ m/(\d+\.\d+)\t(\d+)\t(\d+)\t(\d+)\t(.+)\t(.+)/;
								my $ID=$1;
								my $NID=$2;
								my $AL=$3;
								my $MM=$4;
								my $GA=$5;
								my $SID=$6; chomp $SID;
								my $SID1=substr($SID, 0,80);my $SID2=$SID1."...";
								print NEWB7 "$ID/$NID/$AL/$MM/$GA $SID2\n";
							}
						}

				if($Verbose){say "  =>$BC: running R10_createTxtReport..."; }
						say MAINLOG "  =>$BC: R10_createTxtReport...";
						$Time=localtime; say MAINLOG "$BC: started: ...\t $Time";
						my $INPUTFASTA="$Cdir/1_fasta/$BC.fasta_mode_closest.fasta";
						my  $TotalNberReadsinFasta=qx( grep -c ">" $INPUTFASTA); chomp $TotalNberReadsinFasta;
						my $ID=$BCinfo{$BC};
						LORCAN_subs::R10_createTxtReport({
							Verbose 				=> $Verbose,
							BC 						=> $BC,
							IN1 					=> "$Cdir/2_individual_barcodes/$BC/log.txt",
							ID 						=> $ID,
							IN3 					=> "$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis/blast_nt_output_7_new.txt",
							IN4 					=> "$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis/blast_alignment.txt",
							IN5 					=> "$Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/CONSREF1.header",
							IN6 					=> "$Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/02_phylogeny/mafft.fasta-gb.iqtree",
							OUT 					=> "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report",
							TotalNberReadsinFasta 	=> $TotalNberReadsinFasta,
							REF_FP 					=> $RefDB{$DBchoice},				
							DBchoice 				=> $DBchoice,
						});
						$Time=localtime; say MAINLOG "$BC: finished: ...\t $Time";
						my $R10_test=(-s "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report");
						unless($R10_test != 0 || -e "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report"){
						say "$BC.report (R10) not produced. Exiting..."; exit;}
				}else {if($Verbose){say "  =>$BC: (R10_createTxtReport) already completed"; }}	
				
		}
$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;		

			
#----------------------------R11_CreatePDFfromTXT
foreach my $BC (@BCDIR){
	$pid = $pm->start and next;	
	chdir "$Cdir/2_individual_barcodes/$BC";  	
	my $Allconseqfilename="$BC\_consensus_sequences.fasta";
	if(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename"){
				my $PDFfile="$Cdir/2_individual_barcodes/$BC/05_report/$BC"."_report.pdf";
				unless(-e $PDFfile){
					if($Verbose){say "  =>$BC: running R11_CreatePDFfromTXT..."; }
					say MAINLOG "  =>$BC: R11_CreatePDFfromTXT...";
					$Time=localtime; say MAINLOG "$BC: started: ...\t $Time";
					LORCAN_subs::R11_CreatePDFfromTXT({
								Verbose => $Verbose,
								IN 		=> "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report",
								OUT 	=> $PDFfile,
					});
					$Time=localtime; say MAINLOG "$BC: finished: ...\t $Time";
					unless (-e $PDFfile)
						{say "$BC report.pdf (R12) not produced. Exiting..."; exit;}
				}else {if($Verbose){say "  =>$BC: (R11_CreatePDFfromTXT) already completed"; }}	
		}
$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;					

say MAINLOG "  [",scalar(localtime),"] R1-11 complete\n";
print MAINLOG"\n=> End of steps R1-11 for each barcoded sample in parallel\n";
###############################
# post-processing ---------------------------------------------------
#---------------------------- R12_final check for not completed samples and reporting as PDF too 
#
foreach my $BC (@BCnoCons){
	$pid = $pm->start and next;	
	chdir "$Cdir/2_individual_barcodes/$BC";  	
	unless(-e "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report"){
		if($Verbose){say "  =>$BC: running R12_final check for not completed samples";}
		unless(-e "$Cdir/2_individual_barcodes/$BC/05_report/"){system("mkdir $Cdir/2_individual_barcodes/$BC/05_report/");}
	my $INPUTFASTA="$Cdir/1_fasta/$BC.fasta_mode_closest.fasta";
	my  $TotalNberReadsinFasta=qx( grep -c ">" $INPUTFASTA); chomp $TotalNberReadsinFasta;
	my $ID=$BCinfo{$BC};
	LORCAN_subs::R11_createTxtReportNo({
							Verbose 				=> $Verbose,
							BC 						=> $BC,
							IN 						=> "$Cdir/2_individual_barcodes/$BC/log.txt",
							ID 						=> $ID,
							OUT 					=> "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report",
							TotalNberReadsinFasta 	=> $TotalNberReadsinFasta,
							REF_FP 					=> $RefDB{$DBchoice},				
							DBchoice 				=>  $DBchoice,
						});
	} else {if($Verbose){say "  =>$BC: (R12_final check for not completed samples) already completed";}}
$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;	

#----------------------------CreatePDFfromTXT
foreach my $BC (@BCDIR){
	$pid = $pm->start and next;	
	chdir "$Cdir/2_individual_barcodes/$BC"; 		
		my $PDFfile="$Cdir/2_individual_barcodes/$BC/05_report/$BC"."_report.pdf";
   	unless(-e $PDFfile){
					if($Verbose){say "  =>$BC: running R14b_CreatePDFfromTXT..."; }
					say MAINLOG "  =>$BC: R12_CreatePDFfromTXT...";
					$Time=localtime; say MAINLOG "$BC: started: ...\t $Time";
					LORCAN_subs::R11_CreatePDFfromTXT({
								Verbose => $Verbose,
								IN 		=> "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report",
								OUT 	=> $PDFfile,
					});
					$Time=localtime; say MAINLOG "$BC: finished: ...\t $Time";
					unless (-e $PDFfile)
						{say "$BC report.pdf (R12) not produced. Exiting..."; exit;}
	}else {if($Verbose){say "  =>$BC: (R12_CreatePDFfromTXT) already completed"; }}	
$pm->finish; # do the exit in the child process	
}# end of R12
$pm->wait_all_children;
#------------------
foreach my $BC (@BCnoCons){
	$pid = $pm->start and next;	
	chdir "$Cdir/2_individual_barcodes/$BC"; 		
		my $PDFfile="$Cdir/2_individual_barcodes/$BC/05_report/$BC"."_report.pdf";
   	unless(-e $PDFfile){
					if($Verbose){say "  =>$BC: running R14b_CreatePDFfromTXT..."; }
					say MAINLOG "  =>$BC: R12_CreatePDFfromTXT...";
					$Time=localtime; say MAINLOG "$BC: started: ...\t $Time";
					LORCAN_subs::R11_CreatePDFfromTXT({
								Verbose => $Verbose,
								IN 		=> "$Cdir/2_individual_barcodes/$BC/05_report/$BC.report",
								OUT 	=> $PDFfile,
					});
					$Time=localtime; say MAINLOG "$BC: finished: ...\t $Time";
					unless (-e $PDFfile)
						{say "$BC report.pdf (R12) not produced. Exiting..."; exit;}
	}else {if($Verbose){say "  =>$BC: (R12_CreatePDFfromTXT) already completed"; }}	
$pm->finish; # do the exit in the child process	
}# end of R12
$pm->wait_all_children;

#---------------------------- R13_file cleaning
unless ($args{V}) {
	foreach my $BC (@BCDIR){
		qx(rm -f $Cdir/2_individual_barcodes/$BC/*.fasta);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/*.fai);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/log.txt);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/new.*);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/aln.*);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/list.mapped);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/list.mappedS);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/02_individual_consensus/*/aln1.*);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/02_individual_consensus/*/list_*);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/02_individual_consensus/*/new*);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/02_individual_consensus/*/R7g.Rmd);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/AllBlastREFs.fasta);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/AllBlastREFsUniq.fasta);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/OUTREFs.fasta);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/CONSREF1_simple.fasta);
	}
	foreach my $BC  (@BCnoCons){
		qx(rm -f $Cdir/2_individual_barcodes/$BC/*.fasta);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/*.fai);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/log.txt);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/new.*);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/aln.*);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/list.mapped);
		qx(rm -f $Cdir/2_individual_barcodes/$BC/01_mapping_reads_to_refDB/list.mappedS);
	}
} else {say "R13: no intermediary file removed (verbose mode on)"}
#-----------------------------------------------------

#---------------------------- R14_collecting PDF
unless(-e "$Cdir/3_PDF_reports")  {system("mkdir $Cdir/3_PDF_reports");}	
foreach my $BC (@BCDIR){
	system("cp $Cdir/2_individual_barcodes/$BC/05_report/$BC\_report.pdf $Cdir/3_PDF_reports/");
}
foreach my $BC (@BCnoCons){
	system("cp $Cdir/2_individual_barcodes/$BC/05_report/$BC\_report.pdf $Cdir/3_PDF_reports/");
}
#-----------------------------------------------------


if ($Verbose) {	say "Finished LORCAN workflow on: ",scalar(localtime);}
say MAINLOG "=> Finished LORCAN workflow on:",scalar(localtime);
#########################################
		
		
		
		
#Emailing alert switch off in verbose mode ------------------
unless ($args{V}) {
			if($EMAIL){ say MAINLOG "Email was sent to: $EMAIL";
				qx(	mailx -s \"LORCAN analyses complete\" $EMAIL <<< \"ANALYSES CAN BE FOUND IN DIRECTORY: $Cdir\");
				}
}

