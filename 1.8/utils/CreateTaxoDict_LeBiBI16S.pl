#!/bin/perl
#===============================================================================
#
#         FILE: CreateTaxoDict.pl
#
#
#  DESCRIPTION: Create a taxdict.csv file from a FASTA reference file (input needed for LORCAN)
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#       AUTHOR: Alban Ramette, alban.ramette@ifik.unibe.ch
# ORGANIZATION: IFIK
my      $version="1.0";
#      CREATED: 01/08/2019 
#         BUGS: ---
#         TODO: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use feature qw(say);
#-------------
use Cwd;
my $cwdir = getcwd;
#print "cwd: $cwdir\n"; # where the script was run
#-------------- flags
use Getopt::Std;

my %args;
getopts('i:o:vVh', \%args);


my ($Verbose,$REF_FP,$OUTFILE);

if ($args{h}) {
        print "\n\nWhat it does: parse the referbce fasta headers to format the csv file as required for LORCAN \n";
        print "Arguments:...\n";
        print "\t-i Full path to the reference fasta \n";
        print "\t-o Full path to the output csv file  \n";
		
        print "\t-h help information....\n";
        print "\t-V verbose mode (for debugging)\n";
        print "\t-v version number\n";
        print "\nUsage example: \n";    
         say "\nperl createTaxoDict.pl -i /storage/databases/leBIBI/SSU-rDNA/002/SSU-rDNA-mk37_stringent/16S_stringent_dedup.fasta -o LeBiBiRef.taxdict.csv\n";
}else{ 
	if ($args{V}) {$Verbose=$args{V}; chomp $Verbose;}
	
	if ($args{i}) { #
       if ($Verbose) {print "-i has a value of: $args{i} (reference FASTA file to parse)\n";}
		$REF_FP=$args{i}; chomp $REF_FP;
	}
	if ($args{o}) { #e.g. "/path2/the/out.csv"
			if ($Verbose) {print "-o has a value of: $args{o} (output .csv file)\n";}
			$OUTFILE=$args{o}; chomp $OUTFILE;
	}
	if ($args{v}) { 
			say $version;exit;
	}

}


######################################################################
#### MAIN ############################################################
if ($REF_FP && $OUTFILE){
	if ($Verbose) {say "OK: all arguments are present";}
} else {say "Some input parameters are missing. See -h option for some help. Exiting.";exit;}

open OUT, ">$OUTFILE";

my $FASTAHEADERS=qx(grep  ">"  $REF_FP | cut -f2 --delim=">" );
my @Headers=split "\n", $FASTAHEADERS; chomp @Headers;
my (@Taxo,@AN);
#Abiotrophia_defectiva~v~TT~URS00000B1AF5=Bacteria-Firmicutes-Bacilli-Lactobacillales-Aerococcaceae-Abiotrophia-Abiotrophia_defectiva
foreach my $FullLine (@Headers){
	$FullLine=~ m/^(.+)\~.+\~.+\~(.+)\=/;
	 if($1){ print OUT "$2,";}else {print OUT "NA,";}
	 if($2){ print OUT "$1,";}else {print OUT "NA,";}
	 say OUT $FullLine;

}
