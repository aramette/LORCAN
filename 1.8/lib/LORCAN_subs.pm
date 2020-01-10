package LORCAN_subs;

use 5.013;

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

use strict; 
#---------------------------------------------
$VERSION     = 1.08;
@ISA         = qw(Exporter);
@EXPORT      =  qw(			
					$R0_test_dependencies
					$R1_Seqkit
					$R2_porechop
					$R4_parseMinNberReads
					$R5_Extract_modal_sequences
					$R6_MapWithMinimap2
					$R7_calc_stats_perContig
					$R8_Consensus_format_extract
					$R9_check_seq_blast_taxo
					$R9_2_check_seq_phylogeny
					$R10_createTxtReport
					$R11_createTxtReportNo
					$R11_CreatePDFfromTXT
);   		

#---------------------------------------------
use FindBin;  # locate this script
use lib "$FindBin::Bin/../lib";
use config;
use feature qw(say);

# SUBS ---------------------------------------
sub R0_test_dependencies{
	my ($args)=@_;
	my $Verbose			=$args->{Verbose};
	my $LIBDIRECTORY	=$args->{LIBDIRECTORY};
	my $PERL			=$args->{PERL};
	my $SEQKIT			=$args->{SEQKIT};
	my $PORECHOP		=$args->{PORECHOP};
	my $minimap2		=$args->{minimap2};
	my $samtools		=$args->{samtools};
	my $BLASTN			=$args->{BLASTN};
	my $Rscript			=$args->{Rscript};
	my $REF_FP			=$args->{REF_FP};# $RefDB{$DBchoice}; # %RefDB defined in config.pm file
	#-------------------
	
	
  unless(-e "$LIBDIRECTORY"){
	say "!!! Issue with LIBDIRECTORY. Exiting"; exit;
  }
  my $PERL_test=qx($PERL -e "print 'OK'"); 
  unless($PERL_test eq "OK"){
		say "!!! Could not find Perl binary. Exiting"; exit;
  }
  my $SEQKIT_test=qx($SEQKIT version); chomp $SEQKIT_test;
  if($SEQKIT_test =~/"not found"/){
		say "!!! Could not find Seqkit binary. Exiting"; exit;
  }
  my $PORECHOP_test=qx($PORECHOP --version); chomp $PORECHOP_test;
  if($PORECHOP_test =~/"not found"/){
		say "!!! Could not find porechop binary. Exiting"; exit;
  }
  my $minimap2_test=qx($minimap2 --version); chomp $minimap2_test;
  if($minimap2_test =~/"not found"/){
		say "!!! Could not find minimap2 binary. Exiting"; exit;
  }
  my $samtools_test=qx($samtools --version); chomp $samtools_test;
  if($samtools_test =~/"not found"/){
		say "!!! Could not find samtools binary. Exiting"; exit;
  }
  my $BLASTN_test=qx($BLASTN -version); chomp $BLASTN_test;
  if($BLASTN_test =~/"not found"/){
		say "!!! Could not find BLASTN binary. Exiting"; exit;
  }
  my $Rscript_test=qx($Rscript -e 1); chomp $Rscript_test;
  unless($Rscript_test =~/\[1\]/){
		say "!!! Could not find Rscript binary. Exiting"; exit;
  }
  #checking db
  	#my $REF_FP= $RefDB{$DBchoice}; # %RefDB defined in config.pm file
	unless(-e $REF_FP){say "Database file cannot be found (Test db) Exiting R0!";exit;}
	my @Fparse=File::Basename::fileparse($REF_FP);
	my $blastDBFile= $Fparse[0]; #say "Fasta reference: $blastDBFile";
	my $blastDBFolder=$Fparse[1]; #say "Fasta reference folder: $blastDBFolder";
	my $TotalNberReadsinREF=qx( grep -c ">" $REF_FP); chomp $TotalNberReadsinREF;
	
	if($Verbose) {say "Checking db\n\tblastDBFile: $blastDBFile";say "\tblastDBFolder: $blastDBFolder";}	
	if($Verbose){say "Test db is ok (file found). Continuing..."}

  
  #sub checking for Perl module dependencies
	my @mod=("FindBin","Cwd","Parallel::ForkManager","PDF::API2","Tk::PNG");
	sub CheckModuleNotInstalled{
		my @array = @_;my $ret;
		if($Verbose) {say "Checking for Perl module dependencies:";}
		my @Bad;
		foreach (@array){
			$ret= qx($PERL -M$_ -e \";\" 2>&1);
			if($ret =~ /Can't locate/){
							if($Verbose) {say "not installed!!! ... $_";};
							push(@Bad,$_);
			}else {
							if($Verbose) {say "OK               ... $_";}
			}
		}
		return (@Bad);
	}
	my @NotInstalled=CheckModuleNotInstalled(@mod);
	if (scalar (@NotInstalled) > 0){say "Some Perl modules are not installed: ";
		#foreach (@NotInstalled) {say $_;}
		say "Exiting"; exit;
	} 
	
	say "OK ... All tests for script dependencies sucessfully done.";
	
}
#--------------------------------
sub R1_Seqkit{
	my ($args)=@_;
	my $Verbose=$args->{Verbose};
	my $IN=$args->{IN}; #directory
	my $OUT=$args->{OUT};
	my $NberThreads=$args->{NberThreads};
	my $SEQKIT=$args->{SEQKIT};
	unless($IN=~/.*\/$/){$IN=$IN."/";}
		system("ls $IN* \| parallel --will-cite -j $NberThreads -P $NberThreads \" $SEQKIT stats \{\}\" >> $OUT");
}
#--------------------------------
sub R2_porechop{#porechop demultiplexing		
	my ($args)=@_;
	my $PORECHOP=$args->{PORECHOP};
	my $IN=$args->{IN};
	my $NberThreads=$args->{NberThreads};
	my $OUT=$args->{OUT}; #$Cdir/1_fasta/
	system("$PORECHOP -v 0 -i $IN -t $NberThreads -b $OUT --format fasta --discard_unassigned --require_two_barcodes > /dev/null 2>&1"); 
}
#--------------------------------
sub R4_parseMinNberReads{
	my ($args)=@_;
	my $Cdir					=$args->{Cdir};
	my $IN						= $args->{IN};	#"$Cdir/0_logs/2_stats_fastq_postporechop.txt"
	my $MinRead_threshold		=$args->{MinRead_threshold};
	my $MainLog					=$args->{OUT}; #"$Cdir/0_logs/$LOG"
	open MAINLOG,">>$MainLog"; 
 		# parse the stats in 2_stats_fastq_postporechop.txt 
		#and create directory for barcodes above read threshold
		my (@L,$line,$BCname,$Nreads,@BCnamesOK);
		open STAT,"$IN" or die "cannot open the file 2_stats_fastq_postporechop.txt";
		foreach $line (<STAT>){ #screening for the right nber of reads per Barcode
			unless($line =~/file/){
				chomp $line;
				@L=split(/\s+/,$line);chomp @L;
				$L[0] =~ /\/.+\/(.+)\.fasta/; $BCname=$1; #getting the barcode name 		 #say $BCname;
				$Nreads=$L[3]; # nber of reads in the fasta
				$Nreads=~s/\,//g;
				if ($Nreads >= $MinRead_threshold){push(@BCnamesOK, $BCname);} 
				else {print MAINLOG "\t $BCname: $Nreads < threshold\n";}
			}
		}
		print MAINLOG "(reporting only problematic barcodes)\n";
		foreach my $BC (@BCnamesOK){system("mkdir -p $Cdir/2_individual_barcodes/$BC")}#create barcode directories accordingly
			print MAINLOG "\n=> Created barcode directories accordingly\n";
				foreach my $BC (@BCnamesOK){print MAINLOG "$BC\t"; }#
				print MAINLOG "\n";
	return(@BCnamesOK);
}
sub R5_Extract_modal_sequences{
	my ($args)=@_;
	my $Verbose=$args->{Verbose};
	my $BCDIRref=$args->{IN}; #array reference
	my $MainLog=$args->{OUT}; #"$Cdir/0_logs/$LOG"
	my $Cdir=$args->{Cdir};
	my $NberThreads=$args->{NberThreads};
	my $Rscript=$args->{Rscript};
	our $delta=$args->{delta}; # +- delta bp from the mode
	our $MaxRead=$args->{MaxRead}; # Maximum number of reads to retain per sample
	my $LIBDIRECTORY=$args->{LIBDIRECTORY};
	#----------
	my @BCDIR=@$BCDIRref; #deferencing $ref
	my $NBc=scalar(@BCDIR);
	if($NberThreads gt $NBc ){$NberThreads=$NBc}
	my $pm = new Parallel::ForkManager($NberThreads);
	my $pid;
	foreach my $BC (@BCDIR){  # produce the mode+-delta bp file
					if($Verbose) {say "$BC being processed";}
					print MAINLOG "$BC\t";
				$pid = $pm->start and next; 
					my $Fas="$Cdir/1_fasta/$BC.fasta";
					my $SCMstat="$Cdir/1_fasta/$BC\_stats.txt";
					
					#....
#caution: the text below MUST _NOT_ be indented until END_RMDtxt.
my $RTEXT = <<END_Rtxt;					
FASTAfile="FILENAME"
msg.trap <- capture.output( suppressMessages( library(seqinr)))
cat("============\n",file="OUTF",append=FALSE)
cat(paste0("Sample name:\t",FASTAfile,"\n"),file="OUTF",append=TRUE)
FASTAfile=gsub(pattern = "/",replacement = "//",FASTAfile)
fa <- read.fasta(FASTAfile)# Loading input file
Nseq=length(fa)
cat("============\n",file="OUTF",append=TRUE)
cat(paste0("Total sequences after porechop: ",Nseq,"\n"),file="OUTF",append=TRUE)
Seqs.length<-sapply(fa,length,simplify="array")
Seqs.length.max <- max(Seqs.length)
Seqs.length.mean <- mean(Seqs.length)
Seqs.length.median <- median(Seqs.length)
X=Seqs.length
cat("Summary of sequence lengths:\n",file="OUTF",append=TRUE)
sink(file="OUTF",append=TRUE)
	print(summary(Seqs.length))
sink(NULL)
Mode=as.numeric(names(sort(-table(X)))[1])
Xmode=X[X<=Mode+100 & X>=Mode-DELTA]
Xmode.n=length(Xmode)
cat(paste0("\nNumber of sequences close to the mode (SCM): ",Xmode.n,"/",length(fa)," (",round(Xmode.n*100/length(fa),1),"%)\n"),file="OUTF",append=TRUE)
cat("Summary of SCM lengths:\n",file="OUTF",append=TRUE)
sink(file="OUTF",append=TRUE)
	print(summary(Xmode))
sink(NULL)

if(Xmode.n >= MAXN){Xmode=Xmode[1:MAXN]
					cat(paste0("\nKeeping only first ",MAXN," SCMs\n"),file="OUTF",append=TRUE)
					cat("Summary of ",MAXN," modal sequence lengths:\n",file="OUTF",append=TRUE)
					sink(file="OUTF",append=TRUE)
						print(summary(Xmode))
					sink(NULL)
} #limit only the top number of Xmode reads 

# Extracting sequences closest to the mode
write.fasta(sequences=fa[names(Xmode)], 
			names=names(Xmode),
			file.out=paste0(FASTAfile,"_mode_closest.fasta"),
)
END_Rtxt
			$RTEXT=~s/FILENAME/$Fas/;  # edit what needs to be changed
			$RTEXT=~s/OUTF/$SCMstat/g;  # edit what needs to be changed
			$RTEXT=~s/DELTA/$delta/g;  # edit what needs to be changed
			$RTEXT=~s/MAXN/$MaxRead/g;  # edit what needs to be changed
			my $Rfile="$Cdir/1_fasta/".$BC."_Rscript.R";
						open OUTR,">$Rfile";
						say OUTR $RTEXT;
							system("$Rscript $Rfile >> $MainLog");
							system("rm -f $Rfile");
					system("echo $BC>> $Cdir/0_logs/3_Barcodes_with_modal_sequences.txt");
				$pm->finish; # do the exit in the child process
	}
	say "";
	$pm->wait_all_children;
}

sub R6_MapWithMinimap2{
	my ($args)=@_;
	my $Verbose=$args->{Verbose};
	my $IN=$args->{IN};# 
	my $OUT1=$args->{OUT1};# "new.prim.bam",
	my $OUT2=$args->{OUT2};# "list.mapped",
	my $Wdir=$args->{Wdir};# working directory
	my $REF_FP=$args->{REF_FP};# $REF_FP= $RefDB{$DBchoice} # %RefDB defined in config.pm file
	my $minimap2=$args->{minimap2};
	my $NberThreads=$args->{NberThreads};	
	my $samtools=$args->{samtools};
	my $PERL=$args->{PERL};

			if ($Verbose){say "   Minimap mapping => aln.sam, aln.sorted.bam produced...";}
				my @Fparse=File::Basename::fileparse($REF_FP);
				my $blastDBFile= $Fparse[0]; #say "Fasta reference: $blastDBFile";
			system("cp $REF_FP $Wdir");
				chdir $Wdir;
			system("$minimap2 -t $NberThreads -ax map-ont $blastDBFile $IN > aln.sam");
			system("$samtools faidx $blastDBFile;$samtools import  $blastDBFile.fai  aln.sam  aln.bam;$samtools sort  aln.bam  -o aln.sorted.bam");
		#filtering the secondary alignments (removing multi mappers)
		system("$samtools view -F 256 -b aln.sorted.bam > aln.prim.bam; $samtools index aln.prim.bam");
		system("$samtools  idxstats aln.prim.bam > aln.prim.bam.stats");
			open STATS,"aln.prim.bam.stats";
			open LM, ">list.mapped"; 
			open LMS, ">list.mappedS"; 
			my @STATSL;
			foreach(<STATS>){
					@STATSL=split(/\t/,$_);
					if(int($STATSL[2]) > 0){say LM "$STATSL[2] $STATSL[0]";}
					if(int($STATSL[2]) > 0){say LMS $STATSL[0];}
			
			}
			close STATS;close LM;close LMS;
			#need only the REF for grepping
			system("$samtools view -h aln.prim.bam > aln.prim.sam");
			system("grep \"\@HD\" aln.prim.sam > new.sam");
			system("grep \"\@PG\" aln.prim.sam >> new.sam");
			system("grep -f list.mappedS aln.prim.sam >> new.sam");
			system("$samtools view -b -S  new.sam >  $OUT1 ");
			system("$samtools index $OUT1");
	
}
sub Matchtaxo{
	#reformat list.mapped to a format that is easy to be parsed by R7
	my ($args)=@_;
	my $Verbose=$args->{Verbose};
	my $IN1=$args->{IN1};	#"list.mapped",	 
	my $IN2=$args->{IN2};	#"LeBiBiRef.taxdict.csv",	 
	my $OUT=$args->{OUT}; 	#list.mapped.taxo
	#------------
	if ($IN1 && $IN2 &&  $OUT){
		unless(-e $IN1){say "-i $IN1 does not exist; exiting";exit;}
		unless(-e $IN2){say "-j $IN2 does not exist; exiting";exit;}
	} else {say "missing arguments! exiting"; exit;}
	open LM,"<$IN1";
	my @LM_data= <LM>; chomp @LM_data; 
	@LM_data = grep(!/^$/, @LM_data); # removing empty lines if present
	open TD,"<$IN2";
	my @TD_data= <TD>; chomp @TD_data; 
	open OUT,">$OUT";
	my ($AN,@Info,@Infolm);
	foreach my $td (@TD_data){
		@Info=split ",", $td;
		$AN= shift @Info;
		foreach my $lm (@LM_data){
			 if($lm =~/$AN/) {
				 @Infolm= split "\t", $lm;
				 say OUT  "$Infolm[0],$AN,$Info[0],$Infolm[1]"; 
			 }
		}	
	}
	close OUT;
}

sub ParseMpileupOutput{
	use feature qw(say);
	my $INFILE=shift;
	my $OUTFILE=shift;
	open IN,"<$INFILE" or die "cannot open infile";
	open OUT,">$OUTFILE";

	say OUT "Pos\tA\tG\tC\tT\tdel\tins\tmxins\tmxinseq\tinserted\tambiguous";
		# [Pos] base position
		# [del] number of deletions
		# [ins] number of insertions
		# [mxins] highest number of insertions for the elemetn specified in tMxInsSeq
		# [mxinseq] sequence of the most frequent insertion 
		# [inserted] all insertions sequences, comma separated
		# [ambiguous] all tambiguous sequences, comma separated
	my ($bp,$bases,$ref);
	foreach my $line (<IN>){
        my @data = split(/\t/,$line);
         $bp = $data[1];
         $bases = uc $data[4];
			chomp $bases;
         $ref = uc $data[2];
		 # say "bp: $bp";
		 # say "bases: $bases";
		 # say "ref: $ref";
		 our %types = ("A",0,
					  "G",0,
					  "C",0,
					  "T",0,
					  "-",0,
					  "+","",
					  "X","");
		my @KEYS=keys %types;  #foreach (@KEYS) {say "keys: $_";}
		my @BASES=split(//,$bases);
		my $base;
		my $i=0;		
		while( $i<@BASES){
			$base=$BASES[$i];
			if ($base =~ /\^/ || $base =~ /\$/ ){$i++;}
			elsif($base =~ /\-/){$i++;}
			elsif($base =~ /\*/){$types{"-"} += 1;}
			elsif($base =~ /\+/){
						$i++;
                        my $addNum = int($BASES[$i]);
                        my $addSeq = "";
                        for my $a (1..$addNum){
                                $i++;
                                $addSeq = $addSeq.$BASES[$i];
						}
						if($types{"+"} ne ""){$types{"+"}=$types{"+"}.",".$addSeq;}else{$types{"+"}=$addSeq;}
			}
			elsif($base =~ /\./ || $base =~ /\,/){
			             $types{$ref}++;
			} else{
				if (grep {$_ eq $base} @KEYS) {
					$types{$base}++;
				 } else{
					if($types{"X"} ne ""){
						$types{"X"}=$types{"X"}.",".$base;
					}else{$types{"X"}=$base;
					}
				}
			}			
            $i++;
		}
		# say "A:", $types{"A"};
		# say "T:", $types{"T"};
		# say "C:", $types{"C"};
		# say "G:", $types{"G"};
		# say "-:", $types{"-"};
		# say "+:", $types{"+"};
		# say "X:", $types{"X"};
		my $NberInsertion=	split(/\,/,$types{"+"} );
		if($types{"+"} eq ""){$types{"+"}=".";}
		if($types{"X"} eq ""){$types{"X"}=".";}
		# working on insertion sequences to determine the highest frequency of insertions
		my $MaxInsertSeq="";
		my $MaxInsertNb=0;
		if($types{"+"} ne "."){
			my @INSERT=split(/\,/,$types{"+"});
			my @uniqueINSERT = do { my %seen; grep { !$seen{$_}++ } @INSERT };
			my %INSERThash;
			for my $i (@uniqueINSERT){
				$INSERThash{$i}=grep{/^${i}$/} @INSERT;
			}
			for (keys %INSERThash){
				if($INSERThash{$_} >= $MaxInsertNb){
					$MaxInsertNb=$INSERThash{$_};
					$MaxInsertSeq=$_;
				}
			
			};
		}
		my @out = ($bp,$types{"A"},$types{"G"},$types{"C"},$types{"T"},$types{"-"},$NberInsertion,$MaxInsertNb,$MaxInsertSeq,$types{"+"},$types{"X"});
		say OUT join( "\t", @out );
	}
	close IN;
	close OUT;
#
# Pos	A	G	C	T	del	ins	mxins	mxinsseq	inserted	ambiguous
# 11	0	2	81	1	4	7	2	G	G,AG,A,G,CA,T,CTG	.
# 12	22	4	63	6	1	3	2	TT	TT,TT,T	.
# 13	1	5	18	74	5	1	1	TG	TG	.
# 14	9	82	4	2	6	2	1	AAA	AAA,AAACC	.
# 15	1	96	5	0	6	6	3	C	TT,C,AA,CT,C,C	.
#ParseMpileupOutput("test.mpileup","perl.output");
}
sub R7_calc_stats_perContig{
	my ($args)=@_;
	my $Verbose=$args->{Verbose};
	my $IN0=$args->{IN0};#"$Cdir/2_individual_barcodes/$BC",
	my $IN1=$args->{IN1};#"list.mapped.taxo",	 
	my $IN2=$args->{IN2};#"new.prim.bam",	indexed 
	my $OUT1=$args->{OUT1}; #log.txt
	my $OUT2=$args->{OUT2}; #aln1.mpileup.parsed
	my $BC=$args->{BC}; #BC28
	my $ReadThreshold=$args->{ReadThreshold};	
	my $minimap2=$args->{minimap2};
	my $samtools=$args->{samtools};
	my $REF_FP=$args->{REF_FP};# $REF_FP= $RefDB{$DBchoice} # %RefDB defined in config.pm file
	my $SEQKIT=$args->{SEQKIT};
	#--------------------
	my @Fparse=File::Basename::fileparse($REF_FP);
	my $blastDBFile= $Fparse[0];
	open OUT1, ">$OUT1";

	# foreach line in file get and parse the data
	open IN,"<$IN1" or die "cannot open file: $IN1 in R7";
	my @D=<IN>; chomp @D;
	@D = grep(!/^\n$/, @D);
	@D = grep /\S/, @D;
	@D = grep { $_ ne '' } @D;	#remove empty elements
	# split
	my (@Counts,@AN,@Taxo,@TaxoUniq,@split1,@split2,@FHeader);
	foreach (@D){
		@split1 = split /\s/, $_;
			push (@Counts,int($split1[0]));
			@split2 = split /\,/, $split1[1];
				push(@FHeader,$split2[0]);
				push (@AN,$split2[1]);
				push (@Taxo,$split2[2]);

	}
	#total counts
	my $TotalCounts=0;
	foreach (@Counts){$TotalCounts+=$_;}
	my $TotalCountsThreshold=0;
	my $roundedPercent="???";
	foreach (@Counts){$TotalCountsThreshold+=$_;}
	if($TotalCounts!=0){my $Percent=$TotalCountsThreshold*100/$TotalCounts;$roundedPercent = sprintf("%.1f", $Percent);} 
	
	say OUT1 "== A) Read counts ====================================================";
	say OUT1 "Total: $TotalCounts reads aligned to references,";
	say OUT1 "$roundedPercent% ($TotalCountsThreshold reads) were kept after applying \n  \>$ReadThreshold cutoff reads mapping per taxonomic level";
	
	if (scalar @Counts ==0 ){say OUT1 "  No more reference sequence left for $BC (after step R7).";} else {
	#determine unique taxo
	@TaxoUniq = do { my %seen; grep { !$seen{$_}++ } @Taxo };
	# say OUT1 "\n== B) Species names identified in the sample:";
	# foreach (@TaxoUniq){say OUT1 "\t-$_";}
	
	say OUT1 "\n== B) Selection of taxonomic groups based on mapped reads ============";
	#foreach unique taxo, determine the AN with the highest count; and create taxo group
		my $N=scalar @AN; #say $N;
		my (@TaxoGroup,@TaxoGroupRN);#RN= Read number
		my (@SelectedTaxoGroup,@SelectedTaxoGroupRN);#RN= Read number
		
		foreach my $UT (@TaxoUniq){ # within a taxonomic level
			my $TaxoGroupTotRN=0;
			#say "=== $UT === ";
			for (my $i =0; $i<$N;$i++){
				if($Taxo[$i] =~/^$UT$/){
				#say "$Counts[$i] $FHeader[$i]";
				$TaxoGroupTotRN+=$Counts[$i];
				push (@TaxoGroup,$FHeader[$i]);
				push (@TaxoGroupRN,$Counts[$i]);
				}
			}
			#say "TOT=$TaxoGroupTotRN";
			if($TaxoGroupTotRN>=$ReadThreshold){
				push(@SelectedTaxoGroup,$UT);
				push(@SelectedTaxoGroupRN,$TaxoGroupTotRN);
			}
			@TaxoGroup=();@TaxoGroupRN=();
		} #foreach $UT
				if(scalar(@SelectedTaxoGroup) == 0){ #not enough reads mapped!
			say OUT1 "!!!CAUTION: Not enough reads mapping to any reference for $BC (step R7)."; #"Not enough reads" is used in R8
			say OUT1 "(No taxonomic group obtained in total >$ReadThreshold reads)\n";
			say OUT1 "List of mapped reads to references:";
			foreach my $UT (@TaxoUniq){ # within a taxonomic level
				say OUT1 "=== $UT === ";
				for (my $i =0; $i<$N;$i++){
					if($Taxo[$i] =~/^$UT$/){
						say OUT1 "$Counts[$i] $FHeader[$i]";
					}
				}
			} #foreach $UT
		}else{
		# for printing in decreasing order of reads
				my %hash; # for printing in decreasing order of reads
				@hash{@SelectedTaxoGroup} = @SelectedTaxoGroupRN;#hash of unique with their nber occurrence
			foreach my $k (sort {$hash{$b} <=> $hash{$a}} keys %hash) {
				my $Percent=int($hash{$k})*100/$TotalCounts;$roundedPercent = sprintf("%.1f", $Percent);
				say OUT1 "\t-($hash{$k},$roundedPercent%)  $k"; 
				#writing the header of ref fasta per group to a directory
				my $DIRcons="$IN0/02_individual_consensus/$k"."_consensus";
				unless(-e $DIRcons) {system("mkdir -p $DIRcons");}
						open GSOUT,">$DIRcons/list\_$k"; 
									my ($Index);
									my $SelectedNber=0;
						my $SelREF;	# to produce a space separated ref to extract reads from the bam file
						for (my $i =0; $i<$N;$i++){
							if($Taxo[$i] =~/^$k$/){
								say GSOUT "$Counts[$i]\t$FHeader[$i]";
								$SelREF=$SelREF." $FHeader[$i]";
								#find the highest nber
								if($Counts[$i]>=$SelectedNber){
										$SelectedNber=$Counts[$i];$Index=$i;
								}
							}
						}
						#say "highest: $Counts[$Index] $FHeader[$Index]";
						#say "$SelREF";
						
						system("$samtools view  -b $IN2 $SelREF  > $DIRcons/new1.bam");  # select for the reads
						system("$samtools fasta $DIRcons/new1.bam  > $DIRcons/new1.fasta 2>&1"); # extract the reads 

						#extract ref
						system("cat $blastDBFile | $SEQKIT grep -r -p $FHeader[$Index] > $DIRcons/REF1.fasta");
						system("$minimap2 -ax map-ont $DIRcons/REF1.fasta $DIRcons/new1.fasta > $DIRcons/aln1.sam ");
						system("$samtools faidx $DIRcons/REF1.fasta;$samtools import $DIRcons/REF1.fasta.fai  $DIRcons/aln1.sam  $DIRcons/aln1.bam;$samtools sort  $DIRcons/aln1.bam  -o $DIRcons/aln1.sorted.bam");
						#now calculate consensus from the bam file
						system("$samtools mpileup $DIRcons/REF1.fasta;$samtools import $DIRcons/REF1.fasta.fai  $DIRcons/aln1.sam  $DIRcons/aln1.bam;$samtools sort  $DIRcons/aln1.bam  -o $DIRcons/aln1.sorted.bam");
						system("$samtools mpileup -f $DIRcons/REF1.fasta $DIRcons/aln1.sorted.bam > $DIRcons/aln1.mpileup");
						#parse pileup
						unless(-e "$DIRcons/$OUT2"){
							if($Verbose){say "  =>$BC: running ParseMpileupOutput on $DIRcons..."; }
							ParseMpileupOutput("$DIRcons/aln1.mpileup","$DIRcons/$OUT2");
							
						
						
						} else {
							if($Verbose){say "  =>$BC: (ParseMpileupOutput on $DIRcons) already completed"; }
						}
						
						R7_calculateConsensusPerGroup({
										Verbose 		=> $Verbose,
										IN1 			=> "$DIRcons/aln1.mpileup.parsed",
										BC 				=> $BC,
										OUT1 			=> "$DIRcons/",
										OUT2 			=> "ConsensusPlot.html",
						});
			
			}#foreach my $k
		} #if (scalar @Counts !=0 
	
	}
	
}	

sub R7_calculateConsensusPerGroup{
		my ($args)=@_;
		my $Verbose=$args->{Verbose};
		my $IN1=$args->{IN1};#$DIRcons/aln1.mpileup.parsed
		my $BC=$args->{BC};#BC28
		my $OUT1=$args->{OUT1};#"$DIRcons/"
		my $OUT2=$args->{OUT2};#"ConsensusPlot.html",
		#-----------
		if($Verbose){say "  =>$BC:  running R7_calculateConsensusPerGroup"; }
		open ANLIST,"<$IN1";# 	
		my $ANLIST;
		foreach(<ANLIST>){$ANLIST=$_;}
		my @ALIST=split(" ",$ANLIST); chomp @ALIST; #foreach(@ALIST) {say $_;}
		if (scalar(@ALIST) == 0){say "\@ALIST array is empty (step R7). No accession number found.";}
		#(apply cutoffs for min nber of reads per pos and overall)
		
#caution: the text below MUST _NOT_ be indented until END_RMDtxt.
my $RMDTEXT = <<'END_RMDtxt';
---
title: "LORCAN output"
Author: 'IFIK bioinformatic group'
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  html_document
---

=== Consensus plots === 
```{r echo=FALSE, fig.width=15, fig.height=10}

#version=1
# Alban Ramette
msg.trap <- capture.output( suppressMessages( library(Biostrings)) )
msg.trap <- capture.output( suppressMessages( library(magrittr)) )
msg.trap <- capture.output( suppressMessages( library(dplyr)) )
#---------------------------------------
F="PATHTOFILE"; # F="aln1.mpileup.parsed"
getwd()
CWD=getwd()
Taxo=last(strsplit(CWD,"/")[[1]])

F=gsub(pattern = "/",replacement = "//",F)
DATA=read.delim(F,h=TRUE,as.is = TRUE)
colnames(DATA)[6]<-"-"
DATA$MaxAGCT<-apply(DATA[,2:5],1,max) #max per base per position

# need to integrate also insertions
NotableInsertions <- DATA %>% filter(mxins>MaxAGCT)
if(nrow (NotableInsertions) == 0) { # nothing to add
	M<-as.matrix(t(DATA[,2:6]))
	colnames(M)<- DATA$Pos
	ConsensusSeq<-consensusString(M,ambiguityMap="N",threshold=0.5)
}else { #add inserted sequences to the matrix before computing the consensus
	Di=cbind(DATA[,1:6],ins=rep(0,nrow(DATA)))
		for (R  in 1:nrow(NotableInsertions)){
			IN=NotableInsertions[R,]
			POS1=as.numeric(IN[1])
			Distart<-which(Di$Pos==POS1)-1
			Di1=Di[1:Distart,] #up part
			Di2=Di[(Distart+2):nrow(Di),] #down part
			#get individual insertion bases
			INSbases<-strsplit(x=IN$mxinseq,split="")[[1]]
			INSbasesN <-length(INSbases)
			#update positions for Mi2
			Di2$Pos<-Di2$Pos+INSbasesN-1
			#create a dataframe, and filled it
			INSdf<-data.frame(Pos=seq(from=POS1,to=POS1+(INSbasesN-1)),"A"=rep(0,INSbasesN),"G"=rep(0,INSbasesN),"C"=rep(0,INSbasesN),"T"=rep(0,INSbasesN),del=rep(0,INSbasesN),"ins"=rep(1,INSbasesN))
			colnames(INSdf)[6]="-"
			for (B in 1:INSbasesN){
				 if(INSbases[B] == "A"){INSdf[B,"A"]<-IN$mxins
				 } else if (INSbases[B] == "G"){INSdf[B,"G"]<-IN$mxins
				 } else if ( INSbases[B] == "C"){INSdf[B,"C"]<-IN$mxins
				 } else if ( INSbases[B] == "T"){INSdf[B,"T"]<-IN$mxins
				 }
			}
			#merge with Di1 and Di2
			Di=rbind(Di1,INSdf,Di2)
			rownames(Di)<-1:nrow(Di)
			#update NotableInsertions$Pos
			Nnotab=nrow(NotableInsertions)
			if(R<Nnotab){NotableInsertions$Pos[(R+1):Nnotab]<-NotableInsertions$Pos[(R+1):Nnotab]+(INSbasesN-1)}
		}#end of R  in 1:nrow(NotableInsertions)
		#test with 2 short insertions OK
			M1<- cbind(A=as.integer(Di$A),G=as.integer(Di$G),C=as.integer(Di$C),T=as.integer(Di$T),del=as.integer(Di$"-"))
			colnames(M1)[5]="-"
			M=as.matrix(t(M1))
			colnames(M)<- as.character(Di$Pos)
			ConsensusSeq<-consensusString(M,ambiguityMap="N",threshold=0.5)
}#end of else   #add inserted sequences to the matrix b


##==========================================================================================================================
#PLOTS	
	MIN_NT_DEPTH=10 #number to keep
	TotalLength= width(ConsensusSeq)
	WS=50 #WS= window size (bp) for the plots
	BINS=seq(1,TotalLength,WS)
	DP<-apply(M,2,sum) #depth per position
	MAX= max(DP)
  for(i in BINS){
			PosStart=i
			PosEnd=i+WS-1
			IMG=paste0("Consensus sequence for ",PosStart,"-",PosEnd," bp\n")
				#png(filename=IMG)
				if(i == rev(BINS)[1]){PosEnd=TotalLength} #to accomodate shorter bin for the last one.
				Data <- M[,PosStart:PosEnd]
				if(!is.null(dim(Data))){
				Npos <- ncol(Data) #nber of positions
				par(mar=c(5,4,5,5))
				plot(1:Npos,(1:Npos)/Npos,type="n",ylim=c(0,1.2),xaxt="n",yaxt="n",
					 xlab="Positions",ylab="Base frequency",las=1)
				abline(v=1:Npos,col="grey");abline(h=seq(0.1,1,0.1),col="grey");abline(h=0.5,col="grey",lwd=2) # grid
				axis(1,at=1:Npos,labels=colnames(Data),las=2,cex.axis=0.8)
				axis(2,at=seq(0,1.2,0.2),labels=c(seq(0,1.,0.2),""),las=2,cex.axis=1)
				lines(DP[PosStart:PosEnd]/MAX,col="blue",lwd=2) #coverage line
				abline(h=MIN_NT_DEPTH/MAX,col="blue",lty=2) # cutoff to retain the base
					#mtext(paste0("<- abs cutoff(",MIN_NT_DEPTH,")"), side=4, line=1,col=1)
				axis(4,at=seq(0,1,0.1),labels=round(seq(0,1,0.1) * MAX,0),col.axis="blue",las=2)
				mtext("Total sequence depth per position (N)", side=4, line=3,col="blue")
				abline(h=0.5,lwd=2)
				# add consensus seq above; color the "N" cases in red
				CONSseq <- strsplit(substr(ConsensusSeq,start = PosStart,stop=PosEnd),
									fixed = TRUE,split = "")[[1]]
				CONSseqCol <- rep(1,Npos)
				CONSseqCol[which(CONSseq=="N")] <- 2  
				CONSseqCol[which(CONSseq=="-")] <- 2  
				#indicating inserted sequences
					if(exists("Di")){
						InsPos=Di[PosStart:PosEnd,"ins"]
					} else{InsPos=rep(0,WS)}
				for(P in 1:Npos){
				  Dp <- Data[,P]
				  if(sum(Dp)!=0){
					Dp <- Dp[Dp!=0]
					Dp1<-Dp/MAX
					text(x=P,y=Dp1,labels = names(Dp),col = CONSseqCol[P],cex=1)
					if(CONSseqCol[P]==2){points(x=P,y=0.5,pch=15,col=2,cex=2)}
					if(CONSseqCol[P]==3){points(x=P,y=0.5,pch=15,col=3,cex=2)}
				  }
				 text(x=P,y=1.1,labels = CONSseq[P],col = CONSseqCol[P],cex=1.5) #consensus character
				if(InsPos[P]==1){text(x=P,y=1.15,labels = "v",cex=1.5,col="red")} #indicator of insertion
				 
				}
				#add indication that there is an inserted sequence as compared to consensus
				
			
				title(IMG,line = 2)
				#dev.off()
		}#if null Data
}
### write to single fasta and to one multi fasta
FILENAME <- paste0("BCXX_ref_",Taxo,".fasta")

ConsensusSeqNoGaps=gsub("-","",x=ConsensusSeq)
TotalLengthNoGaps= width(ConsensusSeqNoGaps)
ConSPLIT=strsplit(ConsensusSeq,"")[[1]]
NberN=length(which(ConSPLIT=="N"))
Nberdel=length(which(ConSPLIT=="-"))

cat(paste0(">BCXX_using_",Taxo,"_",TotalLengthNoGaps,"bp_N",NberN,"_D",Nberdel,"_",format(Sys.time(), "%b_%d_%Y"),"\n"),file = FILENAME,append = FALSE)

cat(ConsensusSeqNoGaps,file = FILENAME,append = TRUE)


```
END_RMDtxt

			$RMDTEXT=~s/PATHTOFILE/$IN1/;  # edit what needs to be changed
			$RMDTEXT=~s/BCXX/$BC/g;  	   # edit what needs to be changed
			my $NewRmd="$OUT1/R7g.Rmd";
			#say "NewRmd: $NewRmd";

			unless (-e $OUT1){mkdir $OUT1;}

			open OUT,">$NewRmd";
			say OUT $RMDTEXT;
			
			unless(-e "$OUT1/No_consensus_produced.txt"){
					system("$Rscript -e \"rmarkdown\:\:render\(\'$NewRmd\',output_file=\'$OUT2\',quiet=TRUE\)\"");
			}
			
}
sub R8_Consensus_format_extract{
	my ($args)=@_;
	my $Verbose=$args->{Verbose};
	my $BCDIRref=$args->{IN}; #my $BCDIRref = \@BCDIR;
	my $OUT1=$args->{OUT1}; #log.txt",
	my $OUT2=$args->{OUT2}; #consensus_sequences.fasta",
	my $Cdir=$args->{Cdir}; #$Cdir
	#--------------
	my @BCDIR=@$BCDIRref;
	foreach my $BC (@BCDIR){ #getting the consensus sequences, reformating, extracting info from headers => log.txt of the barcode
		my $Allconseqfilename="$BC\_$OUT2";
	unless(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename"){
		open INLOG,"<$Cdir/2_individual_barcodes/$BC/$OUT1";
		 my @TaxoOrder;#getting the same order of taxo levels as in the log
		 foreach my $line(<INLOG>){
			 if($line=~/\-\(.+\%\)(.+)/){	#say $1;
				 push(@TaxoOrder,$1);
			 } 
		 }
		close INLOG;
		s/^\s+// foreach @TaxoOrder;
		chomp @TaxoOrder;
		#say foreach (@TaxoOrder);
		open OUTLOG,">>$Cdir/2_individual_barcodes/$BC/$OUT1";
		if(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/"){
		my $Cons=qx(ls $Cdir/2_individual_barcodes/$BC/02_individual_consensus/*/$BC\_*consensus.fasta);
		my @CONSFP=split(/\n/,$Cons);
		my %HEADERS;
		my %CONTENTSsamekeys;
		foreach my $F(@CONSFP){
				 open CONS,"<$F";
				 my $header;
				 my $headerTaxo;#to extract only the taxo
				 my $Seq;
				 while(<CONS>){
					 if(/^>/) {$header=$_;
					 }else {$Seq=$_;}
				}
				close CONS;
				chomp $header;
				$header=~/^\>BC.+\_using\_(.+)\_consensus\_.+/;
				$headerTaxo=$1;
				chomp $Seq;
			$HEADERS{$headerTaxo} 	       = $header if defined $header;
			$CONTENTSsamekeys{$headerTaxo} = $Seq if defined $Seq;

		}
		my $Ntaxo=@TaxoOrder; #printing to 
		
		say OUTLOG "\n== C) STATISTICS ABOUT CONSENSUS SEQUENCES ===========================\n";
		printf OUTLOG ("%4s%4s%4s  %4s\n", "Len", "N", "-", "Taxonomic levels");
		say OUTLOG "----------------------------------------";
		for (my $i=0;$i<@TaxoOrder;$i++){
			my $H=$HEADERS{$TaxoOrder[$i]};
			$H=~/\_using_.+_consensus\_(.+)bp\_N(.+)\_D(\d+)\_.+/;
			printf  OUTLOG ("%4s%4s%4s  %4s\n", $1, $2, $3, $TaxoOrder[$i]);
		}
		say OUTLOG "----------------------------------------";
		say OUTLOG "[Len]: length (bp) of the consensus sequences; [N] number of unknown bases,";
		say OUTLOG "i.e. not unique A/T/G/C; [-]: number of bases present in the reference sequence";
		say OUTLOG "but absent in the consensus sequence (deletions).";
		say OUTLOG "NOTE: if too many Ns are present (e.g. >3-5), \nthe respective consensus sequence is not of good quality.";
		say OUTLOG "\n== D) CONSENSUS SEQUENCES =============================================\n";
		for (my $i=0;$i<@TaxoOrder;$i++){
			say OUTLOG "== Reference group: $TaxoOrder[$i]";
			print OUTLOG "$_\n" for unpack '(A60)*', $CONTENTSsamekeys{$TaxoOrder[$i]};
			say OUTLOG "";
			
		}
		close OUTLOG;
		#printing consensus seqs to 1 file
		open OUTCONS,">$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename"; # all consensus sequences there
		for (my $i=0;$i<@TaxoOrder;$i++){
			say OUTCONS "$HEADERS{$TaxoOrder[$i]}\n$CONTENTSsamekeys{$TaxoOrder[$i]}";   
		}
	} #	if(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/")
	} #unless(-e "$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename"
	}	#foreach my $BC 
}
sub R9_check_seq_blast_taxo{
	my ($args)=@_;
	my $Verbose			=$args->{Verbose};
	my $IN				=$args->{IN};#consensus_from_multi_ref_Bam.fasta 
	my $OUT				=$args->{OUT};#$Cdir/2_individual_barcodes/$BC
	my $REF_FP			=$args->{REF_FP}; #$REF_FP=$RefDB{$DBchoice}
	my $NberThreads		=$args->{NberThreads};#same
	my $PERL			=$args->{PERL};
	my $BLASTN			=$args->{BLASTN};

	#----------
	chdir $OUT; 
	my @Fparse=File::Basename::fileparse($REF_FP);
	my $blastDBFile= $Fparse[0]; #say "Fasta reference: $blastDBFile";
	my $blastDBFolder=$Fparse[1]; #say "Fasta reference folder: $blastDBFolder";
	system("export BLASTDB=$blastDBFolder;$BLASTN -query $IN -db $blastDBFile -outfmt 3  -out blast_nt_output_3.txt -num_threads $NberThreads   -num_descriptions 10 -max_hsps 30  -num_alignments 10 "); 
		system("$PERL -ne \'print if \/\^Query=\/ \.\. \/\^Lambda\/\' blast_nt_output_3.txt \| grep -v \"Lambda\" > blast_alignment.txt ");
		unless ($Verbose){system("rm -f blast_nt_output_3.txt");}
	system("export BLASTDB=$blastDBFolder;$BLASTN -query $IN -db $blastDBFile -outfmt \"7 pident nident length mismatch gaps sseqid \" -out blast_nt_output_7.txt -num_threads $NberThreads   -max_target_seqs 10 "); 
	system(" $PERL -pi -e \'s\/\# 10 hits found\\n\/\/\' blast_nt_output_7.txt ");
	system(" $PERL -pi -e \'s\/\^\# BLASTN.+\\n\/\/\' blast_nt_output_7.txt ");
	#system(" $PERL -pi -e \'s\/\^\# Fields.+\\n\/\\% ident\\/Nber ident pos\\/alignt length\\/N mismatches\\/N gaps\\/subject id\\/g\' blast_nt_output_7.txt ");
	system(" $PERL -pi -e \'s\/\^\# BLAST processed.+\\n\/\/\' blast_nt_output_7.txt ");
}

sub R9_2_check_seq_phylogeny{
	my ($args)=@_;
	my $Verbose			=$args->{Verbose};
	my $IN1				=$args->{IN1};#"$Cdir/2_individual_barcodes/$BC/02_individual_consensus/$Allconseqfilename",
	my $IN2				=$args->{IN2};#"$Cdir/2_individual_barcodes/$BC/03_BLAST_analysis/blast_nt_output_7.txt",
	my $OUT				=$args->{OUT};#$Cdir/2_individual_barcodes/$BC/02_individual_consensus
	my $NberThreads		=$args->{NberThreads};#same
	my $PERL			=$args->{PERL};
	my $SEQKIT			=$args->{SEQKIT};
	my $mafft			=$args->{mafft};
	my $Gblocks			=$args->{Gblocks};
	my $iqtree			=$args->{iqtree};
	my $REF_FP			=$args->{REF_FP}; #$REF_FP=$RefDB{$DBchoice}	
	#----------
	chdir $OUT; 
	#Goal: Alignment of all consensus with the REF1.fasta for each consensus
	#getting all Reference across subfolders of 02_individual_consensus/ 
	open IN2,"<$IN2" or die "cannot open IN2 in R9_2_check_seq_phylogeny";
	my (@REFAN,@SPLIT,@OUTGROUP);
	#use AN at the bottom of the blast7 output for each query to use as outgroup
		my $J=0;
	foreach (<IN2>){
		unless(/\#/) {
			chomp;
			@SPLIT=split/\t/,$_;
			my $SPL=pop @SPLIT;
			push(@REFAN,$SPL); #accumulating the ref
			$J++;
			if($J == 10){ push(@OUTGROUP,$SPL);} # based on the output of IN2 with 10 lines per query
		} else{ 
			$J=0;
		}
	}
	#remove outgroups from REF sequences (to provide a new order)
	chomp @OUTGROUP;
	chomp @REFAN;
	foreach my $O (@OUTGROUP){
		@REFAN = grep(!/$O/,@REFAN);
	}
	#retrieve outgroup sequences
	my $OUTGREF="$OUT/OUTREFs.fasta";
	foreach(@OUTGROUP){ #get REF sequences
		system("cat $REF_FP | $SEQKIT grep -n -r -p $_ >> $OUTGREF");	
	}
	my $ALLREF="$OUT/AllBlastREFs.fasta";
	foreach(@REFAN){ #get REF sequences
		system("cat $REF_FP | $SEQKIT grep -n -r -p $_ >> $ALLREF");	
	}
	my $ALLREF1="$OUT/AllBlastREFsUniq.fasta";
	system("$SEQKIT rmdup --quiet -n $ALLREF > $ALLREF1");
	#merge Cons and REF1 
		if(-e $OUTGREF) {system("cat $OUTGREF $IN1 $ALLREF1  > $OUT/CONSREF1.fasta");}
		#system("rm -f $ALLREF");
	#separating the headers from the sequences of the consensus fasta (easier to view)
		my $FASTAFILEHEADER="$OUT/CONSREF1.header";
		my $FASTAFILESimplifiedFASTA="$OUT/CONSREF1_simple.fasta";
		open F,"<$OUT/CONSREF1.fasta" || die "cannot open $OUT/CONSREF1.fasta"; 
		open OUTH,">$FASTAFILEHEADER"; # to save the fasta headers
		open OUTS,">$FASTAFILESimplifiedFASTA"; # to save the simplified fasta sequences
		my $i=1;
		foreach(<F>){
			if(/>/) {
				$_=~m/^>(.+)$/;
				say OUTH "Seq_[$i]: $1";
				print OUTS ">Seq_$i\n";
				$i++;
			}else{
				print OUTS $_;
			}
		}
	# align with clustal output
	unless(-e "$OUT/01_alignment") {system(mkdir "$OUT/01_alignment");}
	unless(-e "$OUT/01_alignment/mafftWithHeaders.clw"){
	 system("$mafft --quiet --clustalout --maxiterate 1000 --localpair --thread 1  $FASTAFILESimplifiedFASTA >  $OUT/mafft.clw");
	 system("$mafft --quiet --maxiterate 1000 --localpair --thread 1  $FASTAFILESimplifiedFASTA >  $OUT/mafft.fasta");# for phylogenic tree
	 system("$PERL -p -i -e 'tr/[atgcn]/[ATGCN]/' $OUT/mafft.clw");
	 system("$PERL -p -i -e 'tr/[atgcn]/[ATGCN]/' $OUT/mafft.fasta");
	 #merge headers and mafft
	 system("cat $FASTAFILEHEADER $OUT/mafft.clw > $OUT/mafftWithHeaders.clw");
	 system("$PERL -p -i -e 's/CLUSTAL.*//' $OUT/mafftWithHeaders.clw");
	 unlink "$OUT/mafft.clw";
	 system("mv $OUT/mafft.* $OUT/01_alignment/");
	 system("mv $OUT/mafftWithHeaders.clw $OUT/01_alignment/");
	}
	 #Gblocks
	 unless(-e "$OUT/01_alignment/mafft.fasta-gb"){
		system("$Gblocks $OUT/01_alignment/mafft.fasta -t=d >  /dev/null 2>&1");
	  }
	 #iqtree  produce phylogeny
my $treeINFO= <<TREEINFO;
Phylogenetic analysis results produced by iqtree:
================================================
*.iqtree  	IQ-TREE report:                
*.treefile  Maximum-likelihood tree:       
*.mldist  	Likelihood distances:          

Ultrafast bootstrap approximation results:
==========================================
*.contree  Consensus tree in Newick format
*.log  Screen log file               
TREEINFO
	 
	unless(-e "$OUT/02_phylogeny") {system(mkdir "$OUT/02_phylogeny");}	 
	unless(-e "$OUT/02_phylogeny/mafft.fasta-gb.treefile"){
		 system("$iqtree -s $OUT/01_alignment/mafft.fasta-gb -m GTR+I+G -bb 1000 -czb  -ntmax 20 -quiet");
			unlink "$OUT/01_alignment/mafft.fasta-gb.ckp.gz";
			unlink "$OUT/01_alignment/mafft.fasta-gb.splits.nex";
		 system("mv $OUT/01_alignment/*.iqtree $OUT/02_phylogeny/");
		 system("mv $OUT/01_alignment/*.treefile $OUT/02_phylogeny/");
		 system("mv $OUT/01_alignment/*.mldist $OUT/02_phylogeny/");
		 system("mv $OUT/01_alignment/*.contree $OUT/02_phylogeny/");
		 system("mv $OUT/01_alignment/*.log $OUT/02_phylogeny/");
		 system("mv $OUT/01_alignment/*.bionj $OUT/02_phylogeny/");
		 open TINFO,">$OUT/02_phylogeny/README_file_info.txt";
		 say TINFO $treeINFO;
		 close TINFO;
	}

}
sub R10_createTxtReport{ #combines log.txt and blast output
	my ($args)=@_;
	my $Verbose					=$args->{Verbose};
	my $IN1						=$args->{IN1}; #"log.txt"
	my $ID						=$args->{ID}; #$ID=$BCinfo{$BC};
	my $IN3						=$args->{IN3}; #"blast_nt_output_7.txt"
	my $IN4						=$args->{IN4}; #"blast_alignment.txt"
	my $IN5						=$args->{IN5}; #"$Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/CONSREF1.header",
	my $IN6						=$args->{IN6}; #$Cdir/2_individual_barcodes/$BC/04_phylogenetic_analysis/02_phylogeny/mafft.fasta-gb.iqtree",
	my $OUT						=$args->{OUT}; #$Cdir/2_individual_barcodes/$BC/$BC.report
	my $BC						=$args->{BC};
	my $TotalNberReadsinFasta	=$args->{TotalNberReadsinFasta};
	my $REF_FP					=$args->{REF_FP};
	my $DBchoice				=$args->{DBchoice}; #"BiBi16S"
	#----------------
	my $TotalNberReadsinREF=qx( grep -c ">" $REF_FP); chomp $TotalNberReadsinREF;
	open OUTBC,">$OUT";
	say OUTBC "===================================================================";	
	say OUTBC "===================================================================";	
	#say OUTBC "=================== !!!! NON ACCREDITED !!!! ======================";	
	say OUTBC "===================================================================";
	say OUTBC "***          REPORT OF SEQUENCED PCR AMPLICONS (NANOPORE)       ***";
	say OUTBC "===================================================================";
	say OUTBC "Analysis of barcoded sample ($BC): $ID";
	say OUTBC "Date of report: ",scalar(localtime);
	say OUTBC "===================================================================";
	say OUTBC "Total number of reads in the input fasta file: \*$TotalNberReadsinFasta\*";
	say OUTBC "The reference database \*$DBchoice\* contains \*$TotalNberReadsinREF\* sequences";

	open(LOG1,"<$IN1") or die "cannot open log.txt for barcode sample $BC";
	foreach (<LOG1>){print OUTBC $_;}
	
	say OUTBC "\n== E) Checking the taxonomic classification of the obtained consensus sequences \n(using Blastn on $DBchoice database)";
	open(R10in1,"<$IN3");foreach(<R10in1>)	{print OUTBC $_;}
	say OUTBC "====================================================================";
	open(R10in2,"<$IN4");foreach(<R10in2>)	{print OUTBC $_;}
	say OUTBC "====================================================================";
	say OUTBC "== F) PHYLOGENETIC TREE ============================================";
	say OUTBC "====================================================================";
	open(R10in3,"<$IN5");#headers
	my %TaxoLabels;
	foreach my $line (<R10in3>)	{
	 
		$line=~s/\[//;
		$line=~s/\]//;
		$line=~ m/(Seq_\d+)\:\s(.+)$/;
	    $TaxoLabels{$1} = $2 if defined $2;
	}
	#say OUTBC "====================================================================\n";
	open(R10in4,"<$IN6");
	my @R10IN4=<R10in4>;
	
	foreach my $R (keys %TaxoLabels) { #replacing the whole array
		s/$R$/$TaxoLabels{$R}/ for @R10IN4;
	}	
	if (grep {$_ eq "Type of analysis"} @R10IN4) {print OUTBC $_;}
	if (grep {$_ eq "Model of substitution"} @R10IN4) {print OUTBC $_;}
	my $TREEPOS1=0;
		foreach(@R10IN4){
			$TREEPOS1++;
			if(/MAXIMUM LIKELIHOOD TREE/){last;}
		}
	my $TREEPOS2=0;
		foreach(@R10IN4){
			if(/^Tree in newick format\:$/){last;}
			$TREEPOS2++;
		}
	for(my $L=$TREEPOS1-1;$L<$TREEPOS2;$L++){
		print OUTBC $R10IN4[$L];
	}			
	say OUTBC "====================================================================\n";
	
}
sub R11_createTxtReportNo{ #print log.txt if no consensus was built
	my ($args)=@_;
	my $Verbose					=$args->{Verbose};
	my $BC						=$args->{BC}; #"BC16"
	my $IN						=$args->{IN}; #"$Cdir/2_individual_barcodes/$BC/log.txt",
	my $ID						=$args->{ID}; #$ID=$BCinfo{$BC};	
	my $OUT						=$args->{OUT}; #$Cdir/2_individual_barcodes/$BC/$BC.report
	my $TotalNberReadsinFasta	=$args->{TotalNberReadsinFasta};
	my $REF_FP					=$args->{REF_FP};
	my $DBchoice				=$args->{DBchoice}; #"BiBi16S"
	#----------------
	my $TotalNberReadsinREF=qx( grep -c ">" $REF_FP); chomp $TotalNberReadsinREF;
	open OUTBC,">$OUT";
	say OUTBC "===================================================================";	
	say OUTBC "===================================================================";	
	#say OUTBC "============= !!!! NON ACCREDITED !!!! ============================";
	say OUTBC "===================================================================";
	say OUTBC "***          REPORT OF SEQUENCED PCR AMPLICONS (NANOPORE)       ***";
	say OUTBC "===================================================================";
	say OUTBC "Analysis of barcoded sample ($BC):  $ID";
	say OUTBC "Date of report: ",scalar(localtime);
	say OUTBC "================================================";
	say OUTBC "No sequence reference was kept with enough reads matching.\nNo consensus was built for this sample.\n";
	say OUTBC "================================================";
	say OUTBC "Total number of reads used from the modal sequences: \*$TotalNberReadsinFasta\* (modal sequences)";
	say OUTBC "The reference database \*$DBchoice\* contains \*$TotalNberReadsinREF\* sequences";
	say OUTBC "================================================";
	open IN,"<$IN";
	foreach(<IN>){print OUTBC $_;}
}
sub R11_CreatePDFfromTXT{
	my ($args)=@_;
	my $Verbose		=$args->{Verbose};
	my $IN			=$args->{IN};
	my $OUT			=$args->{OUT};
	#----------------
	use PDF::API2;
	use strict;
	use warnings;
	use utf8;
	use feature qw(say);
	use POSIX;
	#----------------
	unless($OUT =~ /\.pdf/){say "much choose a .pdf file as output file! Exiting.";exit;}
	open (IN,$IN) or die "cannot open the input file to print as PDF";
	my @TEXT=<IN>; 
	my $NberLines=scalar @TEXT;
	my $SampleBC="??";
#----------------#adding some blanks for the tree to be displayed on a next page
		my $CL=0;
		my $JUMP=0;
	foreach(@TEXT){
	 if($_=~/Numbers in parentheses /){$JUMP=$CL; } #to print on a new page!
	 $CL++;
	 if($_=~ /Analysis of barcoded sample (.+)$/){$SampleBC=$1;}
	}
	# say "NberLines=",$NberLines;
	# say "JUMP=",$JUMP;
	# say "ceiling JUMP=",ceil($JUMP/67);
	# say "Diff=",$NberLines-$JUMP;
	my $Nblanks=ceil($JUMP/67)*67 - $JUMP;
	#say "nber blank lines to add: ",$Nblanks;
	#modifying the @TEXT to include more lines
	my @blanks=("\n");
	for (1..$Nblanks){push(@blanks,"\n");}
	my @TEXT_spaced=@TEXT;
	splice @TEXT_spaced,$JUMP+1,0,@blanks;
	$NberLines=scalar @TEXT_spaced; #updating
	#1+(0:2)*67 for top of the page
	#(1+(0:2)*67 ) -1 for bottom of the page
#----------------
	#$NberLines+=20; #if extra lines are needed to be added
	my $NPages=ceil($NberLines/65);	#65 lines per page
	#say "NPages:$NPages";
	
	 #calculating absolute Y positions to print lines on each page
		my $J=720;
		my @Ypos;
		while($J>50) {
			push(@Ypos,$J);
			$J=$J-10; 
		}
		my $NbYpos=scalar (@Ypos);
	#working on the footer text
		my $SampleDate=scalar(localtime);
	#----------------	
    my $pdf = PDF::API2->new();# Create a blank PDF page
		my $Y=0;
		my $CurrentPageNber=0;
		my ($page,$font,$text);
		for (my $i=0;$i<$NberLines;$i++) {
			if ($Y == 0) {
					$page = $pdf->page();# Add a blank page
					$page->mediabox('Letter');
					$font = $pdf->corefont('Courier');
					$text = $page->text();
					$text->font($font, 11);
			}
			$text->translate(40,$Ypos[$Y]);
			#$text->text($i." ".$TEXT_spaced[$i]); 
			$text->text($TEXT_spaced[$i]);
			$text->translate(150,30); 
			my $FT="$SampleBC \($SampleDate\) - $CurrentPageNber/$NPages";
			$text->text($FT);
		$Y++;
		if ($i % $NbYpos == 0) {$Y=0;$CurrentPageNber++;}
	 }
	$pdf->saveas($OUT);
}
	
1;
