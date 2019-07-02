# Preparation of custom 16S database 

We will use [mothur](https://www.mothur.org/) to subset long database sequences to specific region of interest (typically the region amplified by the primers).   

We need:  
- **the database of long sequences (here 16S) to customize**    
e.g. [db-BiBI](https://umr5558-bibiserv.univ-lyon1.fr/BIBIDOCNEW/db-BIBI.html#availability)    
SSU-rDNA-mk37_stringent/16S_stringent_dedup.fasta    
- **primer sequences to delineate the region of interest**   
	  AGAGTTTGATCNTGGCTCAG\<_sequence of interest_\>THYGTGCCAGCWGCCGCGGTA   
- **reference sequence on which the positions will be defined**   
[J01859.1](http://www.ncbi.nlm.nih.gov/nuccore/174375?report=fasta) Escherichia coli 16S ribosomal RNA, complete sequence      
- **a reference database of aligned sequences (aka seed, or template alignment)**    
[silva/Refv132/seed/silva.seed_v132.align](https://mothur.org/wiki/Silva_reference_files)   

Make sure that all files are located in the same folder and that you are running _mothur_ from that folder.    

### 0) sanity check   
make sure that the primer sequences indeed match the target db sequences. E.g.       
```
DB=16S_stringent_dedup # no .fasta suffix
grep -c "AGAGTTTGATC[ATCG]TGGCTCAG" $DB.fasta  			 #50,860
grep -c "T[ATCG][ATCG]GTGCCAGC[ATCG]GCCGCGGTA" $DB.fasta #155,989 
```


### 1) Define primer sites on the reference sequence    
Can be done manually oruse a primer aligning software (e.g. [Multalin](http://multalin.toulouse.inra.fr/multalin/))     
find start-end positions of primers
	AGAGTTTGATCNTGGCTCAG...THYGTGCCAGCWGCCGCGGTA

>**J01859.1 Escherichia coli 16S ribosomal RNA, complete sequence**
```
       AGAGTTTGATCNTGGCTCAG
       ||||||||||||||||||||    
AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGT
AACAGGAAGAAGCTTGCTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATG
GAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCG
GGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACG
ATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGG
CAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTT
CGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCG
                    THYGTGCCAGCWGCCGCGGTA
                    |||||||||||||||||||||
CAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAAT
TACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAAC
TGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGT
AGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCG
TGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCC
TTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACT
CAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCT
TACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGC
TGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCT
TTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGA
CGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGA
CCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATG
AAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCG
CCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTT
TGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
```
		
>**TargetAmplicon.fasta**    
```
ATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGT
AACAGGAAGAAGCTTGCTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATG
GAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCG
GGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACG
ATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGG
CAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTT
CGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCG
CAGAAGAAGCACCGGCTAAC
```


### 2) Align DB sequences to template sequence db    
```
module load UHTS/mothur/1.40.5   

mothur "#align.seqs(candidate=$DB.fasta,template=silva.seed_v132.align,processors=20)"> alignseq.txt

#Output File Names: 
# $DB.align   # 11 Gb
# $DB.align.report
# $DB.flip.accnos
```
This may take some time (~30 min; 2014 secs to align 233521 sequences).       

### 3) Find the positions of the target amplicon considering the positions of aligned sequences in the db   
```
mothur "#align.seqs(fasta=TargetAmplicon.fasta,reference=$DB.align,processors=20)"
#Output File Names:
# TargetAmplicon.align
# TargetAmplicon.align.report

mothur "#summary.seqs(fasta=TargetAmplicon.align)"
#	                Start   End     NBases  Ambigs  Polymer NumSeqs
#	Minimum:        1163    11885   440     0       5       1
#	2.5%-tile:      1163    11885   440     0       5       1
#	25%-tile:       1163    11885   440     0       5       1
#	Median:         1163    11885   440     0       5       1
#	75%-tile:       1163    11885   440     0       5       1
#	97.5%-tile:     1163    11885   440     0       5       1
#	Maximum:        1163    11885   440     0       5       1
#	Mean:           1046    11885   440     0       5
#	# of Seqs:      1
```
It takes few mins and depends on the size of $DB.align (4.5 min for 233521 seqs with 20 processors).   
Therefore the TargetAmplicon sequence is found from positions 1163 to 11885.     

### 4) Examine with pcr.seqs the coordinates above    
```
mothur "#pcr.seqs(fasta=$DB.align, ecoli=TargetAmplicon.align, keepdots=FALSE, processors=20,nomatch=reject)"
#Output File Names:
# $DB.pcr.align
# $DB.scrap.pcr.align

mothur "#summary.seqs(fasta=$DB.pcr.align)"
#                 Start   End     NBases  Ambigs  Polymer NumSeqs
# Minimum:        1       573     14      0       3       1
# 2.5%-tile:      1       10722   369     0       4       5836
# 25%-tile:       4       10722   415     0       4       58359
# Median:         4       10722   436     0       5       116717
# 75%-tile:       6       10722   447     0       5       175075
# 97.5%-tile:     911     10722   465     2       6       227597
# Maximum:        9143    10722   643     165     165     233432
# Mean:   131     10718   428     0       4
# # of Seqs:      233432

```

Extract the coordinates (median) from the output:
```
St=`mothur "#summary.seqs(fasta=$DB.pcr.align)" | grep "Median"| cut -f2`
En=`mothur "#summary.seqs(fasta=$DB.pcr.align)" | grep "Median"| cut -f3`
```
Here, St=4 and En=10722.   

### 5) Restrict the db at the defined positions
```
mothur "#screen.seqs(fasta=$DB.pcr.align,start=$St, end=$En)"
#Output File Names:
# $DB.pcr.good.align
# $DB.pcr.bad.accnos

# Check here that the stats are OK: 
mothur "#summary.seqs(fasta=$DB.pcr.good.align)"
#                 Start   End     NBases  Ambigs  Polymer NumSeqs
# Minimum:        1       10722   351     0       3       1
# 2.5%-tile:      1       10722   411     0       4       4338
# 25%-tile:       1       10722   433     0       4       43377
# Median:         4       10722   440     0       5       86753
# 75%-tile:       4       10722   448     0       5       130129
# 97.5%-tile:     4       10722   469     2       6       169167
# Maximum:        4       10722   643     165     165     173504
# Mean:           3       10722   438     0       4
# # of Seqs:      173504

```

### 6) Filtering out alignment gaps 
```		
mothur "#filter.seqs(fasta=$DB.pcr.good.align, vertical=T, trump=.)" | grep -A 7 "Length of filtered alignment:"  > Alignment.summary

more Alignment.summary
# Length of filtered alignment: 2404
# Number of columns removed: 8318
# Length of the original alignment: 10722
# Number of sequences used to construct filter: 173504
# 
# Output File Names:
# $DB.filter
# $DB.pcr.good.filter.fasta
```		
$DB.pcr.good.filter.fasta  contains the aligned, truncated sequences.    

### 7) Deduplicating sequences
```
mothur "#unique.seqs(fasta=$DB.pcr.good.filter.fasta)" 
#Output File Names:
# $DB.pcr.good.filter.names  		# identical names are on the same line; separated by a tab
# $DB.pcr.good.filter.unique.fasta  # same nbers of sequences

#grep -c ">" $DB.pcr.good.filter.unique.fasta #63,407 seqs
```

### 8) Renaming main output files   
A) Sequences  
```
mv $DB.pcr.good.filter.unique.fasta custom_aligned.fasta
cp custom_aligned.fasta custom.fasta


# to remove extra info about alignment positions
perl -p -i -e "unless(/>/){s/[\-\.]//g}" custom.fasta 
```

B) Sequence names in a csv file  
```
cut -f2 $DB.pcr.good.filter.names > custom.synonymous.names
```  


C) Extra cleaning of the names file and fasta header; removing the information after "=..."

CleanNameFileLeBiBi.pl:
```		
open IN,"<$ARGV[0]";#custom.synonymous.names
while(<IN>){   
 @D=split(/\,/,$_);   
 s/\=.+//g for @D;  
 print join("\,",@D);  
}     
```		

Removing the first element of each line   
```		
perl CleanNameFileLeBiBi.pl custom.synonymous.names > custom_simplified.names
```		


### 9) Main results
3 files:    
- **custom_aligned.fasta**     
  Aligned sequences with alignment gaps (-)   
- **custom.fasta**      
  Multifasta of unaligned sequences for further downstream applications    
- **custom_simplified.names**    
  Name file that contains all identical sequences, with the first one being the one in the fasta file   
  
  
### 10) Preparation of the BLAST database
```	
module load Blast/blast/2.6.0
makeblastdb -in  custom.fasta -dbtype nucl   
```	  
then edit the lib/config.pm file to point to the fasta file   
```	  
our %RefDB=(
BiBi16S => "/pathto/leBiBicustomDB/custom.fasta",
);
```	  


