#!/usr/bin/perl
# Michael Hiller, 2023
# filters a sam file for mapping quality using the CIGAR string
# uses a build-in two state HMM to segment the polyA tail from the end of a read

use strict;
use warnings;
use diagnostics;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use POSIX qw(Inf NaN);


$| = 1;		# == fflush(stdout)
my $verbose = 0;
my $verboseHMM = 0;
my $statFile = 0;		# flag: if switched on, produce a tsv file with statistics for all reads
my $keepBad5Prime = 0;		# flag
my $polyAReadSuffix = 0;	# if >0, compute and output the length of the polyA tail in the $polyAReadSuffix + 3'clip_len suffix of the read
my $suffixStepSize = 50;
# alignment quality filter parameters
my $minPerID=98;
my $maxClip5=20;
my $maxClip3=20;
my $minPolyAlen = 30;
# HMM parameters
my $P2P = 0.9;		# transition prob for polyA tail
my $emitA = 0.99;	# emission prob for A in polyA state

# HMM matrix for viterbi and traceback
my @V;	# dynamic programming matrix; every cell holds the log score
my @T; 	# traceback matrix. Cell [i,j] holds the index [k,l] of the previous cell that gave the Viterbi maximum
# and other arrays for the HMM
my @state2Name;
my @state2IncomingTransitions;
my @state2IncomingTransitionsScore;


# options
my $usage = "usage: $0 input.sam [-perID int[0-100] -clip5 int -clip3 int -statFile file -polyAReadSuffix int -P2P double -emitA double -v[verbose]  -v1[verboseHMM]]
 Input is a sam file, must be computed with -eqx (X as mismatch in the CIGAR) otherwise the perID value is wrong
 Output are two files
	{input}.good.sam  for sam entries that pass the filters.
	{input}.bad.sam   for sam entries that do not pass. 
	PerID, 5/3 clip and polyA length is added to the read ID 
 \
# mapping filter parameters
 -perID int      min %id (computed without the 5' and 3' clip). Must be [0-100] (percent). Default $minPerID  
 -clip5 int      max 5' soft or hard clip. Default $maxClip5
 -clip3 int      max 3' soft or hard clip. Default $maxClip3
# polyA=TES filters
 -minPolyAlen int min length of the unaligning polyA tail to be considered a real TES. Default $minPolyAlen 
# polyA HMM parameters
 -P2P double     transition probability of looping in the polyA state (default $P2P)
 -emitA double   probability of emitting A in the polyA state (default $emitA)
# other parameters
 -statFile       flag: if set output a {input}.tsv file that contains statistics of all reads
                    readID {tab} perID {tab} 5'clip {tab} 3'clip (non-PolyA part){tab} polyA tail length of 3'clip [optionally]{tab} polyA tail length of read
 -keepBad5Prime  flag: if set, produce an {input}.TESBad5Prime.bed file which contains reads that align well but only have a 5' clip above our threshold. Can be used for PAScaller
 -polyAReadSuffix int  if >0, compute and output the length of the polyA tail in the \$polyAReadSuffix + 3'clip_len suffix of the read (default don't do that)\n";

################################
GetOptions ("v|verbose"  => \$verbose, "v1|verboseHMM"  => \$verboseHMM, "keepBad5Prime" => \$keepBad5Prime, "statFile" => \$statFile, 
   "minPolyAlen=i" => \$minPolyAlen, "polyAReadSuffix=i" => \$polyAReadSuffix, 
   "perID=i"  => \$minPerID, "clip5=i"  => \$maxClip5, "clip3=i"  => \$maxClip3,
   "P2P=f" => \$P2P, "emitA=f" => \$emitA) || die "$usage\n";
die "$usage\n" if ($#ARGV < 0);
die "ERROR: perID ($minPerID) must be [0,100]\n" if ($minPerID < 0 || $minPerID > 100);
die "ERROR: clip5 ($maxClip5) must be >= 0\n" if ($maxClip5 <0);
die "ERROR: clip3 ($maxClip3) must be >= 0\n" if ($maxClip3 <0);
die "ERROR: polyAReadSuffix ($polyAReadSuffix) must be > 0\n" if ($polyAReadSuffix < 0);
die "ERROR: P2P ($P2P) must be [0,1]\n" if ($P2P < 0 || $P2P > 1);
die "ERROR: emitA ($emitA) must be [0,1]\n" if ($emitA < 0 || $emitA > 1);

# pre-compute log emission probs for the HMM
my $logemitA = log($emitA);
my $logemitNotA = log( (1-$emitA) / 3);  # remaining prob
my $logUniform = log(0.25);
print "Transition PolyA->PolyA state: $P2P    Prob emit A in polyA state: $emitA\n" if ($verbose);

################################
# output files
die "ERROR: input file $ARGV[0] does not end with .sam\n" if ($ARGV[0] !~ /.*\.sam$/);
my $filePrefix = substr($ARGV[0], 0, rindex($ARGV[0], ".sam"));
my $outGood = $filePrefix . ".good.sam";
my $outBad = $filePrefix . ".bad.sam";
my $outTES = $filePrefix . ".TES.bed";
my $outTESBad5Prime = $filePrefix . ".TESBad5Prime.bed";
#my $outnoTES = $filePrefix . ".noTES.bed";
my $outStat = $filePrefix . ".tsv";
print "read $ARGV[0] and write output to $outGood + $outBad and TES to $outTES\n";
print "output $outTESBad5Prime which contains reads that align well but only have a 5' clip above our threshold\n" if ($keepBad5Prime);
die "ERROR: output file $outGood already exists. Pls delete it or rename it\n" if (-e $outGood);
die "ERROR: output file $outBad already exists. Pls delete it or rename it\n" if (-e $outBad);

################################
# process the input
open(fileIn, $ARGV[0]) || die "ERROR: cannot open $ARGV[0]\n";
open(fileSAMGood, ">$outGood") || die "ERROR: cannot write to $outGood\n";
open(fileSAMBad, ">$outBad") || die "ERROR: cannot write to $outBad\n";
open(fileBedTES, ">$outTES") || die "ERROR: cannot write to $outTES\n";
open(fileBedTESBad5Prime, ">$outTESBad5Prime") || die "ERROR: cannot write to $outTESBad5Prime\n";
#open(fileBednoTES, ">$outnoTES") || die "ERROR: cannot write to $outnoTES\n";
if ($statFile) {
	open(fileTSVStat, ">$outStat") || die "ERROR: cannot write to $outStat\n";
	if ($polyAReadSuffix > 0) {
		print fileTSVStat "readID\tperID\tclip5\tclip3 (non-PolyA part)\tpolyA tail length of 3'clip\tpolyA tail length of read\n";
	}else{
		print fileTSVStat "readID\tperID\tclip5\tclip3 (non-PolyA part)\tpolyA tail length of 3'clip\n";
	}
}
my $line1;
while ($line1 = <fileIn>) {
	chomp($line1);
	# copy comments into the 2 output files
	if ($line1 =~ /^@.*/) {
		print fileSAMGood "$line1\n";
		print fileSAMBad "$line1\n";
		next;
	}

	my @f = split(/\t/, $line1);
	my $ID = $f[0];
	my $strandFlag = $f[1];
	my $chrom = $f[2];
	my $genomicStart = $f[3];
	my $CIGAR = $f[5];
	my $seq = $f[9];

	print "\n$line1\n--> $ID $strandFlag $CIGAR\n" if ($verbose);
	# 4 == unmapped read
	if ($strandFlag == 4) {	
		print "==> skip unmapped read\n" if ($verbose);
		next;
	}
	# revcomp mapping is indicated by these 2 flags
	my $strand = "+";
	$strand = "-" if ($strandFlag == 16 || $strandFlag == 2064);
	my ($perID, $clip5, $clip3, $hasHardclip3, $TES) = parseCIGAR($CIGAR, $seq, $genomicStart, $strand);

	# revcomp mapping --> revcomp seq and swap 5' and 3' clip
	if ($strand eq "-") {
		$seq = revComp($seq);
		my $x = $clip5;
		$clip5 = $clip3;
		$clip3 = $x;
		print "--> - strand mapping (flag $strandFlag)\n--> $seq\n" if ($verbose);
	}

	# compute length of polyA tail
	my $polyAlen = 0;
	my $polyAseq = "";
	if ($clip3 > 0) {
		if ($hasHardclip3 == 0) {
			my $clipSeq = getSuffix($seq, $clip3);
			($polyAlen, $polyAseq) = segmentPolyA($clipSeq);
			print "--> 3' clip seq $clipSeq\n" if ($verbose);
		}else{
			print "--> has a hardclip at the 3' end\n" if ($verbose);
		}
	}
	# subtract polyA from 3'clip
	$clip3 = $clip3 - $polyAlen;
	print "--> perID: $perID,  clip5: $clip5,  clip3: $clip3,  polyAlen: $polyAlen  $polyAseq\n" if ($verbose);

	# determine the length of the polyA tail in the read by taking the unmapped portion + $polyAReadSuffix 
	# this difference informs us about the intrapriming potential
	my $polyAlenRead = -1; 
	my $polyAseqRead = "";
	if ($polyAReadSuffix > 0) {
		# we take the suffix of the read that has length $polyAReadSuffix + $clip3
		my $suffixLen = $polyAReadSuffix + $clip3 + $polyAlen;   # Note that we subtracted polyAlen above -- need to add it back
		($polyAlenRead, $polyAseqRead) = segmentPolyARead($suffixLen, $seq);

		# rarely the following happens
		# segmentPolyA.perl TGTGCCTATGTAATAAAGTCTATACACTGGCAAAAACACAAAAAAAAAAAAAAAAAAAA
		# TGTGCCTATGTAATAAAGTCTATACACTGGCAAAAACAC,AAAAAAAAAAAAAAAAAAAA,39,20
		# segmentPolyA.perl CAAAAACACAAAAAAAAAAAAAAAAAAAA
		# CAAAAACAC,AAAAAAAAAAAAAAAAAAAA,9,20
		# but segmenting only the unaligning 3'clip, we get 
		# segmentPolyA.perl AAAAACACAAAAAAAAAAAAAAAAAAAA
		# ,AAAAACACAAAAAAAAAAAAAAAAAAAA,0,28
		# --> edge cases. If polyA read < polyA 3'clip, then polyA read = polyA 3'clip
		if ($polyAlenRead < $polyAlen) {
			print "\tWARNING: inferred polyA length of read ($polyAlenRead) < ($polyAlen) inferred polyA length of 3' clip --> correct this\n" if ($verbose);
			$polyAlenRead = $polyAlen;
			$polyAseqRead = $polyAseq;
		}

		print fileTSVStat "$ID\t$perID\t$clip5\t$clip3", ($hasHardclip3==0 ? "" : "Hardclip"), "\t$polyAlen\t$polyAlenRead\t$polyAseq\t$polyAseqRead\n" if ($statFile);
		$ID = sprintf("%s_PerID%1.3f_5Clip%d_3Clip%d_PolyA%d_PolyARead%d", $ID, $perID, $clip5, $clip3, $polyAlen, $polyAlenRead);
	}else{
		print fileTSVStat "$ID\t$perID\t$clip5\t$clip3", ($hasHardclip3==0 ? "" : "Hardclip"), "\t$polyAlen\t$polyAseq\n" if ($statFile);
		$ID = sprintf("%s_PerID%1.3f_5Clip%d_3Clip%d_PolyA%d", $ID, $perID, $clip5, $clip3, $polyAlen);
	}
	
#	$ID = sprintf("%s,%1.3f,%d,%d,polyA%d", $ID, $perID, $clip5, $clip3, $polyAlen);
	$f[0] = $ID;
	$" = "\t";
	$line1 = "@f";
	print "result: $line1\n" if ($verbose);
	
	if ($perID*100 >= $minPerID && $clip5 <= $maxClip5 && $clip3 <= $maxClip3) {
		print "====> GOOD\n" if ($verbose);
		print fileSAMGood "$line1\n";
		print fileBedTES "$chrom\t",$TES-1,"\t$TES\t$ID\t1000\t$strand\n";		
	}else{
		if ($keepBad5Prime && $perID*100 >= $minPerID && $clip3 <= $maxClip3) {
			print fileBedTESBad5Prime "$chrom\t",$TES-1,"\t$TES\t$ID\t1000\t$strand\n";
		}
		print fileSAMBad "$line1\n";
	}
	
}
close fileIn;
close fileSAMGood;
close fileSAMBad;
close fileBedTES;
close fileBedTESBad5Prime;
#close fileBednoTES;
if ($statFile) {
	close fileTSVStat;
}


#########################################################################################################
# take the read suffix and infer the polyA tail
# if the entire suffix is polyA, increase the suffix length by 50 bp
#########################################################################################################
sub segmentPolyARead {
	my ($suffixLen, $seq) = @_;
	
	my $polyAlenRead = -1;
	my $polyAseqRead = "";
	for (;;) {
		($polyAlenRead, $polyAseqRead) = segmentPolyA( getSuffix($seq, $suffixLen) );
		printf "--> polyA length of $suffixLen bp read suffix (%s): %d\n", getSuffix($seq, $suffixLen), $polyAlenRead if ($verbose);

		# continue with a larger suffix if the entire suffix is polyA
		if ($polyAlenRead == $suffixLen) {
			$suffixLen += $suffixStepSize;
		}else{
			return ($polyAlenRead, $polyAseqRead);
		}
	}
}

########################################################################################################
# save fct to get a suffix from a string that avoids underflows
# if the start pos, determined as seqLen - suffix, is <0, then this function returns the entire string.
# Example with perl substr
#   string: 1234567890  length=10 suffix = 12
#   substr($seq, length($seq)-$suffix, $suffix)   returns '90'
########################################################################################################
sub getSuffix {
	my ($seq, $suffix) = @_;

	my $l = length($seq);
	my $start = $l - $suffix;
	$start = 0 if ($start < 0);
	return substr($seq, $start, $suffix);
}


#########################################################################################################
# parse the CIGAR string and return %id, 5'3' clip and the length of the polyA tail
# also returns the TES (the genomic end position of the mapping)
#########################################################################################################
sub parseCIGAR {
	my ($CIGAR, $seq, $genomicStart, $strand) = @_;
	
	die "ERROR: CIGAR string contains 'M', indicating that it was not run with minimap -eqx\n" if ($CIGAR =~ /M/);

	my $hasHardclip3 = 0;   # flag to return if that read has a hard clipped 3' end --> if so

	# CIGAR will look like 216=47196N475=3403N213=2362N105=8540N43=5835N76=2055N219=8342N88=3829N153=2052N98=2086N173=4307N206=8983N2592=1I4178=1X3838=1X3=1X6=10S
	# replace =(match) N(intron) I(ins) D(del) X(mismatch) S(softclip) H(hardclip) with the letter+\t
	# then we can easily split it by tab
	$CIGAR =~ s/=/=\t/g;
	$CIGAR =~ s/X/X\t/g;
	$CIGAR =~ s/I/I\t/g;
	$CIGAR =~ s/D/D\t/g;
	$CIGAR =~ s/N/N\t/g;
	$CIGAR =~ s/S/S\t/g;
	$CIGAR =~ s/H/H\t/g;

	print "\tCIGAR: $CIGAR\n" if ($verbose);
	my @f = split(/\t/, $CIGAR);

	# compute %id
	my $match = 0;
	my $aliLength = 0;
	my $clip5 = 0;
	my $clip3 = 0;
	for (my $i = 0; $i <= $#f; $i++) {
		my $len = length($f[$i]);
		my $c = substr($f[$i], $len-1, 1);	# = X D I N S H
		my $x = substr($f[$i], 0, $len-1);	# integer
		die "ERROR: $x ($f[$i]) is not an integer\n" if ($x !~ /^\d+$/);

		# skip intron
		next if ($c eq "N");
		
		# 5' clip (then $aliLength == 0) or 3' clip
		if ($c eq "S" || $c eq "H") {
			if ($aliLength == 0) {
				$clip5 = $x;
			}else{
				$clip3 = $x;
				$hasHardclip3 = 1 if ($c eq "H");   # flag for 3' hardclip
			}
			next;
		}
		if ($c eq "=") {
			$match += $x;
		}
		die "ERROR: $c is not = X I or D\n" if (! ($c eq "=" || $c eq "X" || $c eq "I" || $c eq "D"));
		$aliLength += $x;
	}
	
	my $perID = $match / $aliLength;
	print "\t==> perID = match/aliLen:  $perID = $match / $aliLength\n" if ($verbose);

	# for a minus strand mapping, the genomic start pos from the sam file is the end position (1 based)
	# for plus strand mapping, we have to compute that from the CIGAR string
	my $TES = $genomicStart;
	$TES = getTES($genomicStart, $CIGAR) if ($strand eq "+");

	return ($perID, $clip5, $clip3, $hasHardclip3, $TES);
}




##############################################
# compute end pos from CIGAR for a plus strand mapping
# CIGAR is 1-based
#   for a minus strand mapping the given genomic coordinate is the end position in 0-based
#   for a plus strand mapping, we need to subtract -1
#
#    012345678      coordinates 0-based
#    123456789      coordinates 1-based
#      >>>>>        genomicStart = 3, say single match block of length 5= --> 3+5=8 which is the pos downstream of the read --> subtract 1
#        <<<        genomicStart = 5
##############################################
sub getTES {
	my ($genomicStart, $CIGAR) = @_;
	my @f = split(/\t/, $CIGAR);

	print "\tcompute TES from $genomicStart  and $CIGAR\n" if ($verbose);
	# walk the CIGAR
	for (my $i = 0; $i <= $#f; $i++) {
		my $len = length($f[$i]);
		my $c = substr($f[$i], $len-1, 1);	# = X D I N S H
		my $x = substr($f[$i], 0, $len-1);	# integer
		die "ERROR: $x ($f[$i]) is not an integer\n" if ($x !~ /^\d+$/);
	
		# do not count hard or soft clips
		next if ($c eq "S" || $c eq "H");
		# skip insertions in the read
		next if ($c eq "I");
		
		$genomicStart += $x;
		print "\t\t$f[$i]  --> genomicPos $genomicStart\n" if ($verbose);
	}
	return $genomicStart - 1;   # why? See above. Note, we only call getTES() for plus strand mappings
}	


##############################################
# return reverse complement seq
##############################################
sub revComp {
	my $seq = shift;
	$seq =~ tr/ATGCatgc-/TACGtacg-/;
	$seq = reverse($seq);
	return $seq;
}




###########################################################################
# HMM functions
###########################################################################

##############################################
# return reverse complement seq
##############################################
sub segmentPolyA {
	my $seq = shift;

	# reverse and append an X at seq[0] --> this allows us substr($seq, $pos, 1) 
	$seq = "X" . uc reverse($seq);
	my $seqLen = length($seq);
	print "\tPolyA: sequence of length: $seqLen\n$seq\n" if ($verboseHMM);

	# generate HMM
	undef @state2Name;
	undef @state2IncomingTransitions;
	undef @state2IncomingTransitionsScore;
	# state 0 --> polyA state
	# state 1 --> Sequence state (rest of the sequence)
	my $numberOfStates = generateHMM();

	# start fresh
	undef @V;	# dynamic programming matrix; every cell holds the log score
	undef @T; 	# traceback matrix. Cell [i,j] holds the index [k,l] of the previous cell that gave the Viterbi maximum

	###################
	# init: both P and S can be start states with equal prob
	$V[0][0] = log(0.5);
	$V[1][0] = log(0.5);

	###################
	# Viterbi
	print "\nstart dynamic programming for [1-$seqLen] sequence\n" if ($verboseHMM);
	for (my $seqPos=1; $seqPos < $seqLen; $seqPos ++) { 
		print "------------------------------------------------------------------------------------------\n" if ($verboseHMM);
		# we only have 2 states
		computeScore(0, $seqPos, $seq);
		computeScore(1, $seqPos, $seq);
	}

	# both P and S are also end states
	print "DONE: Max score: V[P][$seqLen-1] = ", $V[0][$seqLen-1], " vs. V[S][$seqLen-1] = ", $V[1][$seqLen-1],"\n\n\n" if ($verboseHMM);
	my $polyALen = -1;
	my $polyAseq = "";
	if ($V[0][$seqLen-1] > $V[1][$seqLen-1]) {
		print "Traceback from state P\n" if ($verboseHMM);
		# Traceback from polyA state
		($polyALen, $polyAseq) = traceback(0,$seqLen-1, $seq);
	}else{
		print "Traceback from state S\n" if ($verboseHMM);
		# Traceback from seq state
		($polyALen, $polyAseq) = traceback(1,$seqLen-1, $seq);
	}
	#printV($numberOfStates,$seqLen);
	return ($polyALen, $polyAseq);
}

##############################################
# print dynamic programming matrix
##############################################
sub printV {
	my ($stateNumber, $seqPos) = @_;
	for (my $i=0; $i<$stateNumber; $i++) {
		for (my $j=1; $j<$seqPos; $j++) {
			printf "$V[$i][$j]\t";
		}
		print "\n";
	}
}

##############################################
# do the traceback
##############################################
sub traceback {
	# end position in the dynamic programming matrix
	my ($stateNumber, $seqPos, $seq) = @_;

	my @path;
	my $endOfPolyA = -1;

	for (;;) {
		push @path, $state2Name[$stateNumber];

		my ($prevState, $prevPos) = (split(/,/, $T[$stateNumber][$seqPos]))[0,1];
		print "TRACEBACK: Transition $state2Name[$stateNumber] TO $state2Name[$prevState]\tT[$seqPos]:", $T[$stateNumber][$seqPos],"  seqPos $seqPos -> $prevPos\n" if ($verboseHMM);

		# check if we transition from seq state to polyA state
		if ($stateNumber == 1 && $prevState == 0) {
			print "\tFOUND end of polyA tail at position $seqPos\n" if ($verboseHMM);
			$endOfPolyA = $seqPos;
		}

		$stateNumber = $prevState;
		$seqPos = $prevPos;
		last if ($seqPos == 0);
	}

	# end of polyA in reversed seq (the original orientation)
	# if -1, the entire seq is a polyA or the entire seq is not polyA
	if ($endOfPolyA != -1) {
		$endOfPolyA = length($seq) - $endOfPolyA;
	}else{
		if ($stateNumber == 0) {
			$endOfPolyA = 0;    # we ended in the polyA state --> everything is polyA
		}elsif ($stateNumber == 1) {  # we ended in the seq state --> no polyA
			$endOfPolyA = length($seq) -1;   # seq still has the extra X
		}else{
			die "ERROR: stateNumber $stateNumber is neither 0 nor 1\n";
		}
	}
	# reverse again and chop of the added 'X'
	$seq = reverse($seq);
	chop($seq);
	$" = "";
	print "$seq\n@path\n" if ($verboseHMM);
	# this is the final output
	print substr($seq, 0, $endOfPolyA), ",", substr($seq, $endOfPolyA, length($seq)),",$endOfPolyA,", length($seq)-$endOfPolyA, "\n" if ($verboseHMM); 
	return (length($seq)-$endOfPolyA, substr($seq, $endOfPolyA, length($seq)));
}


##############################################
# compute score of a dynamic programming matrix cell
##############################################
sub computeScore {
	# this is [i,j] in the dynamic programming matrix
	my ($stateNumber, $seqPos, $seq) = @_;

	print "\nCompute MATRIX V for [$stateNumber][$seqPos]     character ", substr($seq, $seqPos, 1), "\n" if ($verboseHMM);
	my $max = -Inf;
	my $argMax = "";

	print "\tState $state2Name[$stateNumber]   (stateNumber $stateNumber)\n" if ($verboseHMM);
	my @prevStates = split(/,/, $state2IncomingTransitions[$stateNumber]); 
	my @transitionScores = split(/,/, $state2IncomingTransitionsScore[$stateNumber]); 
	print "\tprevStates @prevStates   transitionScores @transitionScores\n" if ($verboseHMM);

	for (my $k=0; $k<=$#prevStates; $k++) {
		my $prev = $prevStates[$k];
		
		my $emitScore = getEmissionScore($stateNumber, $seqPos, $seq);
		# previous position (where to access the matrix V[k,i-1]) 
		my $prevPos = $seqPos - 1;

		die "ERROR in computeScore($stateNumber, $seqPos): matrix cell V[$prev][$prevPos] is undefined\n" if (! defined $V[$prev][$prevPos]); 
		my $score = $emitScore + $transitionScores[$k] + $V[$prev][$prevPos];
		printf "\t    Emit + Trans (%s-->%s) + V[$state2Name[$prev]][$prevPos] =  %1.3f + %1.3f + %1.3f = %1.3f  %s\n", ,$state2Name[$prev],$state2Name[$stateNumber],$emitScore, $transitionScores[$k], $V[$prev][$prevPos], $score, ($score > $max ? " --> NEW MAX" : "") if ($verboseHMM);
		if ($score > $max) {
			$max = $score;
			$argMax = "$prev,$prevPos";
		}
	}	

	print "\tResult for V[$state2Name[$stateNumber]][$seqPos] = $max  ($argMax)\n" if ($verboseHMM);
	$V[$stateNumber][$seqPos] = $max;
	$T[$stateNumber][$seqPos] = $argMax;
	return;
}	

##############################################
# compute emission score 
##############################################
sub getEmissionScore {
	my ($stateNumber, $seqPos, $seq) = @_;

	# P states emits A's with a high prob
	if ($stateNumber == 0) {
		if (substr($seq, $seqPos, 1) eq "A") {
			return $logemitA;
		}else{	
			return $logemitNotA;
		}
	# S emits uniform
	}elsif ($stateNumber == 1) {
		return $logUniform;
	}

	die "ERROR: $stateNumber $state2Name[$stateNumber] is neither an 0 nor 1\n";
}



#########################################
# setup two-state HMM 
#########################################
sub generateHMM {
	# PolyA state
	my $stateNumber = 0;
	my $name="P";
	$state2Name[$stateNumber] = $name;
	$state2IncomingTransitions[$stateNumber] = 0;				# P -> P
	$state2IncomingTransitionsScore[$stateNumber] = log($P2P);

	# seq state
	$stateNumber ++;
	$name="S";
	$state2Name[$stateNumber] = $name;
	$state2IncomingTransitions[$stateNumber] .= 1 . ",";				# S -> S
	$state2IncomingTransitionsScore[$stateNumber] .= log(1) . ",";
	$state2IncomingTransitions[$stateNumber] .= 0 . ",";				# P -> S
	$state2IncomingTransitionsScore[$stateNumber] .= log(1-$P2P) . ",";
	
	return $stateNumber+1;
}


