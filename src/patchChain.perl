#!/usr/bin/perl

# Michael Hiller, MPI-CBG & MPI-PKS, Dresden, Germany

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use List::Util 'shuffle';

$| = 1;		# == fflush(stdout)
my $verbose = 0;
my $chainMinScore = 0; 	  		# consider only chains with a chainMinScore, min size on query and/or target; default consider all
my $chainMinSizeT = 0;
my $chainMinSizeQ = 0;	
my $minIdentity = 0;				# optional parameters for axt filtering: min %id and or min entropy in a given windowSize (windowSize must always be given if identity or entropy is used)
my $minEntropy = 0;
my $windowSize = 0;
my $gapMinSizeT = 10;			# patch only gaps that are at least that long on the target side
my $gapMinSizeQ = 10;			# patch only gaps that are at least that long on the query sid
my $gapMaxSizeT = 100000;	  	# patch only gaps that are at most that long on the target sid
my $gapMaxSizeQ = 100000;	  	# patch only gaps that are at most that long on the query side
my $lastzParameters = "--format=axt Q=/path/to/HoxD55.q K=1500 L=2500 M=0 T=0 W=5";		# specify as a string --> will be passed as is to lastz
my $unmask = 0;				# initially a flag for unmasking (lower case to upper case) characters from the twobit files. This variable will be replaced by '[unmask]' if the user sets this flag. 
my $numJobs = 1000;				# number of jobs that are created; each job produces one input and one output file
my $jobDir = "doPatchChain";		# directory where job csh files are stored
my $outputDir = "doPatchChain";		# directory where the output files are written to. NOTE: this is an abs path or must be relative to the $jobDir
my $jobList = "jobList"; 		# contains the patch jobs


# options
my $usage = "$0 creates jobs that perform a highly sensitive local pairwise alignment for loci flanked by aligning blocks from a chain.
It considers all chains that pass the score and span filters (parameters below), extracts all the gaps that have a min/max span on the target and query and creates local alignment jobs. 
The output are $numJobs .csh scripts that each call lastz, filter the output for \%id and entropy (if parameters given) and produce one output file per job.
The order of the individual alignments are shuffled to create balanced cluster jobs that should result in similar runtimes. 
The resulting .csh jobs could be executed in parallel on a compute cluster. 

Requirement: axtToPsl from the kent source code and lastz must be in \$PATH

usage: $0 chain.infile  2bit_target_fullPath  2bit_query_fullPath  chrom.sizes_target_fullPath  chrom.sizes_query_fullPath   [options]

Options are: 
         -chainMinScore [0]         consider only chains above that score threshold
         -chainMinSizeT [0]         consider only chains above that span chainMinSizeT bp in the target (reference) genome
         -chainMinSizeQ [0]         consider only chains above that span chainMinSizeQ bp in the query (other) genome
         -minIdentity [0]           for filtering: exclude all lastz hits with less than minIdentity (percent, int) identical bases in a window of windowSize (-windowSize needs to be given)
         -minEntropy [0]            for filtering: exclude all lastz hits with an entropy less than minEntropy (float, [0-2]) in a window of windowSize (-windowSize needs to be given)	
         -windowSize [0]            for filtering: window size	
         -gapMinSizeT [10]          patch only gaps that are at least that long on the target side
         -gapMinSizeQ [10]          patch only gaps that are at least that long on the query side
         -gapMaxSizeT [100000]      patch only gaps that are at most that long on the target side
         -gapMaxSizeQ [100000]      patch only gaps that are at most that long on the query side
         -lastzParameters [$lastzParameters]  
                                    use \"..\" to specify all parameters that are exactly passed to lastz
         -unmask                    flag for unmasking (lower case to upper case) characters from the 2bit files
         -numJobs [1000]            number of generated cluster jobs (.csh files); each job produces one input and one output file
         -jobList [$jobList]        joblist file (will be clobbered if exists)
         -jobDir [$jobDir]          dir for the job csh files; will be created if not exist; existing files will be clobbered
         -outputDir [$outputDir]    dir for output files; This dir has be created outside of $0 . NOTE: this is an abs path or must be relative to the $jobDir
			
         -v                         detailed verbose output
\n";
GetOptions ("v|verbose"  => \$verbose, 
		"chainMinScore=i" => \$chainMinScore, "chainMinSizeT=i" => \$chainMinSizeT, "chainMinSizeQ=i" => \$chainMinSizeQ, 
		"minIdentity=i" => \$minIdentity, "minEntropy=f" => \$minEntropy, "windowSize=i" => \$windowSize, 
		"gapMinSizeT=i" => \$gapMinSizeT, "gapMinSizeQ=i" => \$gapMinSizeQ, "gapMaxSizeT=i" => \$gapMaxSizeT, "gapMaxSizeQ=i" => \$gapMaxSizeQ, 
		"lastzParameters=s" => \$lastzParameters, "unmask" => \$unmask,
		"jobDir=s" => \$jobDir, "outputDir=s" => \$outputDir, "jobList=s" => \$jobList, "numJobs=i" => \$numJobs

) ||	die "$usage\n";
die "$usage\n" if ($#ARGV < 4);

# two bit + chrom.sizes files
my $T2bit = $ARGV[1];					# full path to 2bit file for query (the other genome)
my $Q2bit = $ARGV[2];					# full path to 2bit file for target (the reference genome)
my $tChromSizes = $ARGV[3];			# needed for axtToPsl
my $qChromSizes = $ARGV[4];

# if entropy or %id given, windowSize cannot be 0 --> check
if (($minEntropy != 0 || $minIdentity != 0) && ($windowSize == 0)) {
	die "ERROR: minEntropy or minIdentity given but windowSize is 0\n";
}

# replace the content of the unmask by '[unmask]' if the user sets this flag, otherwise ""
if ($unmask) {
    $unmask = '\[unmask\]';
} else {
    $unmask = "";
}

# extract all blocks that need to be patched
my @regionsToBePatched;
# read all chains
open(file1, $ARGV[0]) || die "ERROR: cannot open $ARGV[0]\n";
my $line1;
while ($line1 = <file1>) {
	chomp($line1);
	next if ($line1 =~ /^#/);	# ignore comment

	# read the chain line 
	# e.g. chain 196228 chr4 62094675 + 12690854 12816143 chr23 24050845 - 20051667 20145391 1252
	if ($line1 =~ /^chain /) {
		my @c = split(/ /, $line1);
		my $ID = $c[12];

		my $score 		= $c[1];
		my $tName 		= $c[2];
		my $tStart 		= $c[5]; 
		my $tEnd 		= $c[6];
		my $qName 		= $c[7];
		my $qStart 		= $c[10];
		my $qEnd 		= $c[11];
		my $qStrand      	= $c[9];
		my $qChromSize          = $c[8];
		print STDERR "ERROR: target strand is not + for chain: $line1\n" if ($c[4] ne "+");

		print "read $line1\n" if ($verbose);
		# check if we consider this chain
		if ($score >= $chainMinScore && ($tEnd - $tStart >= $chainMinSizeT) && ($qEnd - $qStart >= $chainMinSizeQ)) {
			print "VALID chain: score $score  target span: ", $tEnd - $tStart, "  query span: ", $qEnd - $qStart, "\n" if ($verbose);
		}else{
			print "invalid chain: score $score  OR  target span: ", $tEnd - $tStart, "  OR  query span: ", $qEnd - $qStart, "  less than given thresholds\n" if ($verbose);
			# read the rest of the chain without any action
			while ($line1 = <file1>) {
				last if ($line1 eq "\n");
			}
			next;		# read next chain
		}

		my $curTPos = $tStart;
		my $curQPos = $qStart;
		while ($line1 = <file1>) {
			last if ($line1 eq "\n");
			chomp ($line1);
			my @a = split(/\t/, $line1);
			if ($#a == 0) {		# last block
				my $TblockEnd = $curTPos + $a[0];
				my $QblockEnd = $curQPos + $a[0];
				print "last\n\tT: block: $tName $curTPos $TblockEnd  		Q: block: $qName $curQPos $QblockEnd   ($line1)\n" if ($verbose);
			}else{
				my $blockLen = $a[0];
				my $TblockEnd = $curTPos + $blockLen;
				my $TgapEnd   = $curTPos + $blockLen + $a[1];
				my $QblockEnd = $curQPos + $blockLen;
				my $QgapEnd   = $curQPos + $blockLen + $a[2];
				my $qGapSpan = $QgapEnd - $QblockEnd;
				my $tGapSpan = $TgapEnd - $TblockEnd;
				print "\tT: block: $tName $curTPos $TblockEnd   downstream gap: $tName $TblockEnd $TgapEnd  (span $tGapSpan)     Q: block: $qName $curQPos $QblockEnd  downstream gap: $qName $QblockEnd $QgapEnd    (span $qGapSpan)   ($line1)\n" if ($verbose);
				
				# check if should try to patch this gap
				if ($tGapSpan >= $gapMinSizeT  &&  $tGapSpan <= $gapMaxSizeT  &&  $qGapSpan >= $gapMinSizeQ  &&  $qGapSpan <= $gapMaxSizeQ) {
					print "\t--> will be patched\n" if ($verbose);
					# lastz works with coords [a..b] where a is the first and b is the last base (not the position downstream). Everything is 1-based --> add +1 to the gap starts
					$TblockEnd++;
					$QblockEnd++;

					# calculate +strand coordinates if neccessary
					if( $qStrand eq "-" ) { 
					    my $x=$qChromSize-$QgapEnd+1;
					    my $y=$qChromSize-$QblockEnd+1;
					    push @regionsToBePatched, "$tName\t$TblockEnd\t$TgapEnd\t$qName\t$x\t$y";

					} else { 
					    push @regionsToBePatched, "$tName\t$TblockEnd\t$TgapEnd\t$qName\t$QblockEnd\t$QgapEnd";
					}
				}else{
				    print "\t--> will NOT be patched (outside of min-max spans for T or Q)\n" if ($verbose);
				}
				$curQPos = $QgapEnd;
				$curTPos = $TgapEnd;
			}
		}
	}else{
		die "ERROR: expect chain line but got $line1\n";
	}
}
close file1;
print "\n\nDone with reading gaps\n" if ($verbose);


###############################
# create the jobs
my $numToBePatched = $#regionsToBePatched + 1;
my $numPerFile = int($numToBePatched / $numJobs);
$numPerFile = 1 if ($numPerFile < 1);
print "num gaps ToBePatched: $numToBePatched  numJobs: $numJobs  --> $numPerFile gaps per job\n" if ($verbose);
`mkdir -p $jobDir`;



###############################
# shuffle the regions --> this should give job files having a similar runtime
my @shuffled = shuffle(@regionsToBePatched);

###############################
# write $numJobs files in $jobDir each containing $numPerFile coordinate (target, query) pairs
`mkdir -p $outputDir`;
my $curJob = 0;
# create first file
open fileout1, ">$jobDir/job$curJob.csh" || die "ERROR: cannot write job file $jobDir/job$curJob\n";
`chmod +x $jobDir/job$curJob.csh`;   # make executable
print fileout1 "#!/bin/csh -efx\nrm -f $outputDir/patch$curJob.psl $outputDir/patch$curJob.axt\n";
# loop through all to-be-patched regions
for (my $i=0; $i<=$#shuffled; $i++) {
	if ($i % $numPerFile == 0 && $i != 0) {		# we have opened the file for the first line ($i==0) already
		print "line $i --> new file $curJob\n" if ($verbose);
		# add filtering only if either minEntropy or minIdentity is given
		if ($minEntropy != 0 || $minIdentity != 0) {
			my $command .= sprintf "filterAxtIdentityEntropy.py $outputDir/patch$curJob.axt $minIdentity $minEntropy $windowSize $outputDir/patch$curJob.filtered.axt\n";
			$command .= sprintf "axtToPsl $outputDir/patch$curJob.filtered.axt $tChromSizes $qChromSizes $outputDir/patch$curJob.filtered.psl\n";
			$command .= sprintf "rm $outputDir/patch$curJob.filtered.axt $outputDir/patch$curJob.axt\n";
			print "$command" if ($verbose);
			print fileout1 "$command\n";
		}else{
			my $command .= sprintf "axtToPsl $outputDir/patch$curJob.axt $tChromSizes $qChromSizes $outputDir/patch$curJob.psl\n";
			$command .= sprintf "rm $outputDir/patch$curJob.axt\n";
			print "$command" if ($verbose);
			print fileout1 "$command\n";
		}
		close fileout1;
		$curJob ++;
		# create new file
		open fileout1, ">$jobDir/job$curJob.csh" || die "ERROR: cannot write job file $jobDir/job$curJob.csh\n";
		`chmod +x $jobDir/job$curJob.csh`;   # make executable
		print fileout1 "#!/bin/csh -efx\nrm -f $outputDir/patch$curJob.psl $outputDir/patch$curJob.axt\n";		# always delete the previous output file in case a job is re-run
	}	

	# get the command: first get T and Q coords
	my ($tChr, $tStart, $tEnd, $qChr, $qStart, $qEnd) = (split(/\t/, $shuffled[$i]))[0,1,2,3,4,5];
	# print the command for patching this region
	# first lastz (output axt), then Saatviks python script that outputs all full axt's that fulfill the minID & minEntropy requirement, then axtToPsl (keeps the genome coords)
	# add output to the $outputDir/patch$curJob.axt file
	# for the csh (NOT BASH) we need to mask the []
	my $command = sprintf "lastz $T2bit/$tChr\\\[$tStart..$tEnd\\\]$unmask $Q2bit/$qChr\\\[$qStart..$qEnd\\\]$unmask  $lastzParameters >> $outputDir/patch$curJob.axt";
    
	
	print "$command\n" if ($verbose);
	print fileout1 "$command\n";
}
# add filtering only if either minEntropy or minIdentity is given
if ($minEntropy != 0 || $minIdentity != 0) {
	my $command .= sprintf "filterAxtIdentityEntropy.py $outputDir/patch$curJob.axt $minIdentity $minEntropy $windowSize $outputDir/patch$curJob.filtered.axt\n";
	$command .= sprintf "axtToPsl $outputDir/patch$curJob.filtered.axt $tChromSizes $qChromSizes $outputDir/patch$curJob.filtered.psl\n";
	$command .= sprintf "rm $outputDir/patch$curJob.filtered.axt $outputDir/patch$curJob.axt\n";
	print "$command" if ($verbose);
	print fileout1 "$command\n";
}else{
	my $command .= sprintf "axtToPsl $outputDir/patch$curJob.axt $tChromSizes $qChromSizes $outputDir/patch$curJob.psl\n";
	$command .= sprintf "rm $outputDir/patch$curJob.axt\n";
	print "$command" if ($verbose);
	print fileout1 "$command\n";
}

close fileout1;


###############################
# create the jobList; each job is one of the csh scripts we produced above
open fileout1, ">$jobList" || die "ERROR: cannot write jobList file $jobList\n";
opendir(DIR, $jobDir) or die "can't read from $jobDir: $!";
foreach my $file (sort readdir(DIR)) {
	 if ($file =~ /^job.*.csh/) {
		 print fileout1 "$jobDir/$file\n";
	 }
}
closedir(DIR);
close fileout1;


###############################
print "\nDONE\nnum gaps to-be-patched: $numToBePatched   numJobs: $numJobs   jobDir: $jobDir  outputDir: $outputDir  file listing all jobs: $jobList\n";
print "$jobList contains independent jobs that could be executed in parallel on a compute cluster. \n";
