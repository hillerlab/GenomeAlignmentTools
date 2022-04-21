#!/usr/bin/env perl

# Michael Hiller, 2013
# for the new genome alignment pipeline of psl filtering, chain patching, chain and net synteny filtering 
# This script performs the chain filtering, then nets the filtered chains and filters the nets again. 
# 2020: Can now be run on delta (our compute cluster) by running the netClass step (this requires a SQL database on a genome browser server) via ssh on genome (our genome browser server)
# The genome browser server must be ssh-able without password and must be specified in the variable $dbHost below
#
use strict;
use warnings;
use diagnostics;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use POSIX;
use Sys::Hostname;

$| = 1;		# == fflush(stdout)
my $verbose = 0;
my $inChain = "";				# filename + path
my $inNet = "";				# filename + path
my $outFilteredChain = "";	# output filename + path
my $outFilteredNet = "";	# output filename + path
my $Tassembly="";				# ref and query assembly
my $Qassembly="";
my $minScore="";				# can be an array, comma-separated
my $minTsize="";
my $minQsize="";
my $keepSynNetsWithScore = INT_MAX;      # keep nets classified as syn if the are above this score threshold
my $keepInvNetsWithScore = INT_MAX;      # keep nets classified as inv if the are above this score threshold

my $genomePath = $ENV{'genomePath'};       # path to the directory that contains the gbdb-HL/$assembly directories
my $dbHost = "genome.tbg.senckenberg.de";  # must be adjusted if run outside of the hillerlab setting

my $usage = "usage: $0 inChainFile inNetFile outFilteredChain outFilteredNet ReferenceAssembly QueryAssembly minScores minTsizes minQsizes [-v|verbose -keepSynNetsWithScore int, keepInvNetsWithScore int] \n
	inChain, inNet, outFilteredChain and outFilteredNet should just be filenames including the path to the file
	inChain and inNet can be gzipped (DO NOT UNZIP)
	the to-be-produced files outFilteredChain and outFilteredNet must NOT have the .gz suffix ($0 will gzip them afterwards)
	minScores minTsizes minQsizes  can be a comma-separated list
	-keepSynNetsWithScore          keep nets classified as syn if the are above this score threshold
	-keepInvNetsWithScore          keep nets classified as inv if the are above this score threshold

	$0 can be run on delta and on genome
\n";

# options
GetOptions ("v|verbose" => \$verbose, "keepSynNetsWithScore=i" => \$keepSynNetsWithScore, "keepInvNetsWithScore=i" => \$keepInvNetsWithScore) || die "$usage";	
#die "ERROR: you must be on genome to execute $0\n" if (hostname ne "genome");

if ($#ARGV < 7 ) {
	die "$usage";
}
$inChain = $ARGV[0];
$inNet = $ARGV[1];
$outFilteredChain = $ARGV[2];
$outFilteredNet = $ARGV[3];
$Tassembly = $ARGV[4];
$Qassembly = $ARGV[5];
$minScore = $ARGV[6];
$minTsize = $ARGV[7];
$minQsize = $ARGV[8];

my $tmpDir = `set -o pipefail; mktemp -d`;
chomp($tmpDir);
print "set tmpDir to $tmpDir\n";

# split the minscores and sizes
my @minScores = split(/,/, $minScore);
my @minTsizes = split(/,/, $minTsize);
my @minQsizes = split(/,/, $minQsize);
die "ERROR: number of minScores differ from minTsizes\n" if ($#minScores != $#minTsizes);
die "ERROR: number of minScores differ from minQsizes\n" if ($#minScores != $#minQsizes);

my $call="";
# filter the chains
for (my $i=0; $i<=$#minScores; $i++) {
	$call="set -o pipefail; chainFilter $inChain -notQ=chrM -notT=chrM -minScore=$minScores[$i] -qMinSize=$minQsizes[$i] -tMinSize=$minTsizes[$i] > $tmpDir/filtered$i.chain";
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: chainFilter failed\ncall: $call\n";
}

# now extract the chain IDs that are syn/inv according to the minScore from the net
# NOTE: these chains might be below the score/span threshold applied to the chains, so they would be filtered out otherwise
if ($keepSynNetsWithScore < INT_MAX || $keepInvNetsWithScore < INT_MAX) {
	my $call ="set -o pipefail; ExtractSynInvChainsFromNet.perl $inNet $inChain $tmpDir/filteredSynInv.chain -keepSynNetsWithScore $keepSynNetsWithScore -keepInvNetsWithScore $keepInvNetsWithScore";
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: ExtractSynInvChainsFromNet.perl command failed\ncall: $call\n";
}

# now cat all the files and chain
$call = "set -o pipefail; cat $tmpDir/filtered*.chain | chainSort stdin stdout | 
   chainPreNet stdin $genomePath/gbdb-HL/$Tassembly/chrom.sizes $genomePath/gbdb-HL/$Qassembly/chrom.sizes $outFilteredChain";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: chainPreNet command failed\ncall: $call\n";

# net, netSyntenic and netClass
$call ="set -o pipefail; chainNet $outFilteredChain $genomePath/gbdb-HL/$Tassembly/chrom.sizes $genomePath/gbdb-HL/$Qassembly/chrom.sizes stdout /dev/null -minSpace=1 -rescore -linearGap=loose -tNibDir=$genomePath/gbdb-HL/$Tassembly/$Tassembly.2bit -qNibDir=$genomePath/gbdb-HL/$Qassembly/$Qassembly.2bit |
   netSyntenic stdin stdout | netClass stdin $Tassembly $Qassembly $tmpDir/filtered.net -noAr";
# execute as is if on genome
if (hostname eq "genome") {
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: chainNet failed. No $tmpDir/filtered.net is produced\ncall: $call\n";
}else{
# if not, run the netClass step on genome in a new tempdir via ssh
	$call ="set -o pipefail; chainNet $outFilteredChain $genomePath/gbdb-HL/$Tassembly/chrom.sizes $genomePath/gbdb-HL/$Qassembly/chrom.sizes stdout /dev/null -minSpace=1 -rescore -linearGap=loose -tNibDir=$genomePath/gbdb-HL/$Tassembly/$Tassembly.2bit -qNibDir=$genomePath/gbdb-HL/$Qassembly/$Qassembly.2bit | netSyntenic stdin $tmpDir/filtered.net.beforeNetClass";
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: chainNet failed. No $tmpDir/filtered.net.beforeNetClass is produced\ncall: $call\n";

	# now run netClass on genome
	$call = "set -o pipefail; ssh -x -o 'StrictHostKeyChecking = no' -o 'BatchMode = yes' $dbHost nice mktemp -d";
	print "$call\n" if ($verbose);
	my $remoteTempDir=`$call`;
	die "ERROR: $call failed\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
	chomp($remoteTempDir);
	print "created remote temp dir $remoteTempDir on $dbHost\n";

	# copy data
	$call = "set -o pipefail; scp $tmpDir/filtered.net.beforeNetClass $dbHost:$remoteTempDir";
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: $call failed\n";

	$call = "set -o pipefail; ssh -x -o 'StrictHostKeyChecking = no' -o 'BatchMode = yes' $dbHost netClass $remoteTempDir/filtered.net.beforeNetClass $Tassembly $Qassembly $remoteTempDir/filtered.net -noAr";
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: $call failed\n";

	$call = "set -o pipefail; scp $dbHost:$remoteTempDir/filtered.net $tmpDir/";
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: $call failed\n";

	$call = "set -o pipefail; ssh -x -o 'StrictHostKeyChecking = no' -o 'BatchMode = yes' $dbHost nice rm -rf $remoteTempDir";
	print "$call\n" if ($verbose);
	system("$call") == 0 || die "ERROR: $call failed\n";

}

# now filter the nets
# only filter for syn/inv nets if the parameters are given above
my $syninvFilter = "";
if ($keepSynNetsWithScore < INT_MAX || $keepInvNetsWithScore < INT_MAX) {
	$syninvFilter = "-keepSynNetsWithScore $keepSynNetsWithScore -keepInvNetsWithScore $keepInvNetsWithScore";
}
$call ="set -o pipefail; NetFilterNonNested.perl $tmpDir/filtered.net -minScore $minScore -minSizeT $minTsize -minSizeQ $minQsize $syninvFilter > $outFilteredNet";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: NetFilterNonNested.perl failed\ncall $call\n";
# gzip but remove any existing chain.gz and net.gz previous output files before
$call ="set -o pipefail; rm -f ${outFilteredChain}.gz ${outFilteredNet}.gz; gzip $outFilteredChain $outFilteredNet";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: call $call\n";



# cleanup
`set -o pipefail; rm -rf $tmpDir`;

print "DONE: result files: $outFilteredChain and $outFilteredNet\n";
#die "ERROR: $tmpDir\n";
