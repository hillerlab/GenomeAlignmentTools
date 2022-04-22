#!/usr/bin/env perl
# Michael Hiller, 2013

# extracting syn/inv chains from a net file
# 	It reads the net file, keeps a list of all these chain IDs and then extracts them from the chain file
#
# output goes to stdout

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use POSIX;

$| = 1;		# == fflush(stdout)
my $verbose = 0;
my $keepSynNetsWithScore = INT_MAX;      # keep nets classified as syn if the are above this score threshold
my $keepInvNetsWithScore = INT_MAX;      # keep nets classified as inv if the are above this score threshold

my $usage = "$0 input.net input.chain output.chain  [-v -keepSynNetsWithScore int -keepInvNetsWithScore int   OTHERPARAMETERS]
\t-keepSynNetsWithScore          keep nets classified as syn if the are above this score threshold
\t-keepInvNetsWithScore          keep nets classified as inv if the are above this score threshold
\n
";
GetOptions ("v|verbose"  => \$verbose, 
	"keepSynNetsWithScore=i" => \$keepSynNetsWithScore, "keepInvNetsWithScore=i" => \$keepInvNetsWithScore)
	|| die "$usage\n";
die "$usage\n" if ($#ARGV < 2);

# read the net file
print "Read netfile $ARGV[0]\n";
if ($ARGV[0] =~ /.gz$/) {
	open(file1, "zcat $ARGV[0]|") || die "ERROR: cannot open $ARGV[0]\n";
}else{	
	open(file1, $ARGV[0]) || die "ERROR: cannot open $ARGV[0]\n";
}
my @file = <file1>;
chomp(@file);
close file1;

if ($file[0] !~ /^net /) {
	die "ERROR: expect file to start with net\n";
}

# tmp file in shm to keep the chainIDs
my $tmpChainIDFile = `set -o pipefail; mktemp /dev/shm/$ENV{'USER'}.chainIDs.XXXXXXX`;
chomp($tmpChainIDFile);
print "$tmpChainIDFile\n" if ($verbose);
open(CHAINFILE, ">$tmpChainIDFile") || die "ERROR: cannot write to $tmpChainIDFile\n";
my $numChainsKept = 0;	# for the stats

# start reading at line 1 (after the first net)
for (my $i=1; $i<=$#file; $i++) {
	my $line1 = $file[$i];

	next if ($line1 =~ / gap /);		# gap lines are handled when their corresponding net is handled
	if ($line1 =~ /^net /) {
		next;
	}

	# determine the level of the current line = # spaces at the line beginning 	
	if ($line1 !~ /^([ ]+)([fill|gap].*)/) {
		print STDERR "ERROR: expect fill or gap in $line1\n";
		exit -1;
	}
	my $lineRest = $2;
	print "current line: $lineRest\n" if ($verbose);

	# fill line --> check if we keep that net
	if ($line1 =~ / fill /) {
		# fill 61168697 78 chrIV + 24270860 75 id 5545645 score 3359 ali 75 qDup 75 type nonSyn tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
		my @f = split(/ /,$lineRest);
		# extract the type. This can be a variable position
		# fill 24258609 108 chr11 + 76807460 108 id 1579984 score 3220 ali 108 qOver 0 qFar 43429601 qDup 0 type inv tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
		# fill 24260574 62 chr3 - 97726917 62 id 1686589 score 3151 ali 62 qDup 0 type nonSyn tN 0 qN 0 tR 0 qR 0 tTrf 0 qTrf 0
		my $type = "";
		if ($lineRest !~ /type (\w+) /) {
			print STDERR "ERROR: cannot parse type from this fill line: $lineRest\n";
		}else{
			$type = $1;
		}
		
		# here is the filter: pass score, Tsize and Qsize to the function
		if ( passesFilter($f[10], $f[2], $f[6], $type) eq "true" ){
			my $chainID = getChainID($lineRest);
			print CHAINFILE "$chainID\n";
			print "KEEP: chain ID $chainID: $line1\n" if ($verbose);
			$numChainsKept++;
		}else{
			print "SKIP: $line1\n" if ($verbose);
		}
	}
}
close file1;
close(CHAINFILE);

print "DONE\nKeep $numChainsKept from $ARGV[0]\n";
print "Extract these chains from $ARGV[1] into $ARGV[2]\n";

my $call = "set -o pipefail; chainFilter -idFile=$tmpChainIDFile $ARGV[1] > $ARGV[2]\n";
print "$call\n";
system("$call") == 0 || die "ERROR: chainPreNet command failed\ncall: $call\n";
print "DONE\n";

# cleanup
`set -o pipefail; rm -f $tmpChainIDFile`;


######################################################
# extract chain ID from the fill line
######################################################
sub getChainID {
	my $line = shift;
	my $chainID = -1;
	if ($line !~ /id (\d+) score/) {
		cleanDie("ERROR: cannot parse 'id' from this fill line: $line\n");
	}else{
		$chainID = $1;
	}
	return $chainID;
}

######################################################
# batch or individual filtering
######################################################
sub passesFilter {
	my ($score, $Tsize, $Qsize, $type) = @_;
	
	# first check if syn or inv filter applies
	if ($type eq "syn" && $score >= $keepSynNetsWithScore) {
		return "true";
	}
	if ($type eq "inv" && $score >= $keepInvNetsWithScore) {
		return "true";
	}

	return "false";
}


######################################################
# clean up tmp file and then die
######################################################
sub cleanDie {
	my $output = shift;
	`set -o pipefail; rm -f $tmpChainIDFile`;
	die "$output\n";
}
