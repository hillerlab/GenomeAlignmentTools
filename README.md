# Genome Alignment Tools

Tools for improving the sensitivity and specificity of genome alignments


# Requirements
* install [lastz](http://www.bx.psu.edu/~rsharris/lastz/)
* install the Kent (UCSC) source code by doing
```
# go to whereever you want to have the Kent source code
cd /path/to/
# git clone and make the libraries
git clone http://genome-source.cse.ucsc.edu/kent.git
cd kent/src
make USE_SAMTABIX=0 USE_TABIX=0 USE_BAM=0 libs
```
These binaries should be located in a directory that is added to your $PATH environment variable. 

# Installation
The perl and python scripts just need to be copied to a directory that is contained in your $PATH environment variable. 
Only chainCleaner and the modified chainNet (with the -rescore option) need to be compiled. 
Since both needs the Kent (UCSC) source code, you need to edit the KENTSRC_DIR variable in src/makefile and set the path to the kent source code (should end with kent/src/).

Afterwards, compile chainCleaner and chainNet
```
cd src
make 
```

The subdirectory precompiledBinary_x86_64/ already contains Linux x86_64 compiled binaries.



# Highly-sensitive local alignments
patchChain.perl perform a highly sensitive local pairwise alignment for loci flanked by aligning blocks. 
Given an alignment chain [3], it considers all chains that pass the score and span filters (optional parameters), extracts all the unaligning loci and creates local alignment jobs.

__Usage:__
```
patchChain.perl chain.infile  2bit_target_fullPath  2bit_query_fullPath  \
  chrom.sizes_target_fullPath  chrom.sizes_query_fullPath   [options]
```
Call patchChain.perl without any parameters to see the full parameter list.

__Example:__
Here, we use the parameters for highly sensitive local alignments from [1] "Q=/path/to/HoxD55.q K=1500 L=2500 M=0 T=0 W=5"
The HoxD55 scoring matrix is provided by the kent source code (subdirectory src/blatz/HoxD55.q) and also copied to the example directory.
```
patchChain.perl example/hg38.danRer10.chain example/hg38.2bit example/danRer10.2bit example/hg38.chrom.sizes example/danRer10.chrom.sizes \
    -chainMinScore 5000 -gapMaxSizeT 500000 -gapMaxSizeQ 500000 -gapMinSizeT 30 -gapMinSizeQ 30 \
	 -numJobs 10 -jobDir jobs -jobList jobList -outputDir pslOutput \
	 -minEntropy 1.8 -windowSize 30 -minIdentity 60 \
    -lastzParameters "--format=axt K=1500 L=2500 M=0 T=0 W=5 Q=example/HoxD55.q"

# this results in 10 alignment jobs that are located in jobs/ and listed in 'jobList'

# before executing, download the human and zebrafish genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
wget http://hgdownload.cse.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.2bit
mv danRer10.2bit hg38.2bit example/

# now execute these jobs on a compute cluster or run them sequentially by doing 'chmod +x jobList; ./jobList'

# concatenate all new results
find pslOutput -name "*.psl" | xargs -i cat {} > newAlignments.psl

# combine the genome-wide lastz results (the combined psl file that was used to create the input chains) and the newly found psl alignments
cat genomeWide.lastz.psl newAlignments.psl > all.psl 

# use axtChain from the Kent source to compute alignment chains that include the new alignments

