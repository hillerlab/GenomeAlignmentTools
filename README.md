# Genome Alignment Tools

This repository provides tools for improving the sensitivity and specificity of pairwise genome alignments [1,2]. These tools make use of [lastz](http://www.bx.psu.edu/~rsharris/lastz/) and the alignment chain and net concept [3]. 

These tools integrate into the standard lastz/chain/net workflow of genome alignment as follows:
1. genome-wide local alignments with lastz
2. building alignment chains
3. NEW: highly-sensitive local alignments to improve alignment chains by detecting remote homologies
4. NEW: chainCleaner to improve the specificity in alignment chains
5. NEW: chainNet with parameter -rescore to compute exact scores of all nets
6. NEW: non-nested net filtering to keep nets that likely represent orthologous alignments


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
Since both tools need the Kent (UCSC) source code, you need to edit the KENTSRC_DIR variable in src/makefile and set the path to the kent source code (should end with kent/src/).

Afterwards, compile chainCleaner and chainNet
```
cd src
make 
```

The subdirectory precompiledBinary_x86_64/ already contains Linux x86_64 compiled binaries.



# Highly-sensitive local alignments
patchChain.perl perform a highly sensitive local pairwise alignment for loci flanked by aligning blocks [1]. 
Given an alignment chain [3], it considers all chains that pass the score and span filters (optional parameters), extracts all the unaligning loci and creates local alignment jobs.

__Usage:__
```
patchChain.perl chain.infile  2bit_target_fullPath  2bit_query_fullPath  \
  chrom.sizes_target_fullPath  chrom.sizes_query_fullPath   [options]
```
Call patchChain.perl without any parameters to see the full parameter list.

__Example:__
Here, we use the parameters for highly sensitive local alignments from [1] and the HoxD55 scoring matrix (provided by the kent source code (subdirectory src/blatz/HoxD55.q); also contained in the example directory).
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
```

# chainCleaner
chainCleaner improves the specificity in genome alignment chains by detecting and removing local alignments that obscure the evolutionary history of genomic rearrangements [2].
The input is a chain file, ideally after adding alignments found with highly sensitive parameters if distal species are compared.
The output is a chain file that contains re-scored and score-sorted chains after removing the local alignments from the parent chains and adding them as individual chains. 
The resulting output file can be used to get alignment nets by running chainNet [3].


__Usage:__

```
chainCleaner in.chain tNibDir qNibDir out.chain out.bed -net=in.net 
 OR 
   chainCleaner in.chain tNibDir qNibDir out.chain out.bed -tSizes=/dir/to/target/chrom.sizes -qSizes=/dir/to/query/chrom.sizes 
 First option:   you have netted the chains and specify the net file via -net=netFile
 Second option:  you have not netted the chains. Then chainCleaner will net them. In this case, you must specify the chrom.sizes file for the target and query with -tSizes/-qSizes
 tNibDir/qNibDir are either directories with nib files, or the name of a .2bit file

output:
   out.chain      output file in chain format containing the untouched chains, the original broken chain and the modified breaking chains. NOTE: this file is chainSort-ed.
   out.bed        output file in bed format containing the coords and information about the removed chain-breaking alignments.
```

Call chainCleaner without any parameters to see the full parameter list.

__Example run of chainCleaner:__
Before executing, download the human and mouse genome.
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit
mv hg38.2bit mm10.2bit example/
```

Run chainCleaner on the chains provided in example/
``` 
chainCleaner example/hg38.mm10.chr1.chain.gz -tSizes=example/hg38.chrom.sizes -qSizes=example/mm10.chrom.sizes example/hg38.2bit example/mm10.2bit example/hg38.mm10.chr1.cleaned.chain example/removedSuspects.bed -linearGap=loose 


Verbosity level: 1
foldThreshold: 0.000000    LRfoldThreshold: 2.500000   maxSuspectBases: 2147483647  maxSuspectScore: 100000  minBrokenChainScore: 50000  minLRGapSize: 0
0. need to net the input chains example/hg38.mm10.chr1.chain.gz (no net file given) ...
		tempfile for netting: tmp.chainCleaner.XLJ2Vxs.net
Got 455 chroms in example/hg38.chrom.sizes, 66 in example/mm10.chrom.sizes
Finishing nets
writing stdout
writing /dev/null
DONE (nets in tmp.chainCleaner.XLJ2Vxs.net)
1. parsing fills/gaps from tmp.chainCleaner.XLJ2Vxs.net and getting valid breaks ...
1.1 read net file tmp.chainCleaner.XLJ2Vxs.net into memory ...
DONE
1.2 get fills/gaps from tmp.chainCleaner.XLJ2Vxs.net ...
DONE
1.3 get aligning regions from tmp.chainCleaner.XLJ2Vxs.net ...
DONE
1.4 get valid breaks ...
DONE
Remove temporary netfile tmp.chainCleaner.XLJ2Vxs.net
DONE (parsing fills/gaps and getting valid breaks)
2. reading breaking and broken chains from example/hg38.mm10.chr1.chain.gz and write irrelevant chains to example/hg38.mm10.chr1.cleaned.chain.unsorted ...
DONE
3. reading target and query DNA sequences for breaking and broken chains ...
DONE
4. loop over all breaks. Remove suspects if they pass our filters and write out deleted suspects to example/removedSuspects.bed ...
DONE
5. write the (new) breaking and the broken chains to example/hg38.mm10.chr1.cleaned.chain.unsorted ...
DONE
6. chainSort example/hg38.mm10.chr1.cleaned.chain.unsorted example/hg38.mm10.chr1.cleaned.chain ...
DONE
7. free memory ...
DONE
memory usage 5562015744, utime 1058 s/100, stime 256

ALL DONE. New chains are in example/hg38.mm10.chr1.cleaned.chain. Deleted suspects in example/removedSuspects.bed
```


# chainNet
Given a set of alignment chains, chainNet produces alignment nets, which is a hierarchical collection of chains or parts of chains that attempt to capture only orthologous alignments [3]. 
The original chainNet the score of "sub-net" (nets that come from a part of a chain and fill a gap in a higher-level net) by the fraction of aligning bases. This can lead to a bias in case the aligning blocks of a chain are not equally distributed. We implemented a new parameter "-rescore" in chainNet that computes the real score of each subnet [2].

__Usage:__
```
chainNet - Make alignment nets out of chains
usage:
   chainNet in.chain target.sizes query.sizes target.net query.net
where:
   in.chain is the chain file sorted by score
   target.sizes contains the size of the target sequences
   query.sizes contains the size of the query sequences
   target.net is the output over the target genome
   query.net is the output over the query genome
options:
   -rescore                    compute the real score of the sub-net (instead of approximating it based on the fraction of aligning bases in the subnet)
                               The real score will be much more precise especially for imbalanced chains where most aligning blocks are on one side.
                               This flag will set minScore=0. Each subnet with a negative score gets score 1. Afterwards, run a non-nested score filter.
                               Note: Rescoring is only implemented for the target species net.
                               With this flag, you need to give the target and query genome sequence (-tNibDir and -qNibDir) and specify -linearGap
   -tNibDir=fileName           target genome file (2bit or nib format)
   -qNibDir=fileName           query genome file (2bit or nib format)
```
Call chainNet without any parameters to see the full parameter list.

# Non-nested net filtering
Before building a multiple alignment from the pairwise alignment nets, it is recommended to remove low-scoring alignment nets that are unlikely to represent real homologies.
While the netFilter program [3] removes nets that do not satisfy the specified score and size criteria including all nested nets, 
NetFilterNonNested.perl applies a non-nested filtering procedure that considers and filters each net individually [1]. 
This avoids removing nested nets that would satisfy the specified criteria, even if a parent net is removed. 

__Usage:__
In [1], we applied the following filter criteria:

* UCSC "syntenic net" criteria (thresholds: minTopScore=300000, minSynScore=200000, minSynSize=20000, minSynAli=10000, maxFar=200000) were applied to all placental mammal alignments that have well-assembled genomes 
by running 
`
NetFilterNonNested.perl -doUCSCSynFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 ref.query.net.gz > ref.query.filtered.net
`
* keeping nets that score higher than 100000 and keeping all nested nets that align to the same locus (inversions or local translocations) if they score higher than 5000 for placental mammals with less-well assembled genomes 
`
NetFilterNonNested.perl -doScoreFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 -minScore1 100000  ref.query.net.gz > ref.query.filtered.net
`
* keeping nets that score higher than 10000 and keeping all nested nets that align to the same locus if they score higher than 3000 for non-placental mammals 
`
NetFilterNonNested.perl -doScoreFilter -keepSynNetsWithScore 3000 -keepInvNetsWithScore 3000 -minScore1 10000 ref.query.net.gz > ref.query.filtered.net
`

Call NetFilterNonNested.perl without any parameters to see all filtering options.
