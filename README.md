**If you are only interested in generating alignment chains (e.g. as input for [TOGA](https://github.com/hillerlab/TOGA)), we recommend using our [lastz_make_chains](https://github.com/hillerlab/make_lastz_chains) pipeline that fully automates lastz, axtChain, RepeatFiller and chainCleaner.**

# Genome Alignment Tools

This repository provides tools for improving the sensitivity and specificity of pairwise genome alignments [1,2,3]. These tools make use of [lastz](http://www.bx.psu.edu/~rsharris/lastz/) and the alignment chain and net concept [4]. 

These tools integrate into the standard lastz/chain/net workflow of genome alignment as follows:
1. genome-wide local alignments with lastz
2. building alignment chains
3. NEW: highly-sensitive local alignments to improve alignment chains by detecting remote homologies
4. NEW: chainCleaner to improve the specificity in alignment chains
5. NEW: chainNet with parameter -rescore to compute exact scores of all nets
6. NEW: non-nested net filtering to keep nets that likely represent orthologous alignments


# Requirements
* install [lastz](http://www.bx.psu.edu/~rsharris/lastz/)
* install the light-weight Kent (UCSC) source code by doing
```
# go to whereever you want to have the GenomeAlignmentTools including Kent source code
cd /path/to/
# git clone and make the libraries
git clone https://github.com/hillerlab/GenomeAlignmentTools.git

cd GenomeAlignmentTools/kent/src

make

```
These binaries should be located in a directory that is added to your $PATH environment variable. 

If you are using BASH, for example, add following lines to your ~/.bashrc:

```
# Kent binaries
PATH=$PATH:/path/to/GenomeAlignmentTools/kent/bin;export PATH
```

# Installation
The perl and python scripts just need to be copied to a directory that is contained in your $PATH environment variable. 
Only chainCleaner and the modified chainNet (with the -rescore option) need to be compiled. 
Since both tools need the Kent (UCSC) source code, you need to set the KENTSRC_DIR variable to the path to the kent source code (should end with kent/src/).
```
export KENTSRC_DIR=/path/to/GenomeAlignmentTools/kent/src/
```
Afterwards, compile chainCleaner, chainNet and scoreChain
```
cd /path/to/GenomeAlignmentTools/src
export MACHTYPE=x86_64
make 
```

The subdirectory precompiledBinary_x86_64/ already contains Linux x86_64 compiled binaries.

# RepeatFiller
RepeatFiller [5] is a tool to incorporate newly-detected repeat-overlapping alignments into pairwise alignment chains [4]. Its runtime adds little to the computationally more expensive step of generating chains in pairwise whole-genome alignments. RepeatFiller circumvents the problem that considering all repeat-overlapping alignment seeds during whole genome alignment is computationally not feasible. Therefore, RepeatFiller only aligns local genomic regions that are bounded by colinear aligning blocks, as provided in the chains, which makes it feasible to consider all seeds including those that overlap repetitive regions. RepeatFiller application to mammalian genome alignment chains can add between 22 and 84 Mb of previously-undetected alignments that mostly originate from transposable elements [5]. This helps to comprehensively align repetitive regions and improves the annotation of conserved non-coding elements. 

__Usage:__
```
RepeatFiller.py -c in.chain -T2 target.2bit -Q2 query.2bit [options]
	in.chain	alignment file in chain format
	target.2bit	full path to a .2bit file of reference species
	query.2bit	full path to a .2bit file of query species
	[options]	call RepeatFiller.py --help to see the full list of parameters controlling 
	   		size and score thresholds of chains, and lastz and chaining options

example: RepeatFiller.py -c hg38.speTri2.all.chain -T2 hg38.2bit -Q2 speTri2.2bit
```

# Highly-sensitive local alignments
patchChain.perl performs a highly sensitive local pairwise alignment for loci flanked by aligning blocks [1,3]. 
Given an alignment chain [4], it considers all chains that pass the score and span filters (optional parameters), extracts all the unaligning loci and creates local alignment jobs. After executing these alignment jobs, the newly found and the original local alignments are combined and used to produce a new set of improved chains. 

This procedure is recommended for comparisons between species that are separated by >0.75 substitutions per neutral site [1].

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
The resulting output file can be used to get alignment nets by running chainNet [4].


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
chainCleaner example/hg38.mm10.chr1.chain.gz -tSizes=example/hg38.chrom.sizes \
  -qSizes=example/mm10.chrom.sizes example/hg38.2bit example/mm10.2bit \ 
  example/hg38.mm10.chr1.cleaned.chain example/removedSuspects.bed -linearGap=loose 


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
Given a set of alignment chains, chainNet produces alignment nets, which is a hierarchical collection of chains or parts of chains that attempt to capture only orthologous alignments [4]. 
The original chainNet implementation approximates the score of "sub-nets" (nets that come from a part of a chain and fill a gap in a higher-level net) by the fraction of aligning bases. This can lead to a bias in case the aligning blocks of a chain are not equally distributed. We implemented a new parameter "-rescore" in chainNet that computes the real score of each subnet [2].

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
While the netFilter program [4] removes nested nets in case their parent net does not satisfy the specified score and size criteria, NetFilterNonNested.perl applies a non-nested filtering procedure that considers and filters each net individually [1,3]. 
This avoids removing nested nets that would satisfy the specified criteria, even if a parent net is removed. 

__Usage:__
In [1], we applied the following filter criteria:

* UCSC "syntenic net" criteria (thresholds: minTopScore=300000, minSynScore=200000, minSynSize=20000, minSynAli=10000, maxFar=200000) in addition to keeping all nested nets that align 
to the same locus (inversions or local translocations) if they score higher than 5000 were applied to all placental mammal alignments that have well-assembled genomes 
`
NetFilterNonNested.perl -doUCSCSynFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 ref.query.net.gz > ref.query.filtered.net
`
* keeping nets that score higher than 100000 and keeping all nested nets that align to the same locus if they score higher than 5000 for placental mammals with less-well assembled genomes 
`
NetFilterNonNested.perl -doScoreFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 -minScore1 100000  ref.query.net.gz > ref.query.filtered.net
`
* keeping nets that score higher than 10000 and keeping all nested nets that align to the same locus if they score higher than 3000 for non-placental mammals 
`
NetFilterNonNested.perl -doScoreFilter -keepSynNetsWithScore 3000 -keepInvNetsWithScore 3000 -minScore1 10000 ref.query.net.gz > ref.query.filtered.net
`

Call NetFilterNonNested.perl without any parameters to see all filtering options.



# Kent source code (light-weight)

The 'kent' subdirectory contains a light-weight version of the Kent source code, which will only compile tools for generating or working
with genome alignment chains and nets.
To compile the light-weight version of Kent tools run make in the Kent source directory:

```
cd kent/src
make
```

For further details refer to 'kent/README.txt'.



# References
[1] Sharma V, Hiller M. [Increased alignment sensitivity improves the usage of genome alignments for comparative gene annotation](https://academic.oup.com/nar/article/45/14/8369/3875570/Increased-alignment-sensitivity-improves-the-usage?guestAccessKey=52de9f40-1303-4493-96e9-30def93a259c). Nucleic Acids Res., 45(14), 8369–8377, 2017

[2] Suarez H, Langer BE, Ladde P, Hiller M. [chainCleaner improves genome alignment specificity and sensitivity](https://academic.oup.com//bioinformatics/article/33/11/1596/2929344/chainCleaner-improves-genome-alignment-specificity?guestAccessKey=5b9b078a-39e3-49c2-807b-852efe66f366). Bioinformatics, 33(11):1596-1603, 2017

[3] Hiller M, Agarwal S, Notwell JH, Parikh R, Guturu H, Wenger AM, Bejerano G. [Computational methods to detect conserved non-genic elements in phylogenetically isolated genomes: application to zebrafish](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkt557). Nucleic Acids Res, 41(15):e151. 

[4] Kent WJ, Baertsch R, Hinrichs A, Miller W, Haussler D. [Evolution's cauldron: duplication, deletion, and rearrangement in the mouse and human genomes](http://www.pnas.org/content/100/20/11484.long). PNAS, 100(20):11484-9, 2003 

[5] Osipova E, Hecker N, Hiller M. RepeatFiller newly identifies megabases of aligning repetitive sequences and improves annotations of conserved non-exonic elements, submitted

