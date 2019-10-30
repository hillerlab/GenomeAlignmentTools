# pseudoGeneLink.sql was originally generated by the autoSql program, which also 
# generated pseudoGeneLink.c and pseudoGeneLink.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#links a gene/pseudogene prediction to an ortholog or paralog.
CREATE TABLE pseudoGeneLink (
    chrom varchar(255) not null,	# Reference sequence chromosome or scaffold
    chromStart int unsigned not null,	# pseudogene alignment start position
    chromEnd int unsigned not null,	# pseudogene alignment end position
    name varchar(255) not null,	# Name of pseudogene
    score int unsigned not null,	# score of pseudogene with gene
    strand char(2) not null,	# + or -
    thickStart int unsigned not null,	# Start of where display should be thick (start codon)
    thickEnd int unsigned not null,	# End of where display should be thick (stop codon)
    reserved int unsigned not null,	# Always zero for now
    blockCount int not null,	# Number of blocks
    blockSizes longblob not null,	# Comma separated list of block sizes
    chromStarts longblob not null,	# Start positions relative to chromStart
    trfRatio float not null,	# ratio of tandem repeats
    type varchar(255) not null,	# type of evidence
    axtScore int not null,	# blastz score, gene mrna aligned to pseudogene
    gChrom varchar(255) not null,	# Chromosome name
    gStart int not null,	# gene alignment start position
    gEnd int not null,	# gene alignment end position
    gStrand char(2) not null,	# strand of gene
    exonCount int unsigned not null,	# # of exons in gene 
    geneOverlap int unsigned not null,	# bases overlapping
    polyA int unsigned not null,	# count of As in polyA
    polyAstart int not null,	# start of polyA, relative to end of pseudogene
    exonCover int not null,	# number of exons in Gene covered
    intronCount int unsigned not null,	# number of introns in pseudogene
    bestAliCount int unsigned not null,	# number of good mrnas aligning
    matches int unsigned not null,	# matches + repMatches
    qSize int unsigned not null,	# aligning bases in pseudogene
    qEnd int unsigned not null,	# end of cdna alignment
    tReps int unsigned not null,	# repeats in gene
    overlapRhesus int not null,	# percent of retro that breaks net relative to Rhesus
    overlapMouse int not null,	# percent of retro that breaks net relative to Mouse
    coverage int unsigned not null,	# % of bases that align to gene
    label int not null,	# 1=pseudogene,-1 not pseudogene -2 expressed retroGene
    milliBad int unsigned not null,	# milliBad score, pseudogene aligned to genome
    oldScore int unsigned not null,	# another heuristic
    oldIntronCount int not null,	# old simple intron count
    processedIntrons int not null,	# intron prcessed out of retrogene
    conservedSpliceSites int not null,	# conserved splice site count
    maxOverlap int not null,	# largest overlap with another mrna
    refSeq varchar(255) not null,	# Name of closest regSeq to gene
    rStart int not null,	# refSeq alignment start position
    rEnd int not null,	# refSeq alignment end position
    mgc varchar(255) not null,	# Name of closest mgc to gene
    mStart int not null,	# mgc alignment start position
    mEnd int not null,	# mgc alignment end position
    kgName varchar(255) not null,	# Name of closest knownGene to gene
    kStart int not null,	# kg alignment start position
    kEnd int not null,	# kg alignment end position
    overName varchar(255) not null,	# name of overlapping mrna
    overStart int not null,	# overlapping mrna start position
    overEnd int not null,	# overlapping mrna end position
    overStrand char(2) not null,	# strand of overlapping mrna
    overlapDog int not null,	# percent of retro that breaks net relative to dog
    posConf float not null,	# pvalue for positive
    polyAlen int unsigned not null,	# length of polyA
              #Indices
    PRIMARY KEY(chrom)
);