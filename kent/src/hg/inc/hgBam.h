/* hgBam -- interface to binary alignment format files using Heng Li's samtools lib. */
/* these are the components that link to MySQL */

/* Copyright (C) 2011 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef HGBAM_H
#define HGBAM_H

#ifndef BAMFILE_H
#include "bamFile.h"
#endif

#ifndef SAMALIGNMENT_H
#include "samAlignment.h"
#endif

#ifndef JKSQL_H
#include "jksql.h"
#endif

#include "trackDb.h"

char *bamFileNameFromTable(struct sqlConnection *conn, char *table, char *bamSeqName);
/* Return file name from table.  If table has a seqName column, then grab the 
 * row associated with bamSeqName (which can be e.g. '1' not 'chr1' if that is the
 * case in the bam file). */

struct samAlignment *bamFetchSamAlignment(char *fileOrUrl, char *chrom, int start, int end,
	struct lm *lm);
/* Fetch region as a list of samAlignments - which is more or less an unpacked
 * bam record.  Results is allocated out of lm, since it tends to be large... */

struct samAlignment *bamFetchSamAlignmentPlus(char *fileOrUrl, char *chrom, int start, int end,
	struct lm *lm, char *refUrl, char *cacheDir );
/* Fetch region as a list of samAlignments - which is more or less an unpacked
 * bam record.  Results is allocated out of lm, since it tends to be large... */

struct samAlignment *bamAndIndexFetchSamAlignmentPlus(char *fileOrUrl, char *baiUrl, char *chrom, int start, int end,
	struct lm *lm,  char *refUrl, char *cacheDir);
/* Like bamFetchSamAlignmentPlus but can specify bai index file url in addition to the bam file */

struct samAlignment *bamReadNextSamAlignments(samfile_t *fh, bam_hdr_t *header,  int count, struct lm *lm);
/* Read next count alignments in SAM format, allocated in lm.  May return less than
 * count at end of file. */

struct ffAli *bamToFfAli(const bam1_t *bam, struct dnaSeq *target, int targetOffset,
			 boolean useStrand, char **retQSeq);
/* Convert from bam to ffAli format.  If retQSeq is non-null, set it to the 
 * query sequence into which ffAli needle pointers point. */

#endif//ndef HGBAM_H
