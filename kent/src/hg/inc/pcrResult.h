/* pcrResult -- support for internal track of hgPcr results. */

/* Copyright (C) 2008 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef PCRRESULT_H
#define PCRRESULT_H

#include "targetDb.h"
#include "cart.h"

/* Definitions for track name and labels: */
#define PCR_RESULT_TRACK_NAME "hgPcrResult"
#define PCR_RESULT_TRACK_LABEL "PCR Results"
#define PCR_RESULT_TRACK_LONGLABEL "Your Sequence from PCR Search"

/* UI settings: */
#define PCR_RESULT_TARGET_STYLE PCR_RESULT_TRACK_NAME "_targetStyle"
#define PCR_RESULT_TARGET_STYLE_TRIM "trim"
#define PCR_RESULT_TARGET_STYLE_TALL "tall"
#define PCR_RESULT_TARGET_STYLE_DEFAULT PCR_RESULT_TARGET_STYLE_TRIM

char *pcrResultCartVar(char *db);
/* Returns the cart variable name for PCR result track info for db. 
 * Don't free the statically allocated result. */

boolean pcrResultParseCart(char *db, struct cart *cart, char **retPslFile,
			   char **retPrimerFile,
			   struct targetDb **retTarget);
/* Parse out hgPcrResult cart variable into components and make sure
 * they are valid.  If so, set *ret's and return TRUE.  Otherwise, null out 
 * *ret's and return FALSE.  ret's are ignored if NULL. */

void pcrResultGetPrimers(char *fileName, char **retFPrimer, char **retRPrimer);
/* Given a file whose first line is 2 words (forward primer, reverse primer)
 * set the ret's to upper-cased primer sequences.  
 * Do not free the statically allocated ret's. */

void pcrResultGetPsl(char *fileName, struct targetDb *target, char *item,
		     char *chrom, int itemStart, int itemEnd,
		     struct psl **retItemPsl, struct psl **retOtherPsls);
/* Read in psl from file.  If a psl matches the given item position, set 
 * retItemPsl to that; otherwise add it to retOtherPsls.  Die if no psl
 * matches the given item position. */

struct trackDb *pcrResultFakeTdb();
/* Construct a trackDb record for PCR Results track. */

char *pcrResultItemAccName(char *acc, char *name);
/* If a display name is given in addition to the acc, concatenate them
 * into a single name that must match a non-genomic target item's name
 * in the targetDb .2bit.  Do not free the result. */

char *pcrResultItemAccession(char *nameIn);
/* If nameIn contains a concatenate accession and display name, returns
 * just the accession.  Do not free the result.*/

char *pcrResultItemName(char *nameIn);
/* If nameIn contains a concatenated accession and display name, returns
 * just the name.  If accession only, returns NULL.  Do not free the result.*/

#endif /* PCRRESULT_H */
