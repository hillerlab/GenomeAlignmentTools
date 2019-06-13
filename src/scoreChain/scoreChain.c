/* scoreChain - (re)score existing chains */
/* Michael Hiller, 2017 */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "dnaseq.h"
#include "nib.h"
#include "twoBit.h"
#include "fa.h"
#include "axt.h"
#include "portable.h"
#include "gapCalc.h"
#include "chainConnect.h"

static char const rcsid[] = "$Id: scoreChain.c,v 2.00 2017/08/16 Michael Hiller$";

/* Variables set via command line. */
char *gapFileName = NULL;
struct gapCalc *gapCalc = NULL;	/* Gap scoring scheme to use. */
struct axtScoreScheme *scoreScheme = NULL;
boolean doLocalScore = FALSE;
boolean forceLocalScore = FALSE;
boolean returnOnlyScore = FALSE;
boolean returnOnlyScoreAndCoords = FALSE;


/* for the target and query DNA seqs */
char *t2bit = NULL;    /* t and q nib or 2bit file */
char *q2bit = NULL;
/* hash keyed by target chromname, holding pointers to the seq */
struct hash *tSeqHash = NULL;
/* hash keyed by query chromname, holding pointers to the seq */
struct hash *qSeqHash = NULL;
/* hash keyed by query chromname, will hold a pointer to the reverse complemented seq. But we fill it on demand */
struct hash *qSeqMinusStrandHash = NULL;
/* file pointers to t and q 2bit file */
struct twoBitFile *ttbf;
struct twoBitFile *qtbf;


static struct optionSpec options[] = {
   {"scoreScheme", OPTION_STRING},
   {"linearGap", OPTION_STRING},
   {"doLocalScore", OPTION_BOOLEAN},     /* compute and return the local score only if the overall score is negative . */
   {"forceLocalScore", OPTION_BOOLEAN},  /* only compute and return local score. */
   {"returnOnlyScore", OPTION_BOOLEAN},  /* do not write the chain. Just return score. */
   {"returnOnlyScoreAndCoords", OPTION_BOOLEAN},  /* do not write the chain. Just return score and the chain span (tStart - tEnd). */
   {NULL, 0},
};

void usage()
/* Explain usage and exit. */
{
errAbort(
  "scoreChain - (re)score existing chains\n"
  "usage:\n"
  "   scoreChain in.chainFile reference.2bit query.2bit out.chain  -linearGap=loose|medium|filename\n"
  "Where reference.2bit and query.2bit are the names of a .2bit files for the reference and query\n"
  "options:\n"
  " Local score = we set score = 0 if score < 0 and return the max of the score that we reach for a chain\n" 
  ""
  "   -returnOnlyScore             default=FALSE. Just return chain ID{tab}globalScore{tab}localScore{tab}totalAligningBases, not the entire chain\n"
  "   -returnOnlyScoreAndCoords    default=FALSE. Just return chain ID{tab}chainStartInRef{tab}chainEndInRef{tab}localScore{tab}totalAligningBases, not the entire chain\n"
  "   -doLocalScore                default=FALSE. Only if the global score of a chain is negative, compute and output the local score in the chain file.\n"
  "   -forceLocalScore             default=FALSE. Always output the local score in the chain file.\n"
  ""
  "   -scoreScheme=fileName        Read the scoring matrix from a blastz-format file\n"
  "   -linearGap=<medium|loose|filename>    Specify type of linearGap to use.\n"
  "              *Must* specify this argument to one of these choices.\n"
  "              loose is chicken/human linear gap costs.\n"
  "              medium is mouse/human linear gap costs.\n"
  "              Or specify a piecewise linearGap tab delimited file.\n"
  "   sample linearGap file (loose)\n"
  "%s"
  , gapCalcSampleFileContents()
  );
  printf("%s\n", rcsid);
}



/****************************************************************************/
/* Return gaps from file, or default if fileName is NULL. */
/****************************************************************************/
struct gapCalc *gapCalcReadOrDefault(char *fileName)
{
if (fileName == NULL)
    errAbort("Must specify linear gap costs.  Use 'loose' or 'medium' for defaults\n");
return gapCalcFromFile(fileName);
}


/* ################################################################################### */
/* functions for loading and getting DNA seqs                                          */
/* ################################################################################### */

/****************************************************************
   Load sequence and add to hash, unless it is already loaded.  Do not reverse complement.
****************************************************************/
void loadSeq(char *seqPath, boolean isTarget, char *newName, struct hash *hash) {
   struct dnaSeq *seq;

   if (hashFindVal(hash, newName) == NULL)  {
      /* need to load */
      if (isTarget) {
         seq = twoBitReadSeqFrag(ttbf, newName, 0, 0);
         verbose(3, "\t\tLoaded %d bases of %s from %s\n", seq->size, newName, seqPath);
      } else {
         seq = twoBitReadSeqFrag(qtbf, newName, 0, 0);
         verbose(3, "\t\tLoaded %d bases of %s from %s\n", seq->size, newName, seqPath);
      }
      hashAdd(hash, newName, seq);
   }
}


/****************************************************************
   return dnaSeq struct from hash 
   abort if the dnaSeq is not in the hash
   if strand is - (must be a query seq), look up in qSeqMinusStrandHash and fill if not already present 
****************************************************************/
struct dnaSeq *getSeqFromHash(char *chrom, char strand, struct hash *hash) {
   struct dnaSeq *seq = hashFindVal(hash, chrom);

   if (seq == NULL)
      errAbort("ERROR: seq %s is not loaded in hash\n", chrom);

   if (strand == '+') {
      verbose(6, "\t\t\t\tgetSeqFromHash: load %s %c from hash\n", chrom, strand); fflush(stdout);
      return seq;
   } else {
      /* must be query seq */
      /* check if we have this seq in the qSeqMinusStrandHash */
      struct dnaSeq *rcSeq = hashFindVal(qSeqMinusStrandHash, chrom);
      if (rcSeq != NULL) {
         verbose(6, "\t\t\t\tgetSeqFromHash: load %s %c from revcomp hash\n", chrom, strand); fflush(stdout);
         return rcSeq;
      } else {
         /* get rc seq, store in hash and return */
         rcSeq = cloneDnaSeq(seq);
         reverseComplement(rcSeq->dna, rcSeq->size);
         hashAdd(qSeqMinusStrandHash, chrom, rcSeq);
         verbose(6, "\t\t\t\tgetSeqFromHash: add %s %c to revcomp hash\n", chrom, strand); fflush(stdout);
         return rcSeq;
      }
   }
   return NULL;
}


/****************************************************************
   free all dnaSeq struct from hash 
****************************************************************/
void freeDnaSeqHash(struct hash *hash) {
   struct hashEl *hel, *helList = hashElListHash(hash);
   
   for (hel = helList; hel != NULL; hel = hel->next)  {
      freeDnaSeq(hel->val);
   }
}



/* ################################################################################### */
/* functions for calculating chain scores                                              */
/* ################################################################################### */

/****************************************************************
  Calculate chain score locally. 
  That means, we reset score = 0 if score < 0 and return the max of the score that we reach for this chain.
  This is nearly identical to chainCalcScore() from chainConnect.c just a few modifications for the local score. 

  Since we traverse all blocks here, we also compute the total number of aligning bases
****************************************************************/
double chainCalcScoreLocal(struct chain *chain, struct axtScoreScheme *ss, struct gapCalc *gapCalc, struct dnaSeq *query, struct dnaSeq *target, int *retAliBases) {
   struct cBlock *b1, *b2;
   double score = 0;
   double maxScore = 0;
	int aliBases = 0;		/* this version also counts the #ali bases (total size of all blocks) */
   for (b1 = chain->blockList; b1 != NULL; b1 = b2) {
		aliBases += b1->tEnd - b1->tStart;

      score += chainScoreBlock(query->dna + b1->qStart, target->dna + b1->tStart, b1->tEnd - b1->tStart, ss->matrix);

      if (score > maxScore) 
         maxScore = score; 

      b2 = b1->next;
      if (b2 != NULL) {
         score -=  gapCalcCost(gapCalc, b2->qStart - b1->qEnd, b2->tStart - b1->tEnd);
         if (score < 0) 
            score = 0;  
      }
   }
   *retAliBases = aliBases;
	return maxScore;
}



/****************************************************************
   calculate the score of the given chain, both the local (always >0) and the global score
   sets chain->score to the global score (needed for final rescoring of the breaking chains)
   returns both global and local score in the double pointers and the total number of aligning bases (aliBases)
****************************************************************/
double getChainScore (struct chain *chain, double *globalScore, double *localScore, int *aliBases) {
   struct dnaSeq *qSeq = NULL, *tSeq = NULL;
//   printf("\tcalc score for chain with ID %d\n", chain->id);fflush(stdout);

   /* load the seqs */
   qSeq = getSeqFromHash(chain->qName, chain->qStrand, qSeqHash);
   tSeq = getSeqFromHash(chain->tName, '+', tSeqHash);

   chain->score = chainCalcScore(chain, scoreScheme, gapCalc, qSeq, tSeq);
   *globalScore = chain->score;
   *localScore = chainCalcScoreLocal(chain, scoreScheme, gapCalc, qSeq, tSeq, aliBases);

   return chain->score;
}







/****************************************************************************/
/* main */
/****************************************************************************/
int main(int argc, char *argv[])
/* Process command line. */
{

optionInit(&argc, argv, options);

char *scoreSchemeName = NULL;
optionHash(&argc, argv);
gapFileName = optionVal("linearGap", NULL);
scoreSchemeName = optionVal("scoreScheme", NULL);
doLocalScore = optionExists("doLocalScore");
forceLocalScore = optionExists("forceLocalScore");
returnOnlyScore = optionExists("returnOnlyScore");
returnOnlyScoreAndCoords = optionExists("returnOnlyScoreAndCoords");

if (argc != 5)
	usage();

if (returnOnlyScore && returnOnlyScoreAndCoords)
	errAbort("ERROR: You cannot specify both returnOnlyScore and returnOnlyScoreAndCoords\n");


/* load score scheme */
if (scoreSchemeName != NULL) {
   verbose(2, "Reading scoring matrix from %s\n", scoreSchemeName);
   scoreScheme = axtScoreSchemeRead(scoreSchemeName);
} else {
   scoreScheme = axtScoreSchemeDefault();
}
/* load gap costs */
if (gapFileName == NULL)
    errAbort("Must specify linear gap costs.  Use 'loose' or 'medium' for defaults\n");
gapCalc = gapCalcFromFile(gapFileName);

dnaUtilOpen();


/* open the 2bit files */
t2bit = argv[2];
q2bit = argv[3];
/* test if the seq files exist.  */
if (! fileExists(t2bit))
   errAbort("ERROR: target 2bit file or nib directory %s does not exist\n", t2bit);
if (! fileExists(q2bit))
   errAbort("ERROR: query 2bit file or nib directory %s does not exist\n", q2bit);


/* open the 2bit files if t or q seq is given as a 2bit file */
if (twoBitIsFile(t2bit)) {
   ttbf = twoBitOpen(t2bit);
}else{
	errAbort("ERROR: only 2bit files are supported, not %s\n", t2bit);
}
if (twoBitIsFile(q2bit)) {
   qtbf = twoBitOpen(q2bit);
}else{
	errAbort("ERROR: only 2bit files are supported, not %s\n", q2bit);
}
tSeqHash = newHash(0);
qSeqHash = newHash(0);
qSeqMinusStrandHash = newHash(0);


/*  open output file */
FILE *f = mustOpen(argv[4], "w");

/********************/
/* score all chains */
struct lineFile *lf = lineFileOpen(argv[1], TRUE);
struct chain *chain;
while ((chain = chainRead(lf)) != NULL) {

	/* load the t and q seq from the 2bit file, if not already in our hash */
   loadSeq(t2bit, TRUE, chain->tName, tSeqHash);
   loadSeq(q2bit, FALSE, chain->qName, qSeqHash);

	verbose(2, "adjust score for chain %d (t %s %d-%d  q %c %s %d-%d) from %f to ", chain->id, chain->tName, chain->tStart, chain->tEnd, chain->qStrand, chain->qName, chain->qStart, chain->qEnd, chain->score);
   
	double globalScore, localScore;
	int aliBases; 
	getChainScore(chain, &globalScore, &localScore, &aliBases);

	if (forceLocalScore) {
		chain->score = localScore;
	} else {
		chain->score = globalScore;
		if (chain->score <= 0 && doLocalScore) {
			verbose(2, "\tSCORE IS NEGATIVE --> doLocal is set --> set global score %f to local score %f\n ", chain->score, localScore);
			chain->score = localScore;
		}
	}
	verbose(2, "Final score for chainID %d: %f\n ", chain->id, chain->score);

	if (returnOnlyScore) {
		fprintf(f, "%d\t%1.0f\t%1.0f\t%d\n", chain->id, globalScore, localScore, aliBases);
	}else if (returnOnlyScoreAndCoords) {
		fprintf(f, "%d\t%d\t%d\t%1.0f\t%1.0f\t%d\n", chain->id, chain->tStart, chain->tEnd, globalScore, localScore, aliBases);
	}else{
		chainWrite(chain, f);
	}
}
carefulClose(&f);
lineFileClose(&lf);

/* free some memory */
freeDnaSeqHash(tSeqHash);
freeDnaSeqHash(qSeqHash);
freeDnaSeqHash(qSeqMinusStrandHash); 
freeHash(&tSeqHash);
freeHash(&qSeqHash);
freeHash(&qSeqMinusStrandHash);


return 0;
}
