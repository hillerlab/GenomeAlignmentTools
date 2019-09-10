/* trix - text retrieval index.  Stuff for fast two level index
 * of text for fast word searches.  Generally you use the ixIxx program
 * to make the indexes. */

#ifndef TRIX_H
#define TRIX_H

struct trix
/* A two level index */
    {
    struct lineFile *lf;	/* Open file on first level index. */
    struct trixIxx *ixx;	/* Second level index in memory. */
    int ixxSize;		/* Size of second level index. */
    int ixxAlloc;	        /* Space allocated for index. */
    struct hash *wordHitHash;	/* Hash of word hitsLists, so search on "the the the" works fast. */
    boolean useUdc;            /* are we using UDC or lineFile */
    };

struct trixSearchResult
/* Result of a trix search. */
    {
    struct trixSearchResult *next;
    char *itemId;               /* ID of matching item */
    int unorderedSpan;          /* Minimum span in single doc with words in any order. */
    int orderedSpan;            /* Minimum span in single doc with words in search order. */
    int wordPos;		/* Position of word in doc more or less. */
    int leftoverLetters;	/* Number of leftover letters in words. */
    };

enum trixSearchMode
/* How stringent is the search? */
    {
    tsmExact,                   /* Require whole-word matches. */
    tsmExpand,                  /* Match words that differ from the search term only in the
                                 * last two letters stopping at a word boundary, or that are
                                 * the search word plus "ing". */
    tsmFirstFive                /* Like tsmExpand, but also match words that have the same
                                 * first 5 letters. */
    };

#define trixPrefixSize 5	/* Size of prefix in second level index. */

struct trix *trixOpen(char *ixFile);
/* Open up index.  Load second level index in memory. */

void trixClose(struct trix **pTrix);
/* Close up index and free up associated resources. */

struct trixSearchResult *trixSearch(struct trix *trix, int wordCount, char **words,
                                    enum trixSearchMode mode);
/* Return a list of items that match all words.  This will be sorted so that
 * multiple-word matches where the words are closer to each other and in the
 * right order will be first.  Single word matches will be prioritized so that those
 * closer to the start of the search text will appear before those later.
 * Do a trixSearchResultFreeList when done.  If mode is tsmExpand or tsmFirstFive then
 * this will match not only the input words, but also additional words that start with
 * the input words. */

void trixSearchResultFree(struct trixSearchResult **pTsr);
/* Free up data associated with trixSearchResult. */

void trixSearchResultFreeList(struct trixSearchResult **pList);
/* Free up a list of trixSearchResults. */

int trixSearchResultCmp(const void *va, const void *vb);
/* Compare two trixSearchResult in such a way that most relevant searches tend to be first. */

#endif //ndef TRIX_H
