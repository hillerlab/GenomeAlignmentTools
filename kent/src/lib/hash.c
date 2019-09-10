/* Hash.c - implements hashing.  See hash.h for usage comments.
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include "localmem.h"
#include "hash.h"
#include "obscure.h"
#include "dystring.h"


/*
 * Hash a string key.  This code is taken from Tcl interpreter. I was borrowed
 * after discovering a lot of collisions and poor utilization of the table
 * when hashing accessions.
 *
 * This function was compared to Bob Jenkins' lookup2 hash function and
 * (http://burtleburtle.net/bob/hash/) and Paul Hsieh's SuperFast
 * hash function (http://www.azillionmonkeys.com/qed/hash.html).
 * Both of those functions provided better utilization of the table,
 * but were also more expensive, so the Tcl function was used.
 * If hashing of binary keys is implemented, SuperFast hash should
 * be considered.
 *
 * for an explanation of this function, see HashStringKey() in the
 * Tcl source file, generic/tclHash.c, available from
 * http://tcl.sourceforge.net/.
 *
 * The Tcl code is:
 * Copyright (c) 1991-1993 The Regents of the University of California.
 * Copyright (c) 1994 Sun Microsystems, Inc.
 *
 * See the file "license.terms" (in the Tcl distribution) for complete
 * license (which is a BSD-style license).
 *
 * Since hashCrc() is in use elsewhere, 
 * a new function hashString() was created for use in hash table.
 * -- markd
 */
bits32 hashString(char *string)
/* Compute a hash value of a string. */
{
char *keyStr = string;
unsigned int result = 0;
int c;

while ((c = *keyStr++) != '\0')
    {
    result += (result<<3) + c;
    }
return result;
}

bits32 hashCrc(char *string)
/* Returns a CRC value on string. */
{
unsigned char *us = (unsigned char *)string;
unsigned char c;
bits32 shiftAcc = 0;
bits32 addAcc = 0;

while ((c = *us++) != 0)
    {
    shiftAcc <<= 2;
    shiftAcc += c;
    addAcc += c;
    }
return shiftAcc + addAcc;
}

struct hashEl *hashLookup(struct hash *hash, char *name)
/* Looks for name in hash table. Returns associated element,
 * if found, or NULL if not.  If there are multiple entries
 * for name, the last one added is returned (LIFO behavior).
 */
{
struct hashEl *el = hash->table[hashString(name)&hash->mask];
while (el != NULL)
    {
    if (strcmp(el->name, name) == 0)
        break;
    el = el->next;
    }
return el;
}

struct hashEl *hashLookupUpperCase(struct hash *hash, char *name)
/* Lookup upper cased name in hash. (Assumes all elements of hash
 * are themselves already in upper case.) */
{
char s[256];
safef(s, sizeof(s), "%s", name);
touppers(s);
return hashLookup(hash, s);
}


struct hashEl *hashLookupNext(struct hashEl *hashEl)
/* Find the next occurance of name that may occur in the table multiple times,
 * or NULL if not found.  Use hashLookup to find the first occurrence.  Elements
 * are returned in LIFO order.
 */
{
struct hashEl *el = hashEl->next;
while (el != NULL)
    {
    if (strcmp(el->name, hashEl->name) == 0)
        break;
    el = el->next;
    }
return el;
}

struct hashEl *hashAddN(struct hash *hash, char *name, int nameSize, void *val)
/* Add name of given size to hash (no need to be zero terminated) */
{
struct hashEl *el;
if (hash->lm) 
    el = lmAlloc(hash->lm, sizeof(*el));
else
    AllocVar(el);
el->hashVal = hashString(name);
int hashVal = el->hashVal & hash->mask;
if (hash->lm)
    {
    el->name = lmAlloc(hash->lm, nameSize+1);
    memcpy(el->name, name, nameSize);
    }
else
    el->name = cloneStringZ(name, nameSize);
el->val = val;
el->next = hash->table[hashVal];
hash->table[hashVal] = el;
hash->elCount += 1;
if (hash->autoExpand && hash->elCount > (int)(hash->size * hash->expansionFactor))
    {
    /* double the size */
    hashResize(hash, digitsBaseTwo(hash->size));
    }
return el;
}

struct hashEl *hashAdd(struct hash *hash, char *name, void *val)
/* Add new element to hash table.  If an item with name, already exists, a new
 * item is added in a LIFO manner.  The last item added for a given name is
 * the one returned by the hashLookup functions.  hashLookupNext must be used
 * to find the preceding entries for a name.
 */
{
return hashAddN(hash, name, strlen(name), val);
}

boolean hashMayRemove(struct hash *hash, char *name)
/* Remove item of the given name from hash table, if present.
 * Return true if it was present */
{
return (hashRemove(hash, name) != NULL);
}

void hashMustRemove(struct hash *hash, char *name)
/* Remove item of the given name from hash table, or error
 * if not present */
{
if (hashRemove(hash, name) == NULL)
    errAbort("attempt to remove non-existant %s from hash", name);
}

void freeHashEl(struct hashEl *hel)
/* Free hash element. Use only on non-local memory version. */
{
freeMem(hel->name);
freeMem(hel);
}

void *hashRemove(struct hash *hash, char *name)
/* Remove item of the given name from hash table. 
 * Returns value of removed item, or NULL if not in the table.
 * If their are multiple entries for name, the last one added
 * is removed (LIFO behavior).
 */
{
struct hashEl *hel;
void *ret;
struct hashEl **pBucket = &hash->table[hashString(name)&hash->mask];
for (hel = *pBucket; hel != NULL; hel = hel->next)
    if (sameString(hel->name, name))
        break;
if (hel == NULL)
    return NULL;
ret = hel->val;
if (slRemoveEl(pBucket, hel))
    {
    hash->elCount -= 1;
    if (!hash->lm)
	freeHashEl(hel);
    }
return ret;
}

struct hashEl *hashAddUnique(struct hash *hash, char *name, void *val)
/* Add new element to hash table. Squawk and die if not unique */
{
if (hashLookup(hash, name) != NULL)
    errAbort("%s duplicated, aborting", name);
return hashAdd(hash, name, val);
}

struct hashEl *hashAddSaveName(struct hash *hash, char *name, void *val, char **saveName)
/* Add new element to hash table.  Save the name of the element, which is now
 * allocated in the hash table, to *saveName.  A typical usage would be:
 *    AllocVar(el);
 *    hashAddSaveName(hash, name, el, &el->name);
 */
{
struct hashEl *hel = hashAdd(hash, name, val);
*saveName = hel->name;
return hel;
}

struct hashEl *hashStore(struct hash *hash, char *name)
/* If element in hash already return it, otherwise add it
 * and return it. */
{
struct hashEl *hel;
if ((hel = hashLookup(hash, name)) != NULL)
    return hel;
return hashAdd(hash, name, NULL);
}

char  *hashStoreName(struct hash *hash, char *name)
/* If element in hash already return it, otherwise add it
 * and return it. */
{
struct hashEl *hel;
if (name == NULL)
    return NULL;
if ((hel = hashLookup(hash, name)) != NULL)
    return hel->name;
return hashAdd(hash, name, NULL)->name;
}

int hashIntVal(struct hash *hash, char *name)
/* Return integer value associated with name in a simple 
 * hash of ints. */
{
void *val = hashMustFindVal(hash, name);
return ptToInt(val);
}

int hashIntValDefault(struct hash *hash, char *name, int defaultInt)
/* Return integer value associated with name in a simple 
 * hash of ints or defaultInt if not found. */
{
struct hashEl *hel = hashLookup(hash, name);
if(hel == NULL)
    return defaultInt;
return ptToInt(hel->val);
}

void *hashMustFindVal(struct hash *hash, char *name)
/* Lookup name in hash and return val.  Abort if not found. */
{
struct hashEl *hel = hashLookup(hash, name);
if (hel == NULL)
    errAbort("hashMustFindVal: '%s' not found", name);
return hel->val;
}

void *hashFindVal(struct hash *hash, char *name)
/* Look up name in hash and return val or NULL if not found. */
{
struct hashEl *hel = hashLookup(hash, name);
if (hel == NULL)
    return NULL;
return hel->val;
}

void *hashOptionalVal(struct hash *hash, char *name, void *usual)
/* Look up name in hash and return val, or usual if not found. */
{
struct hashEl *hel = hashLookup(hash, name);
if (hel == NULL)
    return usual;
else
    return hel->val;
}

void *hashFindValUpperCase(struct hash *hash, char *name)
/* Lookup upper cased name in hash and return val or return NULL if not found.
 * (Assumes all elements of hash are themselves already in upper case.) */
{
struct hashEl *hel = hashLookupUpperCase(hash, name);
if (hel == NULL)
    return NULL;
return hel->val;
}

char *hashMustFindName(struct hash *hash, char *name)
/* Return name as stored in hash table (in hel->name). 
 * Abort if not found. */
{
struct hashEl *hel = hashLookup(hash, name);
if (hel == NULL)
    errAbort("hashMustFindName: '%s' not found", name);
return hel->name;
}

struct hashEl *hashAddInt(struct hash *hash, char *name, int val)
/* Store integer value in hash */
{
char *pt = NULL;
return hashAdd(hash, name, pt + val);
}


void hashIncInt(struct hash *hash, char *name)
/* Increment integer value in hash */
{
struct hashEl *hel = hashLookup(hash, name);
if (hel == NULL)
  {
  hashAddInt(hash, name, 1);
  }
else
  {
  hel->val = ((char *)hel->val)+1;
  /* The much simpler ++hel->val works for gnu C, but really adding one to a void pointer
   * I think is not well defined. */
  }
}

long long hashIntSum(struct hash *hash)
/* Return sum of all the ints in a hash of ints. */
{
long long sum = 0;
int i;
struct hashEl *hel;
for (i=0; i<hash->size; ++i)
    {
    for (hel = hash->table[i]; hel != NULL; hel = hel->next)
	{
	int num = ptToInt(hel->val);
	sum += (long long)num;
	}
    }
return sum;
}

void hashIntReset(struct hash *hash)
/* Reset all values in hash of ints to 0.  Reset element count to 0. */
{
memset(hash->table, 0, hash->size * sizeof hash->table[0]);
hash->elCount = 0;
}

static int adjustPowerOfTwoSize(int powerOfTwoSize)
/* If powerOfTwoSize is 0, return a reasonable default to use instead.  If powerOfTwoSize
 * is out of range, errAbort.  Otherwise return powerOfTwoSize. */
{
if (powerOfTwoSize == 0)
    powerOfTwoSize = 12;
if (powerOfTwoSize > hashMaxSize || powerOfTwoSize < 0)
    errAbort("hash powerOfTwoSize must be >= 0 and <= %d, but %d was passed in.",
             hashMaxSize, powerOfTwoSize);
return powerOfTwoSize;
}

struct hash *newHashLm(int powerOfTwoSize, struct lm *lm)
/* Returns new hash table using the given lm.  Recommended lm block size is 256B to 64kB,
 * depending on powerOfTwoSize. */
{
struct hash *hash = lm ? lmAlloc(lm, sizeof(*hash)) : needMem(sizeof(*hash));
powerOfTwoSize = adjustPowerOfTwoSize(powerOfTwoSize);
hash->powerOfTwoSize = powerOfTwoSize;
hash->size = (1<<powerOfTwoSize);
hash->lm = lm;
hash->mask = hash->size-1;
if (lm)
    lmAllocArray(hash->lm, hash->table, hash->size);
else
    AllocArray(hash->table, hash->size);
hash->autoExpand = TRUE;
hash->expansionFactor = defaultExpansionFactor;   /* Expand when elCount > size*expansionFactor */
return hash;
}

struct hash *newHashExt(int powerOfTwoSize, boolean useLocalMem)
/* Returns new hash table. Uses local memory optionally. */
{
struct hash *hash = NULL;
if (useLocalMem)
    {
    int memBlockPower = 16;
    powerOfTwoSize = adjustPowerOfTwoSize(powerOfTwoSize);
    /* Make size of memory block for allocator vary between
     * 256 bytes and 64k depending on size of table. */
    if (powerOfTwoSize < 8)
        memBlockPower = 8;
    else if (powerOfTwoSize < 16)
        memBlockPower = powerOfTwoSize;
    struct lm *ownLm = lmInit(1<<memBlockPower);
    hash = newHashLm(powerOfTwoSize, ownLm);
    hash->ownLm = TRUE;
    }
else
    hash = newHashLm(powerOfTwoSize, NULL);
return hash;
}

void hashReverseAllBucketLists(struct hash *hash)
/* Reverse all hash bucket list.  You might do this to
 * get them back in the same order things were added to the hash */
{
int i;
for (i=0; i<hash->size; ++i)
    {
    struct hashEl *hel = hash->table[i];
    if (hel != NULL && hel->next != NULL)	    
	slReverse(&hash->table[i]);
    }
}

void hashResize(struct hash *hash, int powerOfTwoSize)
/* Resize the hash to a new size */
{
int oldHashSize = hash->size;
struct hashEl **oldTable = hash->table;

if (powerOfTwoSize > hashMaxSize)
    powerOfTwoSize =  hashMaxSize;
powerOfTwoSize = adjustPowerOfTwoSize(powerOfTwoSize);
if (hash->powerOfTwoSize == powerOfTwoSize)
    return;

assert(powerOfTwoSize <= hashMaxSize && powerOfTwoSize > 0);
hash->powerOfTwoSize = powerOfTwoSize;
hash->size = (1<<powerOfTwoSize);
hash->mask = hash->size-1;

if (hash->lm)
    lmAllocArray(hash->lm, hash->table, hash->size);
else
    AllocArray(hash->table, hash->size);

int i;
struct hashEl *hel, *next;
for (i=0; i<oldHashSize; ++i)
    {
    for (hel = oldTable[i]; hel != NULL; hel = next)
	{
	next = hel->next;
	int hashVal = hel->hashVal & hash->mask;
	hel->next = hash->table[hashVal];
	hash->table[hashVal] = hel;
	}
    }
/* restore original list order */
hashReverseAllBucketLists(hash);

if (!hash->lm)
    freeMem(oldTable);
hash->numResizes++;
}


struct hash *hashFromSlNameList(void *list)
/* Create a hash out of a list of slNames. */
{
struct hash *hash = NULL;
struct slName *namedList = list, *item;
if (!list)
    return NULL;
hash = newHash(0);
for (item = namedList; item != NULL; item = item->next)
    hashAdd(hash, item->name, item);
return hash;
}

struct hash *hashSetFromSlNameList(void *list)
/* Create a hashSet (hash with only keys) out of a list of slNames. */
{
struct hash *hash = NULL;
struct slName *namedList = list, *item;
if (!list)
    return NULL;
hash = newHash(0);
for (item = namedList; item != NULL; item = item->next)
    hashAdd(hash, item->name, NULL);
return hash;
}

void hashTraverseEls(struct hash *hash, void (*func)(struct hashEl *hel))
/* Apply func to every element of hash with hashEl as parameter. */
{
int i;
struct hashEl *hel;
for (i=0; i<hash->size; ++i)
    {
    for (hel = hash->table[i]; hel != NULL; hel = hel->next)
	func(hel);
    }
}

void hashTraverseVals(struct hash *hash, void (*func)(void *val))
/* Apply func to every element of hash with hashEl->val as parameter. */
{
int i;
struct hashEl *hel;
for (i=0; i<hash->size; ++i)
    {
    for (hel = hash->table[i]; hel != NULL; hel = hel->next)
	func(hel->val);
    }
}

int hashElCmp(const void *va, const void *vb)
/* Compare two hashEl by name. */
{
const struct hashEl *a = *((struct hashEl **)va);
const struct hashEl *b = *((struct hashEl **)vb);
return strcmp(a->name, b->name);
}

int hashElCmpWithEmbeddedNumbers(const void *va, const void *vb)
/* Compare two hashEl by name sorting including numbers within name,
 * suitable for chromosomes, genes, etc. */
{
const struct hashEl *a = *((struct hashEl **)va);
const struct hashEl *b = *((struct hashEl **)vb);
return cmpStringsWithEmbeddedNumbers(a->name, b->name);
}

int hashElCmpIntValDesc(const void *va, const void *vb)
/* Compare two hashEl from a hashInt type hash, with highest integer values
 * comingFirst. */
{
struct hashEl *a = *((struct hashEl **)va);
struct hashEl *b = *((struct hashEl **)vb);
return b->val - a->val;
}

void *hashElFindVal(struct hashEl *list, char *name)
/* Look up name in hashEl list and return val or NULL if not found. */
{
struct hashEl *el;
for (el = list; el != NULL; el = el->next)
    {
    if (strcmp(el->name, name) == 0)
        return el->val;
    }
return NULL;
}

struct hashEl *hashElListHash(struct hash *hash)
/* Return a list of all elements of hash.   Free return with hashElFreeList. */
{
int i;
struct hashEl *hel, *dupe, *list = NULL;
for (i=0; i<hash->size; ++i)
    {
    for (hel = hash->table[i]; hel != NULL; hel = hel->next)
	{
	dupe = CloneVar(hel);
	slAddHead(&list, dupe);
	}
    }
return list;
}


void hashElFree(struct hashEl **pEl)
/* Free hash el list returned from hashListAll.  (Don't use
 * this internally.) */
{
freez(pEl);
}

void hashElFreeList(struct hashEl **pList)
/* Free hash el list returned from hashListAll.  (Don't use
 * this internally. */
{
slFreeList(pList);
}

struct hashCookie hashFirst(struct hash *hash)
/* Return an object to use by hashNext() to traverse the hash table.
 * The first call to hashNext will return the first entry in the table. */
{
struct hashCookie cookie;
cookie.hash = hash;
cookie.idx = 0;
cookie.nextEl = NULL;

/* find first entry */
for (cookie.idx = 0;
     (cookie.idx < hash->size) && (hash->table[cookie.idx] == NULL);
     cookie.idx++)
    continue;  /* empty body */
if (cookie.idx < hash->size)
    cookie.nextEl = hash->table[cookie.idx];
return cookie;
}

struct hashEl* hashNext(struct hashCookie *cookie)
/* Return the next entry in the hash table, or NULL if no more. Do not modify
 * hash table while this is being used. */
{
/* NOTE: if hashRemove were coded to track the previous entry during the
 * search and then use it to do the remove, it would be possible to
 * remove the entry returned by this method */
struct hashEl *retEl = cookie->nextEl;
if (retEl == NULL)
    return NULL;  /* no more */

/* find next entry */
cookie->nextEl = retEl->next;
if (cookie->nextEl == NULL)
    {
    for (cookie->idx++; (cookie->idx < cookie->hash->size)
             && (cookie->hash->table[cookie->idx] == NULL); cookie->idx++)
        continue;  /* empty body */
    if (cookie->idx < cookie->hash->size)
        cookie->nextEl = cookie->hash->table[cookie->idx];
    }
return retEl;
}

void* hashNextVal(struct hashCookie *cookie)
/* Return the next value in the hash table, or NULL if no more. Do not modify
 * hash table while this is being used. */
{
struct hashEl *hel = hashNext(cookie);
if (hel == NULL)
    return NULL;
else
    return hel->val;
}

char *hashNextName(struct hashCookie *cookie)
/* Return the next name in the hash table, or NULL if no more. Do not modify
 * hash table while this is being used. */
{
struct hashEl *hel = hashNext(cookie);
if (hel == NULL)
    return NULL;
else
    return hel->name;
}

void freeHash(struct hash **pHash)
/* Free up hash table. */
{
struct hash *hash = *pHash;
if (hash == NULL)
    return;
if (hash->lm)
    {
    if (hash->ownLm)
        lmCleanup(&hash->lm);
    *pHash = NULL;
    }
else
    {
    int i;
    struct hashEl *hel, *next;
    for (i=0; i<hash->size; ++i)
	{
	for (hel = hash->table[i]; hel != NULL; hel = next)
	    {
	    next = hel->next;
	    freeHashEl(hel);
	    }
	}
    freeMem(hash->table);
    freez(pHash);
    }
}


void freeHashAndVals(struct hash **pHash)
/* Free up hash table and all values associated with it.
 * (Just calls freeMem on each hel->val) */
{
struct hash *hash;
if ((hash = *pHash) != NULL)
    {
    hashTraverseVals(hash, freeMem);
    freeHash(pHash);
    }
}

void hashFreeWithVals(struct hash **pHash, void (freeFunc)())
/* Free up hash table and all values associated with it. freeFunc is a
 * function to free an entry, should take a pointer to a pointer to an
 * entry. */
{
struct hash *hash = *pHash;
if (hash != NULL)
    {
    struct hashCookie cookie = hashFirst(hash);
    struct hashEl *hel;
    while ((hel = hashNext(&cookie)) != NULL)
        freeFunc(&hel->val);
    hashFree(pHash);
    }
}

void hashFreeList(struct hash **pList)
/* Free up a list of hashes. */
{
struct hash *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    hashFree(&el);
    }
*pList = NULL;
}

static int bucketLen(struct hashEl *hel)
/* determine how many elements are in a hash bucket */
{
int nel = 0;
for (; hel != NULL; hel = hel->next)
    nel++;
return nel;
}

void hashHisto(struct hash *hash, char *fname)
/* Output bucket usage counts to a file for producing a histogram  */
{
FILE* fh = mustOpen(fname, "w");
int i;

for (i=0; i<hash->size; ++i)
    fprintf(fh, "%d\n", bucketLen(hash->table[i]));
carefulClose(&fh);
}

void hashPrintStats(struct hash *hash, char *label, FILE *fh)
/* print statistic about a hash table */
{
// count up usage
int i, occupiedCnt = 0, maxBucket = 0;
for (i=0; i<hash->size; ++i)
    {
    if (hash->table[i] != NULL)
        occupiedCnt++;
    int sz = bucketLen(hash->table[i]);
    maxBucket = max(maxBucket, sz);
    }

fprintf(fh, "hashTable\t%s\n", label);
fprintf(fh, "tableSize\t%d\t%d\n", hash->size, hash->powerOfTwoSize);
fprintf(fh, "numElements\t%d\n", hash->elCount);
fprintf(fh, "occupied\t%d\t%0.4f\n", occupiedCnt, ((hash->size == 0) ? 0.0 : ((float)occupiedCnt)/hash->size));
fprintf(fh, "maxBucket\t%d\n", maxBucket);
fprintf(fh, "numResizes\t%d\n", hash->numResizes);
fprintf(fh, "\n");
}

struct hashEl *hashReplace(struct hash *hash, char *name, void *val)
/* Replace an existing element in hash table, or add it if not present. */
{
if (hashLookup(hash, name))
    hashRemove(hash, name);
return hashAdd(hash, name, val);
}

char *hashToRaString(struct hash *hash)
/* Convert hash to string in ra format. */
{
struct hashEl *el, *list = hashElListHash(hash);
struct dyString *dy = dyStringNew(0);
slSort(&list, hashElCmp);
for (el = list; el != NULL; el = el->next)
   {
   dyStringAppend(dy, el->name);
   dyStringAppendC(dy, ' ');
   dyStringAppend(dy, el->val);
   dyStringAppendC(dy, '\n');
   }
hashElFreeList(&list);
return dyStringCannibalize(&dy);
}

int hashNumEntries(struct hash *hash)
/* count the number of entries in a hash */
{
int n = 0, i;
for (i=0; i<hash->size; ++i)
    n += bucketLen(hash->table[i]);
return n;
}

struct hash *hashFromString(char *string)
/* parse a whitespace-separated string with tuples in the format name=val or
 * name="val" to a hash name->val */
{
if (string==NULL)
    return NULL;

struct slPair *keyVals = slPairListFromString(string, TRUE);
if (keyVals==NULL)
    return NULL;

struct hash *nameToVal = newHash(0);
struct slPair *kv;
for (kv = keyVals; kv != NULL; kv = kv->next)
    hashAdd(nameToVal, kv->name, kv->val);
return nameToVal;
}

struct hash *hashFromNameArray(char **nameArray, int nameCount)
/* Create a NULL valued hash on all names in array */
{
struct hash *hash = hashNew(0);
int i;
for (i=0; i<nameCount; ++i)
    hashAdd(hash, nameArray[i], NULL);
return hash;
}

struct hash *hashFromNameValArray(char *nameVal[][2], int nameValCount)
/* Make up a hash from nameVal array */
{
struct hash *hash = newHash(0);
int i;
for (i=0; i<nameValCount; ++i)
    hashAdd(hash, nameVal[i][0], nameVal[i][1]);
return hash;
}


