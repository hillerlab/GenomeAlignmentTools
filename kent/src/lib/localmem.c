/* LocalMem.c - local memory routines. 
 * 
 * These routines are meant for the sort of scenario where
 * a lot of little to medium size pieces of memory are
 * allocated, and then disposed of all at once.
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */


#include "common.h"
#include "localmem.h"


struct lm
    {
    struct lmBlock *blocks;
    size_t blockSize;
    size_t allignMask;
    size_t allignAdd;
    boolean doMemoryAllocs; // if true, do our own memory allocs, otherwise use passed in pointer
    };

struct lmBlock
    {
    struct lmBlock *next;
    char *free;
    char *end;
    char *extra;
    };

static struct lmBlock *newBlock(struct lm *lm, size_t reqSize)
/* Allocate a new block of at least reqSize */
{
size_t size = (reqSize > lm->blockSize ? reqSize : lm->blockSize);
size_t fullSize = size + sizeof(struct lmBlock);
struct lmBlock *mb = needLargeZeroedMem(fullSize);
if (mb == NULL)
    errAbort("Couldn't allocate %lld bytes", (long long)fullSize);
mb->free = (char *)(mb+1);
mb->end = ((char *)mb) + fullSize;
mb->next = lm->blocks;
lm->blocks = mb;
return mb;
}

struct lm *lmGuts(int blockSize, void *mem)
/* Create a local memory pool. */
{
struct lm *lm;
int aliSize = sizeof(long);
if (aliSize < sizeof(double))
    aliSize = sizeof(double);
if (aliSize < sizeof(void *))
    aliSize = sizeof(void *);
lm = needMem(sizeof(*lm));
lm->blocks = NULL;
if (blockSize <= 0)
    blockSize = (1<<14);    /* 16k default. */
lm->blockSize = blockSize;
lm->allignAdd = (aliSize-1);
lm->allignMask = ~lm->allignAdd;
if (mem != NULL)
    {
    lm->doMemoryAllocs = FALSE;
    struct lmBlock *mb = mem;
    mb->free = (char *)(mb+1);
    mb->end = ((char *)mb) + blockSize;
    mb->next = lm->blocks;
    lm->blocks = mb;
    }
else
    {
    lm->doMemoryAllocs = TRUE;
    newBlock(lm, blockSize);
    }

return lm;
}

struct lm *lmInit(int blockSize)
/* Create a local memory pool. */
{
return lmGuts(blockSize, NULL);
}

struct lm *lmInitWMem(void *mem, int blockSize)
/* Create a local memory pool. */
{
return lmGuts(blockSize, mem);
}

void lmCleanup(struct lm **pLm)
/* Clean up a local memory pool. */
{
    struct lm *lm = *pLm;
    if (lm == NULL)
        return;
    *pLm = NULL;
    slFreeList(&lm->blocks);
    freeMem(lm);
}

unsigned int lmBlockHeaderSize()
// Return the size of an lmBlock.
{
return sizeof(struct lmBlock);
}

size_t lmUsed(struct lm *lm)
// Returns amount of memory allocated
{
size_t used = 0;

struct lmBlock *mb = lm->blocks;
for (;mb != NULL;mb = mb->next)
    used += (mb->free - (char *)(mb+1));

return used;
}

size_t lmAvailable(struct lm *lm)
// Returns currently available memory in pool
{
struct lmBlock *mb = lm->blocks;
return (mb->end - mb->free);
}

size_t lmSize(struct lm *lm)
// Returns current size of pool, even for memory already allocated
{
size_t fullSize = 0;

struct lmBlock *mb = lm->blocks;
for (;mb != NULL;mb = mb->next)
    fullSize += (mb->end - (char *)(mb+1));

return fullSize;
}

void *lmAlloc(struct lm *lm, size_t size)
/* Allocate memory from local pool. */
{
struct lmBlock *mb = lm->blocks;
void *ret;
size_t memLeft = mb->end - mb->free;
if (memLeft < size)
    {
    if (lm->doMemoryAllocs)
        mb = newBlock(lm, size);
    else
        errAbort("attempted local memory alloc in fixed size allocator");
    }

ret = mb->free;
mb->free += ((size+lm->allignAdd)&lm->allignMask);
if (mb->free > mb->end)
    mb->free = mb->end;
return ret;
}

void *lmAllocMoreMem(struct lm *lm, void *pt, size_t oldSize, size_t newSize)
/* Adjust memory size on a block, possibly relocating it.  If block is grown,
 * new memory is zeroed. */
{
struct lmBlock *mb = lm->blocks;
// rare case that pointer is to last lm alloc, but still try.
// Note this is the one place where the pointer gets reused and it is known to be in this lm
if ((char *)pt + oldSize == mb->free
&&  (char *)pt + newSize <= mb->end)
    {
    if (newSize > oldSize) // only move the free pointer on more mem
        mb->free = pt + newSize;
    return pt;
    }
void *new = lmAlloc(lm, newSize);
memcpy(new, pt, oldSize);
return new;
}

void *lmCloneMem(struct lm *lm, void *pt, size_t size)
/* Return a local mem copy of memory block. */
{
void *d = lmAlloc(lm, size);
memcpy(d, pt, size);
return d;
}

char *lmCloneStringZ(struct lm *lm, const char *string, int size)
/* Return local mem copy of string. */
{
if (string == NULL)
    return NULL;
else
    {
    char *s = lmAlloc(lm, size+1);
    memcpy(s, string, size);
    return s;
    }
}

char *lmCloneString(struct lm *lm, const char *string)
/* Return local mem copy of string. */
{
if (string == NULL)
    return NULL;
else
    return lmCloneStringZ(lm, string, strlen(string));
}

char *lmCloneFirstWord(struct lm *lm, const char *line)
/* Clone first word in line */
{
char *startFirstWord = skipLeadingSpaces(line);
if (startFirstWord == NULL)
    return NULL;
char *endFirstWord = skipToSpaces(startFirstWord);
if (endFirstWord == NULL)
    return lmCloneString(lm, startFirstWord);
else
    return lmCloneStringZ(lm, startFirstWord, endFirstWord - startFirstWord);
}
    
char *lmCloneSomeWord(struct lm *lm, const char *line, int wordIx)
/* Return a clone of the given space-delimited word within line.  Returns NULL if
 * not that many words in line. */
{
if (wordIx < 0)
    return NULL;
int i;
for (i=0; i<wordIx; ++i)
    {
    line = skipLeadingSpaces(line);
    line = skipToSpaces(line);
    if (line == NULL)
        return NULL;
    }
return lmCloneFirstWord(lm, line);
}


struct slName *lmSlName(struct lm *lm, const char *name)
/* Return slName in memory. */
{
struct slName *n;
int size = sizeof(*n) + strlen(name) + 1;
n = lmAlloc(lm, size);
strcpy(n->name, name);
return n;
}

char **lmCloneRowExt(struct lm *lm, char **row, int rowOutSize, int rowInSize)
/* Allocate an array of strings with rowOutSize elements.  Clone the first rowInSize elements
 * of row into the new array, leaving others NULL if rowOutSize is greater than rowInSize.
 * rowOutSize must be greater than or equal to rowInSize. */
{
if (rowOutSize < rowInSize)
    errAbort("lmCloneRowExt: rowOutSize (%d) must be greater than or equal to rowInSize (%d)",
             rowOutSize, rowInSize);
char **rowClone = NULL;
lmAllocArray(lm, rowClone, rowOutSize);
int i;
for (i = 0;  i < rowInSize;  i++)
    rowClone[i] = lmCloneString(lm, row[i]);
return rowClone;
}

char **lmCloneRow(struct lm *lm, char **row, int rowSize)
/* Allocate an array of strings and its contents cloned from row. */
{
return lmCloneRowExt(lm, row, rowSize, rowSize);
}

void lmRefAdd(struct lm *lm, struct slRef **pRefList, void *val)
/* Add reference to list. */
{
struct slRef *ref;
lmAllocVar(lm, ref);
ref->val = val;
slAddHead(pRefList, ref);
}
