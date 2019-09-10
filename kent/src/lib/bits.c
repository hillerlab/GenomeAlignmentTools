/* bits - handle operations on arrays of bits. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include "bits.h"



static Bits oneBit[8] = { 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1};
static Bits leftMask[8] = {0xFF, 0x7F, 0x3F, 0x1F,  0xF,  0x7,  0x3,  0x1,};
static Bits rightMask[8] = {0x80, 0xC0, 0xE0, 0xF0, 0xF8, 0xFC, 0xFE, 0xFF,};
int bitsInByte[256];

static boolean inittedBitsInByte = FALSE;

void bitsInByteInit()
/* Initialize bitsInByte array. */
{
int i;

if (!inittedBitsInByte)
    {
    inittedBitsInByte = TRUE;
    for (i=0; i<256; ++i)
        {
	int count = 0;
	if (i&1)
	    count = 1;
	if (i&2)
	    ++count;
	if (i&4)
	    ++count;
	if (i&8)
	    ++count;
	if (i&0x10)
	    ++count;
	if (i&0x20)
	    ++count;
	if (i&0x40)
	    ++count;
	if (i&0x80)
	    ++count;
	bitsInByte[i] = count;
	}
    }
}

Bits *bitAlloc(int bitCount)
/* Allocate bits. */
{
int byteCount = ((bitCount+7)>>3);
return needLargeZeroedMem(byteCount);
}

Bits *bitRealloc(Bits *b, int bitCount, int newBitCount)
/* Resize a bit array.  If b is null, allocate a new array */
{
int byteCount = ((bitCount+7)>>3);
int newByteCount = ((newBitCount+7)>>3);
return needLargeZeroedMemResize(b, byteCount, newByteCount);
}

Bits *bitClone(Bits* orig, int bitCount)
/* Clone bits. */
{
int byteCount = ((bitCount+7)>>3);
Bits* bits = needLargeZeroedMem(byteCount);
memcpy(bits, orig, byteCount);
return bits;
}

void bitFree(Bits **pB)
/* Free bits. */
{
freez(pB);
}

Bits *lmBitAlloc(struct lm *lm,int bitCount)
// Allocate bits.  Must supply local memory.
{
assert(lm != NULL);
int byteCount = ((bitCount+7)>>3);
return lmAlloc(lm,byteCount);
}

Bits *lmBitRealloc(struct lm *lm,Bits *b, int bitCount, int newBitCount)
// Resize a bit array.  If b is null, allocate a new array.  Must supply local memory.
{
assert(lm != NULL);
int byteCount = ((bitCount+7)>>3);
int newByteCount = ((newBitCount+7)>>3);
return lmAllocMoreMem(lm, b ,byteCount, newByteCount);
}

Bits *lmBitClone(struct lm *lm,Bits* orig, int bitCount)
// Clone bits.  Must supply local memory.
{
assert(lm != NULL);
int byteCount = ((bitCount+7)>>3);
Bits* bits = lmAlloc(lm,byteCount);
memcpy(bits, orig, byteCount);
return bits;
}

void bitSetOne(Bits *b, int bitIx)
/* Set a single bit. */
{
b[bitIx>>3] |= oneBit[bitIx&7];
}

void bitClearOne(Bits *b, int bitIx)
/* Clear a single bit. */
{
b[bitIx>>3] &= ~oneBit[bitIx&7];
}

void bitSetRange(Bits *b, int startIx, int bitCount)
/* Set a range of bits. */
{
if (bitCount <= 0)
    return;
int endIx = (startIx + bitCount - 1);
int startByte = (startIx>>3);
int endByte = (endIx>>3);
int startBits = (startIx&7);
int endBits = (endIx&7);
int i;

if (startByte == endByte)
    {
    b[startByte] |= (leftMask[startBits] & rightMask[endBits]);
    return;
    }
b[startByte] |= leftMask[startBits];
for (i = startByte+1; i<endByte; ++i)
    b[i] = 0xff;
b[endByte] |= rightMask[endBits];
}


boolean bitReadOne(Bits *b, int bitIx)
/* Read a single bit. */
{
return (b[bitIx>>3] & oneBit[bitIx&7]) != 0;
}

int bitCountRange(Bits *b, int startIx, int bitCount)
/* Count number of bits set in range. */
{
if (bitCount <= 0)
    return 0;
int endIx = (startIx + bitCount - 1);
int startByte = (startIx>>3);
int endByte = (endIx>>3);
int startBits = (startIx&7);
int endBits = (endIx&7);
int i;
int count = 0;

if (!inittedBitsInByte)
    bitsInByteInit();
if (startByte == endByte)
    return bitsInByte[b[startByte] & leftMask[startBits] & rightMask[endBits]];
count = bitsInByte[b[startByte] & leftMask[startBits]];
for (i = startByte+1; i<endByte; ++i)
    count += bitsInByte[b[i]];
count += bitsInByte[b[endByte] & rightMask[endBits]];
return count;
}

int bitFind(Bits *b, int startIx, boolean val, int bitCount)
/* Find the index of the the next set bit. */
{
unsigned char notByteVal = val ? 0 : 0xff;
int iBit = startIx;
int endByte = ((bitCount-1)>>3);
int iByte;

/* scan initial byte */
while (((iBit & 7) != 0) && (iBit < bitCount))
    {
    if (bitReadOne(b, iBit) == val)
        return iBit;
    iBit++;
    }

/* scan byte at a time, if not already in last byte */
iByte = (iBit >> 3);
if (iByte < endByte)
    {
    while ((iByte < endByte) && (b[iByte] == notByteVal))
        iByte++;
    iBit = iByte << 3;
    }

/* scan last byte */
while (iBit < bitCount)
    {
    if (bitReadOne(b, iBit) == val)
        return iBit;
    iBit++;
    }
 return bitCount;  /* not found */
}

int bitFindSet(Bits *b, int startIx, int bitCount)
/* Find the index of the the next set bit. */
{
return bitFind(b, startIx, TRUE, bitCount);
}

int bitFindClear(Bits *b, int startIx, int bitCount)
/* Find the index of the the next clear bit. */
{
return bitFind(b, startIx, FALSE, bitCount);
}

void bitClear(Bits *b, int bitCount)
/* Clear many bits (possibly up to 7 beyond bitCount). */
{
int byteCount = ((bitCount+7)>>3);
zeroBytes(b, byteCount);
}

void bitClearRange(Bits *b, int startIx, int bitCount)
/* Clear a range of bits. */
{
if (bitCount <= 0)
    return;
int endIx = (startIx + bitCount - 1);
int startByte = (startIx>>3);
int endByte = (endIx>>3);
int startBits = (startIx&7);
int endBits = (endIx&7);
int i;

if (startByte == endByte)
    {
    b[startByte] &= ~(leftMask[startBits] & rightMask[endBits]);
    return;
    }
b[startByte] &= ~leftMask[startBits];
for (i = startByte+1; i<endByte; ++i)
    b[i] = 0x00;
b[endByte] &= ~rightMask[endBits];
}

void bitAnd(Bits *a, Bits *b, int bitCount)
/* And two bitmaps.  Put result in a. */
{
int byteCount = ((bitCount+7)>>3);
while (--byteCount >= 0)
    {
    *a = (*a & *b++);
    a++;
    }
}

int bitAndCount(Bits *a, Bits *b, int bitCount)
// Without altering 2 bitmaps, count the AND bits.
{
int byteCount = ((bitCount+7)>>3);
int count = 0;
if (!inittedBitsInByte)
    bitsInByteInit();
while (--byteCount >= 0)
    count += bitsInByte[(*a++ & *b++)];

return count;
}

void bitOr(Bits *a, Bits *b, int bitCount)
/* Or two bitmaps.  Put result in a. */
{
int byteCount = ((bitCount+7)>>3);
while (--byteCount >= 0)
    {
    *a = (*a | *b++);
    a++;
    }
}

int bitOrCount(Bits *a, Bits *b, int bitCount)
// Without altering 2 bitmaps, count the OR'd bits.
{
int byteCount = ((bitCount+7)>>3);
int count = 0;
if (!inittedBitsInByte)
    bitsInByteInit();
while (--byteCount >= 0)
    count += bitsInByte[(*a++ | *b++)];

return count;
}

void bitXor(Bits *a, Bits *b, int bitCount)
{
int byteCount = ((bitCount+7)>>3);
while (--byteCount >= 0)
    {
    *a = (*a ^ *b++);
    a++;
    }
}

int bitXorCount(Bits *a, Bits *b, int bitCount)
// Without altering 2 bitmaps, count the XOR'd bits.
{
int byteCount = ((bitCount+7)>>3);
int count = 0;
if (!inittedBitsInByte)
    bitsInByteInit();
while (--byteCount >= 0)
    count += bitsInByte[(*a++ ^ *b++)];

return count;
}

void bitNot(Bits *a, int bitCount)
/* Flip all bits in a. */
{
int byteCount = ((bitCount+7)>>3);
while (--byteCount >= 0)
    {
    *a = ~*a;
    a++;
    }
}

void bitReverseRange(Bits *bits, int startIx, int bitCount)
// Reverses bits in range (e.g. 110010 becomes 010011)
{
if (bitCount <= 0)
    return;
int ixA = startIx;
int ixB = (startIx + bitCount - 1);
for ( ;ixA < ixB; ixA++, ixB--)
    {
    boolean bitA = bitReadOne(bits, ixA);
    boolean bitB = bitReadOne(bits, ixB);
    if (!bitA && bitB)
        {
        bitSetOne(  bits, ixA);
        bitClearOne(bits, ixB);
        }
    else if (bitA && !bitB)
        {
        bitClearOne(bits, ixA);
        bitSetOne(  bits, ixB);
        }
    }
}


void bitPrint(Bits *a, int startIx, int bitCount, FILE* out)
/* Print part or all of bit map as a string of 0s and 1s.  Mostly useful for
 * debugging */
{
int i;
for (i = startIx; i < bitCount; i++)
    {
    if (bitReadOne(a, i))
        fputc('1', out);
    else
        fputc('0', out);
    }
fputc('\n', out);
}

void bitsOut(FILE* out, Bits *bits, int startIx, int bitCount, boolean onlyOnes)
// Print part or all of bit map as a string of 0s and 1s.
// If onlyOnes, enclose result in [] and use ' ' instead of '0'.
{
if (onlyOnes)
    fputc('[', out);

int ix = startIx;
for ( ; ix < bitCount; ix++)
    {
    if (bitReadOne(bits, ix))
        fputc('1', out);
    else
        {
        if (onlyOnes)
            fputc(' ', out);
        else
            fputc('0', out);
        }
    }
if (onlyOnes)
    fputc(']', out);
}

Bits *bitsIn(struct lm *lm,char *bitString, int len)
// Returns a bitmap from a string of 1s and 0s.  Any non-zero, non-blank char sets a bit.
// Returned bitmap is the size of len even if that is longer than the string.
// Optionally supply local memory.  Note does NOT handle enclosing []s printed with bitsOut().
{
if (bitString == NULL || len == 0)
    return NULL;

Bits *bits = NULL;
if (lm != NULL)
    bits = lmBitAlloc(lm,len);
else
    bits = bitAlloc(len);

int ix = 0;
for ( ;ix < len && bitString[ix] != '\0'; ix++)
    {
    if (bitString[ix] != '0' && bitString[ix] != ' ')
        bitSetOne(bits, ix);
    }
return bits;
}

