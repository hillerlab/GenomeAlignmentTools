/* mouseOrtho.c was originally generated by the autoSql program, which also 
 * generated mouseOrtho.h and mouseOrtho.sql.  This module links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "mouseOrtho.h"


void mouseOrthoStaticLoad(char **row, struct mouseOrtho *ret)
/* Load a row from mouseOrtho table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->score = sqlUnsigned(row[4]);
ret->strand = row[5];
ret->thickStart = sqlUnsigned(row[6]);
ret->thickEnd = sqlUnsigned(row[7]);
}

struct mouseOrtho *mouseOrthoLoad(char **row)
/* Load a mouseOrtho from row fetched with select * from mouseOrtho
 * from database.  Dispose of this with mouseOrthoFree(). */
{
struct mouseOrtho *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
ret->strand = cloneString(row[5]);
ret->thickStart = sqlUnsigned(row[6]);
ret->thickEnd = sqlUnsigned(row[7]);
return ret;
}

struct mouseOrtho *mouseOrthoLoadAll(char *fileName) 
/* Load all mouseOrtho from a tab-separated file.
 * Dispose of this with mouseOrthoFreeList(). */
{
struct mouseOrtho *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[8];

while (lineFileRow(lf, row))
    {
    el = mouseOrthoLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct mouseOrtho *mouseOrthoLoadWhere(struct sqlConnection *conn, char *table, char *where)
/* Load all mouseOrtho from table that satisfy where clause. The
 * where clause may be NULL in which case whole table is loaded
 * Dispose of this with mouseOrthoFreeList(). */
{
struct mouseOrtho *list = NULL, *el;
struct dyString *query = dyStringNew(256);
struct sqlResult *sr;
char **row;

sqlDyStringPrintf(query, "select * from %s", table);
if (where != NULL)
    dyStringPrintf(query, " where %s", where);
sr = sqlGetResult(conn, query->string);
while ((row = sqlNextRow(sr)) != NULL)
    {
    el = mouseOrthoLoad(row);
    slAddHead(&list, el);
    }
slReverse(&list);
sqlFreeResult(&sr);
dyStringFree(&query);
return list;
}

struct mouseOrtho *mouseOrthoCommaIn(char **pS, struct mouseOrtho *ret)
/* Create a mouseOrtho out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new mouseOrtho */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->score = sqlUnsignedComma(&s);
ret->strand = sqlStringComma(&s);
ret->thickStart = sqlUnsignedComma(&s);
ret->thickEnd = sqlUnsignedComma(&s);
*pS = s;
return ret;
}

void mouseOrthoFree(struct mouseOrtho **pEl)
/* Free a single dynamically allocated mouseOrtho such as created
 * with mouseOrthoLoad(). */
{
struct mouseOrtho *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->strand);
freez(pEl);
}

void mouseOrthoFreeList(struct mouseOrtho **pList)
/* Free a list of dynamically allocated mouseOrtho's */
{
struct mouseOrtho *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    mouseOrthoFree(&el);
    }
*pList = NULL;
}

void mouseOrthoOutput(struct mouseOrtho *el, FILE *f, char sep, char lastSep) 
/* Print out mouseOrtho.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->chromStart);
fputc(sep,f);
fprintf(f, "%u", el->chromEnd);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->score);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->strand);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->thickStart);
fputc(sep,f);
fprintf(f, "%u", el->thickEnd);
fputc(lastSep,f);
}
