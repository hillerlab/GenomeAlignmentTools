/* Copyright (C) 2011 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

/* errCatch - help catch errors so that errAborts aren't
 * fatal, and warn's don't necessarily get printed immediately. 
 * Note that error conditions caught this way will tend to
 * leak resources unless there are additional wrappers. 
 *
 * Typical usage is
 * errCatch = errCatchNew();
 * if (errCatchStart(errCatch))
 *     doFlakyStuff();
 * errCatchEnd(errCatch);
 * if (errCatch->gotError)
 *     warn("%s", errCatch->message->string);
 * errCatchFree(&errCatch); 
 * cleanupFlakyStuff();
 */

#include "common.h"
#include "errAbort.h"
#include "dystring.h"
#include "hash.h"
#include <pthread.h>
#include "errCatch.h"



struct errCatch *errCatchNew()
/* Return new error catching structure. */
{
struct errCatch *errCatch;
AllocVar(errCatch);
errCatch->message = dyStringNew(0);
return errCatch;
}

void errCatchFree(struct errCatch **pErrCatch)
/* Free up resources associated with errCatch */
{
struct errCatch *errCatch = *pErrCatch;
if (errCatch != NULL)
    {
    dyStringFree(&errCatch->message);
    freez(pErrCatch);
    }
}

static struct errCatch **getStack()
/* Return a pointer to the errCatch object stack for the current pthread. */
{
static pthread_mutex_t getStackMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_lock( &getStackMutex );
static struct hash *perThreadStacks = NULL;
pthread_t pid = pthread_self(); //  can be a pointer or a number
// A true integer has function would be nicer, but this will do.  
// Don't safef, theoretically that could abort.
char key[64];
snprintf(key, sizeof(key), "%lld",  ptrToLL(pid));
key[ArraySize(key)-1] = '\0';
if (perThreadStacks == NULL)
    perThreadStacks = hashNew(0);
struct hashEl *hel = hashLookup(perThreadStacks, key);
if (hel == NULL)
    hel = hashAdd(perThreadStacks, key, NULL);
pthread_mutex_unlock( &getStackMutex );
return (struct errCatch **)(&hel->val);
}

static void errCatchAbortHandler()
/* semiAbort */
{
struct errCatch **pErrCatchStack = getStack(), *errCatchStack = *pErrCatchStack;
errCatchStack->gotError = TRUE;
longjmp(errCatchStack->jmpBuf, -1);
}

static void errCatchWarnHandler(char *format, va_list args)
/* Write an error to top of errCatchStack. */
{
struct errCatch **pErrCatchStack = getStack(), *errCatchStack = *pErrCatchStack;
errCatchStack->gotWarning = TRUE;
dyStringVaPrintf(errCatchStack->message, format, args);
dyStringAppendC(errCatchStack->message, '\n');
}

boolean errCatchPushHandlers(struct errCatch *errCatch)
/* Push error handlers.  Not usually called directly. */
{
pushAbortHandler(errCatchAbortHandler);
pushWarnHandler(errCatchWarnHandler);
struct errCatch **pErrCatchStack = getStack();
slAddHead(pErrCatchStack, errCatch);
return TRUE;
}

void errCatchEnd(struct errCatch *errCatch)
/* Restore error handlers and pop self off of catching stack. */
{
popWarnHandler();
popAbortHandler();
struct errCatch **pErrCatchStack = getStack(), *errCatchStack = *pErrCatchStack;
if (errCatch != errCatchStack)
   errAbort("Mismatch between errCatch and errCatchStack");
*pErrCatchStack = errCatch->next;
}

void errCatchReWarn(struct errCatch *errCatch)
/* Re-warn any warnings that happened even though no abort happened 
 * to make them visible. */
{
if (errCatch->gotWarning && !errCatch->gotError)
    {
    warn("%s", errCatch->message->string);
    }
}

boolean errCatchFinish(struct errCatch **pErrCatch)
/* Finish up error catching.  Report error if there is a
 * problem and return FALSE.  If no problem return TRUE.
 * This handles errCatchEnd and errCatchFree. */
{
struct errCatch *errCatch = *pErrCatch;
boolean ok = TRUE;
if (errCatch != NULL)
    {
    errCatchEnd(errCatch);
    if (errCatch->gotError)
	{
	ok = FALSE;
	warn("%s", errCatch->message->string);
	}
    else
	errCatchReWarn(errCatch);
    errCatchFree(pErrCatch);
    }
return ok;
}

