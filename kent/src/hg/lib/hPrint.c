/* hPrint - turning html printing on and off, which is useful
 * when postscript and PDF images are being drawn  */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "hPrint.h"
#include "htmshell.h"


static boolean suppressHtml = FALSE;
/* If doing PostScript output we'll suppress most of HTML output. */

boolean hPrintStatus()
/* is html printing on or off ?
   return TRUE for print is on, FALSE for printing is off */
{
return ! suppressHtml;
}

void hPrintDisable()
/* turn html printing off */
{
suppressHtml = TRUE;
}

void hPrintEnable()
/* turn html printing on */
{
suppressHtml = FALSE;
}

void hvPrintf(char *format, va_list args)
/* Suppressable variable args printf. Check for write error so we can
 * terminate if http connection breaks. */
{
if (suppressHtml)
    return;
vprintf(format, args);
if (ferror(stdout))
    noWarnAbort();
}

void hPrintf(char *format, ...)
/* Printf that can be suppressed if not making html. */
{
va_list(args);
va_start(args, format);
hvPrintf(format, args);
va_end(args);
}

void hPrintNonBreak(char *s)
/* Print out string but replace spaces with &nbsp; */
{
char c;

if (suppressHtml)
    return;
while ((c = *s++) != '\0')
    if (c == ' ')
	fputs("&nbsp;", stdout);
    else
        putchar(c);
}

void hPrintEncodedNonBreak(char *s)
/* Print with htmlEncode and non-break */
{
char *encoded = htmlEncode(s);
hPrintNonBreak(encoded);
freeMem(encoded);
}

void hPuts(char *string)
/* Puts that can be suppressed if not making html. */
{
if (!suppressHtml)
    puts(string);
}

void hPutc(char c)
/* putc that can be suppressed if not making html. */
{
if (!suppressHtml)
    fputc(c, stdout);
}

void hWrites(char *string)
/* Write string with no '\n' if not suppressed. */
{
if (!suppressHtml)
    fputs(string, stdout);
}

void hButtonMaybePressed(char *name, char *label, char *msg, char *onClick, boolean pressed)
/* If not suppresed, write out button optionally with onclick javascript, message and 
   styled to indicate modal state (button pressed)
 */
{
if (!suppressHtml)
    cgiMakeSubmitButtonMaybePressed(name, label, msg, onClick, pressed);
}

void hButton(char *name, char *label)
/* Write out button if not suppressed. */
{
hButtonMaybePressed(name, label, NULL, NULL, FALSE);
}

void hButtonWithMsg(char *name, char *label, char *msg)
/* Write out button with msg if not suppressed. */
{
hButtonMaybePressed(name, label, msg, NULL, FALSE);
}

void hButtonWithOnClick(char *name, char *label, char *msg, char *onClick)
/* Write out button with onclick javascript if not suppressed. */
{
hButtonMaybePressed(name, label, msg, onClick, FALSE);
}

void hOnClickButton(char *id, char *command, char *label)
/* Write out button with onclick command */
{
hButtonMaybePressed(id, label, NULL, command, FALSE);
}

void hTextVar(char *varName, char *initialVal, int charSize)
/* Write out text entry field if not suppressed. */
{
if (!suppressHtml)
    cgiMakeTextVar(varName, initialVal, charSize);
}

void hIntVar(char *varName, int initialVal, int maxDigits)
/* Write out numerical entry field if not supressed. */
{
if (!suppressHtml)
    cgiMakeIntVar(varName, initialVal, maxDigits);
}

void hDoubleVar(char *varName, double initialVal, int maxDigits)
/* Write out numerical entry field if not supressed. */
{
if (!suppressHtml)
    cgiMakeDoubleVar(varName, initialVal, maxDigits);
}

void hCheckBox(char *varName, boolean checked)
/* Make check box if not suppressed. */
{
if (!suppressHtml)
    cgiMakeCheckBox(varName, checked);
}

void hDropListClassWithStyle(char *name, char *menu[], int menuSize, 
                                char *checked, char *class, char *style)
/* Make a drop-down list with names if not suppressed, 
 * using specified class and style */
{
if (!suppressHtml)
    cgiMakeDropListClassWithStyle(name, menu, menuSize, checked, class, style);
}

void hDropListClass(char *name, char *menu[], int menuSize, char *checked,
                        char *class)
/* Make a drop-down list with names if not suppressed, using specified class. */
{
hDropListClassWithStyle(name, menu, menuSize, checked, class, NULL);
}

void hDropList(char *name, char *menu[], int menuSize, char *checked)
/* Make a drop-down list with names if not suppressed. */
{
hDropListClass(name, menu, menuSize, checked, NULL);
}

void hPrintComment(char *format, ...)
/* Function to print output as a comment so it is not seen in the HTML
 * output but only in the HTML source. */
{
va_list(args);
va_start(args, format);
hWrites("\n<!-- DEBUG: ");
hvPrintf(format, args);
hWrites(" -->\n");
fflush(stdout); /* USED ONLY FOR DEBUGGING BECAUSE THIS IS SLOW - MATT */
va_end(args);
}
