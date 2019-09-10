/* Memgfx - stuff to do graphics in memory buffers.
 * Typically will just write these out as .gif or .png files.
 * This stuff is byte-a-pixel for simplicity.
 * It can do 256 colors.
 *
 * This file is copyright 2000 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef MEMGFX_H
#define MEMGFX_H

#ifndef GFXPOLY_H
#include "gfxPoly.h"
#endif

typedef unsigned int Color;

#if defined(__sgi__) || defined(__sgi) || defined(__powerpc__) || defined(sparc) || defined(__ppc__) || defined(__s390__) || defined(__s390x__)

// BIGENDIAN machines:

#define MEMGFX_BIGENDIAN	1
#define MG_WHITE   0xffffffff
#define MG_BLACK   0x000000ff
#define MG_RED     0xff0000ff
#define MG_GREEN   0x00ff00ff
#define MG_BLUE    0x0000ffff
#define MG_CYAN    0x00ffffff
#define MG_MAGENTA 0xff00ffff
#define MG_YELLOW  0xffff00ff
#define MG_GRAY    0x808080ff

#define MAKECOLOR_32(r,g,b) (((unsigned int)0xff) | ((unsigned int)b<<8) | ((unsigned int)g << 16) | ((unsigned int)r << 24))
#define COLOR_32_RED(c) (((c)>>24)&0xff)
#define COLOR_32_GREEN(c) (((c)>>16)&0xff)
#define COLOR_32_BLUE(c) (((c)>>8)&0xff)

#else

// LITTLE ENDIAN machines:

#define MG_WHITE   0xffffffff
#define MG_BLACK   0xff000000
#define MG_RED     0xff0000ff
#define MG_GREEN   0xff00ff00
#define MG_BLUE    0xffff0000
#define MG_CYAN    0xffffff00
#define MG_MAGENTA 0xffff00ff
#define MG_YELLOW  0xff00ffff
#define MG_GRAY    0xff808080

#define MAKECOLOR_32(r,g,b) (((unsigned int)0xff<<24) | ((unsigned int)b<<16) | ((unsigned int)g << 8) | (unsigned int)r)
#define COLOR_32_RED(c) ((c)&0xff)
#define COLOR_32_GREEN(c) (((c)>>8)&0xff)
#define COLOR_32_BLUE(c) (((c)>>16)&0xff)
#endif

#define MG_WRITE_MODE_NORMAL    0
#define MG_WRITE_MODE_MULTIPLY  (1 << 0)

struct rgbColor
    {
    unsigned char r, g, b;
    };

/* HSV and HSL structs can be used for changing lightness, darkness, or
 * color of RGB colors. Convert RGB->HS[LV], modify hue, saturation, or
 * value/lightness, then convert back to RGB.
 * The datatypes were chosen to be fast but also give accurate conversion
 * back to RGB.
 * Hue is a float [0,360) degrees 0=red, 120=green, 240=blue
 * S/V/L are integers [0,1000]
 */

struct hsvColor
    {
    double h;
    unsigned short s, v;
    };

struct hslColor
    {
    double h;
    unsigned short s, l;
    };

extern struct rgbColor mgFixedColors[9];  /* Contains MG_WHITE - MG_GRAY */

struct memGfx
    {
    Color *pixels;
    int width, height;
    struct rgbColor colorMap[256];
    int colorsUsed;
    int clipMinX, clipMaxX;
    int clipMinY, clipMaxY;
    struct colHash *colorHash;	/* Hash for fast look up of color. */
    unsigned int writeMode;
    };

struct memGfx *mgNew(int width, int height);
/* Get a new thing to draw on in memory. */

void mgFree(struct memGfx **pmg);
/* Free up memory raster. */

void mgClearPixelsTrans(struct memGfx *mg);
/* Set all pixels to transparent. */

void mgClearPixels(struct memGfx *mg);
/* Set all pixels to background. */

void mgSetClip(struct memGfx *mg, int x, int y, int width, int height);
/* Set clipping rectangle. */

void mgUnclip(struct memGfx *mg);
/* Set clipping rect cover full thing. */

Color mgFindColor(struct memGfx *mg, unsigned char r, unsigned char g, unsigned char b);
/* Returns closest color in color map to rgb values.  If it doesn't
 * already exist in color map and there's room, it will create
 * exact color in map. */

Color mgClosestColor(struct memGfx *mg, unsigned char r, unsigned char g, unsigned char b);
/* Returns closest color in color map to r,g,b */

Color mgAddColor(struct memGfx *mg, unsigned char r, unsigned char g, unsigned char b);
/* Adds color to end of color map if there's room. */

int mgColorsFree(struct memGfx *mg);
/* Returns # of unused colors in color map. */


#define _mgBpr(mg) ((mg)->width)
/* Get what to add to get to next line */

#define _mgPixAdr(mg,x,y) ((mg)->pixels+_mgBpr(mg) * (y) + (x))
/* Get pixel address */

#define _mgPutDot(mg, x, y, color) (*_mgPixAdr(mg,x,y) = (color))
/* Unclipped plot a dot */

#define _mgGetDot(mg, x, y) (*_mgPixAdr(mg,x,y))
/* Unclipped get a dot, you do not want to use this, this is special for
 * verticalText only */

void _mgPutDotMultiply(struct memGfx *mg, int x, int y,Color color);

INLINE void mgPutDot(struct memGfx *mg, int x, int y,Color color)
{
if ((x)>=(mg)->clipMinX && (x) < (mg)->clipMaxX && (y)>=(mg)->clipMinY  && (y) < (mg)->clipMaxY) 
    {
        switch(mg->writeMode)
        {
        case MG_WRITE_MODE_NORMAL:
            {
            _mgPutDot(mg,x,y,color);
            }
            break;
        case MG_WRITE_MODE_MULTIPLY:
            {
            _mgPutDotMultiply(mg,x,y,color);
            }
            break;
        }
    }
}
/* Clipped put dot */

#define mgGetDot(mg,x,y) ((x)>=(mg)->clipMinX && (x) < (mg)->clipMaxX && (y)>=(mg)->clipMinY  && (y) < (mg)->clipMaxY) ? _mgGetDot(mg,x,y) : 0
/* Clipped get dot, you do not want to use this, this is special for
 * verticalText only */

void mgPutSeg(struct memGfx *mg, int x, int y, int width, Color *dots);
/* Put a series of dots starting at x, y and going to right width pixels. */

void mgPutSegZeroClear(struct memGfx *mg, int x, int y, int width, Color *dots);
/* Put a series of dots starting at x, y and going to right width pixels.
 * Don't put zero dots though. */

void mgDrawBox(struct memGfx *mg, int x, int y, int width, int height, Color color);
/* Draw a (horizontal) box */

void mgDrawLine(struct memGfx *mg, int x1, int y1, int x2, int y2, Color color);
/* Draw a line from one point to another. */

void mgDrawHorizontalLine(struct memGfx *mg, int y1, Color color);
/*special case of mgDrawLine*/

void mgLineH(struct memGfx *mg, int y, int x1, int x2, Color color);
/* Draw horizizontal line width pixels long starting at x/y in color */

void mgSavePng(struct memGfx *mg, char *filename, boolean useTransparency);
/* Save memory bitmap to filename as a PNG.
 * If useTransparency, then the first color in memgfx's colormap/palette is
 * assumed to be the image background color, and pixels of that color
 * are made transparent. */

boolean mgSaveToPng(FILE *png_file, struct memGfx *mg, boolean useTransparency);
/* Save PNG to an already open file.
 * If useTransparency, then the first color in memgfx's colormap/palette is
 * assumed to be the image background color, and pixels of that color
 * are made transparent. */

typedef void (*TextBlit)(int bitWidth, int bitHeight, int bitX, int bitY,
	unsigned char *bitData, int bitDataRowBytes, 
	struct memGfx *dest, int destX, int destY, 
	Color color, Color backgroundColor);
/* This defines the type of a function that takes a rectangular
 * area of a bitplane and expands it into a rectangular area
 * of a full color screen. */

void mgTextBlit(int bitWidth, int bitHeight, int bitX, int bitY,
	unsigned char *bitData, int bitDataRowBytes, 
	struct memGfx *dest, int destX, int destY, 
	Color color, Color backgroundColor);
/* This function leaves the background as it was. */

void mgTextBlitSolid(int bitWidth, int bitHeight, int bitX, int bitY,
	unsigned char *bitData, int bitDataRowBytes, 
	struct memGfx *dest, int destX, int destY, 
	Color color, Color backgroundColor);
/* This function sets the background to the background color. */

typedef struct font_hdr MgFont;
/* Type of our font.  */

/* Collection of fonts from here and there.  The mgTinyFont() and mgSmallFont() are uniq
 * here.  The rest are synonyms at this point for the adobe fonts below. */
MgFont *mgTinyFont();
MgFont *mgSmallFont();
MgFont *mgMediumFont();
MgFont *mgLargeFont();
MgFont *mgHugeFont();
MgFont *mgTinyBoldFont();
MgFont *mgSmallBoldFont();
MgFont *mgMediumBoldFont();
MgFont *mgLargeBoldFont();
MgFont *mgHugeBoldFont();
MgFont *mgTinyFixedFont();
MgFont *mgSmallFixedFont();
MgFont *mgMediumFixedFont();
MgFont *mgLargeFixedFont();
MgFont *mgHugeFixedFont();

/* Adobe fonts from xfree project. */
MgFont *mgCourier8Font();
MgFont *mgCourier10Font();
MgFont *mgCourier12Font();
MgFont *mgCourier14Font();
MgFont *mgCourier18Font();
MgFont *mgCourier24Font();
MgFont *mgCourier34Font();
MgFont *mgHelvetica8Font();
MgFont *mgHelvetica10Font();
MgFont *mgHelvetica12Font();
MgFont *mgHelvetica14Font();
MgFont *mgHelvetica18Font();
MgFont *mgHelvetica24Font();
MgFont *mgHelvetica34Font();
MgFont *mgHelveticaBold8Font();
MgFont *mgHelveticaBold10Font();
MgFont *mgHelveticaBold12Font();
MgFont *mgHelveticaBold14Font();
MgFont *mgHelveticaBold18Font();
MgFont *mgHelveticaBold24Font();
MgFont *mgHelveticaBold34Font();
MgFont *mgTimes8Font();
MgFont *mgTimes10Font();
MgFont *mgTimes12Font();
MgFont *mgTimes14Font();
MgFont *mgTimes18Font();
MgFont *mgTimes24Font();
MgFont *mgTimes34Font();

/* free Meslo font */
MgFont *mgMenloMediumFont();

void mgText(struct memGfx *mg, int x, int y, Color color, 
	MgFont *font, char *text);
/* Draw a line of text with upper left corner x,y. */

void mgTextCentered(struct memGfx *mg, int x, int y, int width, int height, 
	Color color, MgFont *font, char *text);
/* Draw a line of text centered in box defined by x/y/width/height */

void mgTextRight(struct memGfx *mg, int x, int y, int width, int height, 
	Color color, MgFont *font, char *text);
/* Draw a line of text right justified in box defined by x/y/width/height */

int mgFontPixelHeight(MgFont *font);
/* How high in pixels is font? */

int mgFontLineHeight(MgFont *font);
/* How many pixels to next line ideally? */

int mgFontWidth(MgFont *font, char *chars, int charCount);
/* How wide are a couple of letters? */

int mgFontStringWidth(MgFont *font, char *string);
/* How wide is a string? */

int mgFontCharWidth(MgFont *font, char c);
/* How wide is a character? */

char *mgFontSizeBackwardsCompatible(char *size);
/* Given "size" argument that may be in old tiny/small/medium/big/huge format,
 * return it in new numerical string format. Do NOT free the return string*/

MgFont *mgFontForSizeAndStyle(char *textSize, char *fontType);
/* Get a font of given size and style.  Abort with error message if not found.
 * The textSize should be 6,8,10,12,14,18,24 or 34.  For backwards compatibility
 * textSizes of "tiny" "small", "medium", "large" and "huge" are also ok.
 * The fontType should be "medium", "bold", or "fixed" */

MgFont *mgFontForSize(char *textSize);
/* Get a font of given size and style.  Abort with error message if not found.
 * The textSize should be 6,8,10,12,14,18,24 or 34.  For backwards compatibility
 * textSizes of "tiny" "small", "medium", "large" and "huge" are also ok. */

void mgFillUnder(struct memGfx *mg, int x1, int y1, int x2, int y2, 
	int bottom, Color color);
/* Draw a 4 sided filled figure that has line x1/y1 to x2/y2 at
 * it's top, a horizontal line at bottom at it's bottom,  and
 * vertical lines from the bottom to y1 on the left and bottom to
 * y2 on the right. */

struct memGfx *mgRotate90(struct memGfx *in);
/* Create a copy of input that is rotated 90 degrees clockwise. */

void mgCircle(struct memGfx *mg, int xCen, int yCen, int rad, 
	Color color, boolean filled);
/* Draw a circle using a stepping algorithm.  Doesn't correct
 * for non-square pixels. */

void mgDrawPoly(struct memGfx *mg, struct gfxPoly *poly, Color color,
	boolean filled);
/* Draw polygon, possibly filled in color. */

void mgEllipse(struct memGfx *mg, int x0, int y0, int x1, int y1, Color color,
                        int mode, boolean isDashed);
/* Draw an ellipse (or limit to top or bottom) specified by rectangle, using Bresenham algorithm.
 * Optionally, alternate dots.
 * Point 0 is left, point 1 is top of rectangle
 * Adapted trivially from code posted at http://members.chello.at/~easyfilter/bresenham.html
 * Author: Zingl Alois, 8/22/2016
 */

/* Ellipse drawing modes */
#define ELLIPSE_FULL    0
#define ELLIPSE_TOP     1
#define ELLIPSE_BOTTOM  2

int mgCurve(struct memGfx *mg, int x0, int y0, int x1, int y1, int x2, int y2, Color color,
                        boolean isDashed);
/* Draw a segment of an anti-aliased curve within 3 points (quadratic Bezier)
 * Return max y value. Optionally draw curve as dashed line.
 * Adapted trivially from code posted at http://members.chello.at/~easyfilter/bresenham.html
 * Author: Zingl Alois, 8/22/2016
 */
/* TODO: allow specifying a third point on the line
 *  P(t) = (1-t)^2 * p0 + 2 * (1-t) * t * p1 + t^2 * p2
 */

struct hslColor mgRgbToHsl(struct rgbColor rgb);
/* Convert RGB to HSL colorspace (see http://en.wikipedia.org/wiki/HSL_and_HSV) 
 * In HSL, Hue is the color in the range [0,360) with 0=red 120=green 240=blue,
 * Saturation goes from a shade of grey (0) to fully saturated color (1000), and
 * Lightness goes from black (0) through the hue (500) to white (1000). */

struct hsvColor mgRgbToHsv(struct rgbColor rgb);
/* Convert RGB to HSV colorspace (see http://en.wikipedia.org/wiki/HSL_and_HSV)
 * In HSV, Hue is the color in the range [0,360) with 0=red 120=green 240=blue,
 * Saturation goes from white (0) to fully saturated color (1000), and
 * Value goes from black (0) through to the hue (1000). */
#define hsvValMax 1000
#define hsvSatMax 1000

struct rgbColor mgHslToRgb(struct hslColor hsl);
/* Convert HSL to RGB colorspace (see http://en.wikipedia.org/wiki/HSL_and_HSV) */

struct rgbColor mgHsvToRgb(struct hsvColor hsv);
/* Convert HSV to RGB colorspace (see http://en.wikipedia.org/wiki/HSL_and_HSV) */

struct rgbColor mgRgbTransformHsl(struct rgbColor in, double h, double s, double l);
/* Transform rgb 'in' value using
 *   hue shift 'h' (0..360 degrees), 
 *   saturation scale 's', and 
 *   lightness scale 'l'
 * Returns the transformed rgb value 
 * Use H=0, S=L=1 for identity transformation
 */

struct rgbColor mgRgbTransformHsv(struct rgbColor in, double h, double s, double v);
/* Transform rgb 'in' value using
 *   hue shift 'h' (0..360 degrees), 
 *   saturation scale 's', and 
 *   value scale 'v'
 * Returns the transformed rgb value 
 * Use H=0, S=V=1 for identity transformation
 */

struct rgbColor mgColorIxToRgb(struct memGfx *mg, int colorIx);
/* Return rgb value at color index. */

struct rgbColor colorIxToRgb(int colorIx);
/* Return rgb value at color index. */

#endif /* MEMGFX_H */
