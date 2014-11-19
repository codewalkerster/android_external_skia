/* NEON optimized code (C) COPYRIGHT 2009 Motorola
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 */

#include "SkBitmapProcState.h"
#include "SkPerspIter.h"
#include "SkShader.h"
#include "SkUtils.h"
#include "SkUtilsArm.h"
#include "SkColorPriv.h"

// Helper to ensure that when we shift down, we do it w/o sign-extension
// so the caller doesn't have to manually mask off the top 16 bits
//
static unsigned SK_USHIFT16(unsigned x) {
    return x >> 16;
}

/*  returns 0...(n-1) given any x (positive or negative).

    As an example, if n (which is always positive) is 5...

          x: -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8
    returns:  2  3  4  0  1  2  3  4  0  1  2  3  4  0  1  2  3
 */
static inline int sk_int_mod(int x, int n) {
    SkASSERT(n > 0);
    if ((unsigned)x >= (unsigned)n) {
        if (x < 0) {
            x = n + ~(~x % n);
        } else {
            x = x % n;
        }
    }
    return x;
}

/*
 *  The decal_ functions require that
 *  1. dx > 0
 *  2. [fx, fx+dx, fx+2dx, fx+3dx, ... fx+(count-1)dx] are all <= maxX
 *
 *  In addition, we use SkFractionalInt to keep more fractional precision than
 *  just SkFixed, so we will abort the decal_ call if dx is very small, since
 *  the decal_ function just operates on SkFixed. If that were changed, we could
 *  skip the very_small test here.
 */
static inline bool can_truncate_to_fixed_for_decal(SkFractionalInt frX,
                                                   SkFractionalInt frDx,
                                                   int count, unsigned max) {
    SkFixed dx = SkFractionalIntToFixed(frDx);

    // if decal_ kept SkFractionalInt precision, this would just be dx <= 0
    // I just made up the 1/256. Just don't want to perceive accumulated error
    // if we truncate frDx and lose its low bits.
    if (dx <= SK_Fixed1 / 256) {
        return false;
    }

    // We cast to unsigned so we don't have to check for negative values, which
    // will now appear as very large positive values, and thus fail our test!
    SkFixed fx = SkFractionalIntToFixed(frX);
    return (unsigned)SkFixedFloorToInt(fx) <= max &&
           (unsigned)SkFixedFloorToInt(fx + dx * (count - 1)) < max;
}

void decal_nofilter_scale(uint32_t dst[], SkFixed fx, SkFixed dx, int count);
void decal_filter_scale(uint32_t dst[], SkFixed fx, SkFixed dx, int count);

// Compile neon code paths if needed
#if !SK_ARM_NEON_IS_NONE

// These are defined in src/opts/SkBitmapProcState_matrixProcs_neon.cpp
extern const SkBitmapProcState::MatrixProc ClampX_ClampY_Procs_neon[];
extern const SkBitmapProcState::MatrixProc RepeatX_RepeatY_Procs_neon[];

#endif // !SK_ARM_NEON_IS_NONE

// Compile non-neon code path if needed
#define MAKENAME(suffix)        ClampX_ClampY ## suffix
#define TILEX_PROCF(fx, max)    SkClampMax((fx) >> 16, max)
#define TILEY_PROCF(fy, max)    SkClampMax((fy) >> 16, max)
#define TILEX_LOW_BITS(fx, max) (((fx) >> 12) & 0xF)
#define TILEY_LOW_BITS(fy, max) (((fy) >> 12) & 0xF)
#define CHECK_FOR_DECAL
#if	defined(__ARM_HAVE_NEON)
    #include "src/opts/SkBitmapProcState_matrix_clamp_neon_t.h"
#else
    #include "SkBitmapProcState_matrix.h"
#endif

#define MAKENAME(suffix)        RepeatX_RepeatY ## suffix
#define TILEX_PROCF(fx, max)    SK_USHIFT16(((fx) & 0xFFFF) * ((max) + 1))
#define TILEY_PROCF(fy, max)    SK_USHIFT16(((fy) & 0xFFFF) * ((max) + 1))
#define TILEX_LOW_BITS(fx, max) ((((fx) & 0xFFFF) * ((max) + 1) >> 12) & 0xF)
#define TILEY_LOW_BITS(fy, max) ((((fy) & 0xFFFF) * ((max) + 1) >> 12) & 0xF)
#if	defined(__ARM_HAVE_NEON)
    #include "src/opts/SkBitmapProcState_matrix_repeat_neon_t.h"
#else
    #include "SkBitmapProcState_matrix.h"
#endif

#define MAKENAME(suffix)        GeneralXY ## suffix
#define PREAMBLE(state)         SkBitmapProcState::FixedTileProc tileProcX = (state).fTileProcX; (void) tileProcX; \
                                SkBitmapProcState::FixedTileProc tileProcY = (state).fTileProcY; (void) tileProcY; \
                                SkBitmapProcState::FixedTileLowBitsProc tileLowBitsProcX = (state).fTileLowBitsProcX; (void) tileLowBitsProcX; \
                                SkBitmapProcState::FixedTileLowBitsProc tileLowBitsProcY = (state).fTileLowBitsProcY; (void) tileLowBitsProcY
#define PREAMBLE_PARAM_X        , SkBitmapProcState::FixedTileProc tileProcX, SkBitmapProcState::FixedTileLowBitsProc tileLowBitsProcX
#define PREAMBLE_PARAM_Y        , SkBitmapProcState::FixedTileProc tileProcY, SkBitmapProcState::FixedTileLowBitsProc tileLowBitsProcY
#define PREAMBLE_ARG_X          , tileProcX, tileLowBitsProcX
#define PREAMBLE_ARG_Y          , tileProcY, tileLowBitsProcY
#define TILEX_PROCF(fx, max)    SK_USHIFT16(tileProcX(fx) * ((max) + 1))
#define TILEY_PROCF(fy, max)    SK_USHIFT16(tileProcY(fy) * ((max) + 1))
#define TILEX_LOW_BITS(fx, max) tileLowBitsProcX(fx, (max) + 1)
#define TILEY_LOW_BITS(fy, max) tileLowBitsProcY(fy, (max) + 1)
#include "SkBitmapProcState_matrix.h"

static inline U16CPU fixed_clamp(SkFixed x)
{
    if (x < 0) {
        x = 0;
    }
    if (x >> 16) {
        x = 0xFFFF;
    }
    return x;
}

static inline U16CPU fixed_repeat(SkFixed x)
{
    return x & 0xFFFF;
}

// Visual Studio 2010 (MSC_VER=1600) optimizes bit-shift code incorrectly.
// See http://code.google.com/p/skia/issues/detail?id=472
#if defined(_MSC_VER) && (_MSC_VER >= 1600)
#pragma optimize("", off)
#endif

static inline U16CPU fixed_mirror(SkFixed x)
{
    SkFixed s = x << 15 >> 31;
    // s is FFFFFFFF if we're on an odd interval, or 0 if an even interval
    return (x ^ s) & 0xFFFF;
}

#if defined(_MSC_VER) && (_MSC_VER >= 1600)
#pragma optimize("", on)
#endif

static SkBitmapProcState::FixedTileProc choose_tile_proc(unsigned m)
{
    if (SkShader::kClamp_TileMode == m)
        return fixed_clamp;
    if (SkShader::kRepeat_TileMode == m)
        return fixed_repeat;
    SkASSERT(SkShader::kMirror_TileMode == m);
    return fixed_mirror;
}

static inline U16CPU fixed_clamp_lowbits(SkFixed x, int) {
    return (x >> 12) & 0xF;
}

static inline U16CPU fixed_repeat_or_mirrow_lowbits(SkFixed x, int scale) {
    return ((x * scale) >> 12) & 0xF;
}

static SkBitmapProcState::FixedTileLowBitsProc choose_tile_lowbits_proc(unsigned m) {
    if (SkShader::kClamp_TileMode == m) {
        return fixed_clamp_lowbits;
    } else {
        SkASSERT(SkShader::kMirror_TileMode == m ||
                 SkShader::kRepeat_TileMode == m);
        // mirror and repeat have the same behavior for the low bits.
        return fixed_repeat_or_mirrow_lowbits;
    }
}

static inline U16CPU int_clamp(int x, int n) {
    if (x >= n) {
        x = n - 1;
    }
    if (x < 0) {
        x = 0;
    }
    return x;
}

static inline U16CPU int_repeat(int x, int n) {
    return sk_int_mod(x, n);
}

static inline U16CPU int_mirror(int x, int n) {
    x = sk_int_mod(x, 2 * n);
    if (x >= n) {
        x = n + ~(x - n);
    }
    return x;
}

#if 0
static void test_int_tileprocs() {
    for (int i = -8; i <= 8; i++) {
        SkDebugf(" int_mirror(%2d, 3) = %d\n", i, int_mirror(i, 3));
    }
}
#endif

static SkBitmapProcState::IntTileProc choose_int_tile_proc(unsigned tm) {
    if (SkShader::kClamp_TileMode == tm)
        return int_clamp;
    if (SkShader::kRepeat_TileMode == tm)
        return int_repeat;
    SkASSERT(SkShader::kMirror_TileMode == tm);
    return int_mirror;
}

//////////////////////////////////////////////////////////////////////////////

void decal_nofilter_scale(uint32_t dst[], SkFixed fx, SkFixed dx, int count)
{
    int i;

    for (i = (count >> 2); i > 0; --i)
    {
        *dst++ = pack_two_shorts(fx >> 16, (fx + dx) >> 16);
        fx += dx+dx;
        *dst++ = pack_two_shorts(fx >> 16, (fx + dx) >> 16);
        fx += dx+dx;
    }
    count &= 3;

    uint16_t* xx = (uint16_t*)dst;
    for (i = count; i > 0; --i) {
        *xx++ = SkToU16(fx >> 16); fx += dx;
    }
}

void decal_filter_scale(uint32_t dst[], SkFixed fx, SkFixed dx, int count)
{


    if (count & 1)
    {
        SkASSERT((fx >> (16 + 14)) == 0);
        *dst++ = (fx >> 12 << 14) | ((fx >> 16) + 1);
        fx += dx;
    }
    while ((count -= 2) >= 0)
    {
        SkASSERT((fx >> (16 + 14)) == 0);
        *dst++ = (fx >> 12 << 14) | ((fx >> 16) + 1);
        fx += dx;

        *dst++ = (fx >> 12 << 14) | ((fx >> 16) + 1);
        fx += dx;
    }
}

///////////////////////////////////////////////////////////////////////////////
// stores the same as SCALE, but is cheaper to compute. Also since there is no
// scale, we don't need/have a FILTER version

static void fill_sequential(uint16_t xptr[], int start, int count) {
#if 1
    if (reinterpret_cast<intptr_t>(xptr) & 0x2) {
        *xptr++ = start++;
        count -= 1;
    }
    if (count > 3) {
        uint32_t* xxptr = reinterpret_cast<uint32_t*>(xptr);
        uint32_t pattern0 = PACK_TWO_SHORTS(start + 0, start + 1);
        uint32_t pattern1 = PACK_TWO_SHORTS(start + 2, start + 3);
        start += count & ~3;
        int qcount = count >> 2;
        do {
            *xxptr++ = pattern0;
            pattern0 += 0x40004;
            *xxptr++ = pattern1;
            pattern1 += 0x40004;
        } while (--qcount != 0);
        xptr = reinterpret_cast<uint16_t*>(xxptr);
        count &= 3;
    }
    while (--count >= 0) {
        *xptr++ = start++;
    }
#else
    for (int i = 0; i < count; i++) {
        *xptr++ = start++;
    }
#endif
}

static int nofilter_trans_preamble(const SkBitmapProcState& s, uint32_t** xy,
                                   int x, int y) {
    SkPoint pt;
    s.fInvProc(s.fInvMatrix, SkIntToScalar(x) + SK_ScalarHalf,
               SkIntToScalar(y) + SK_ScalarHalf, &pt);
    **xy = s.fIntTileProcY(SkScalarToFixed(pt.fY) >> 16,
                           s.fBitmap->height());
    *xy += 1;   // bump the ptr
    // return our starting X position
    return SkScalarToFixed(pt.fX) >> 16;
}

void clampx_nofilter_trans(const SkBitmapProcState& s,
                                  uint32_t xy[], int count, int x, int y) {
    SkASSERT((s.fInvType & ~SkMatrix::kTranslate_Mask) == 0);

    int xpos = nofilter_trans_preamble(s, &xy, x, y);
    const int width = s.fBitmap->width();
    if (1 == width) {
        // all of the following X values must be 0
        memset(xy, 0, count * sizeof(uint16_t));
        return;
    }

    uint16_t* xptr = reinterpret_cast<uint16_t*>(xy);
    int n;

    // fill before 0 as needed
    if (xpos < 0) {
        n = -xpos;
        if (n > count) {
            n = count;
        }
        memset(xptr, 0, n * sizeof(uint16_t));
        count -= n;
        if (0 == count) {
            return;
        }
        xptr += n;
        xpos = 0;
    }

    // fill in 0..width-1 if needed
    if (xpos < width) {
        n = width - xpos;
        if (n > count) {
            n = count;
        }
        fill_sequential(xptr, xpos, n);
        count -= n;
        if (0 == count) {
            return;
        }
        xptr += n;
    }

    // fill the remaining with the max value
    sk_memset16(xptr, width - 1, count);
}

/*
 * tao.zeng@amlogic.com, combine operation of clampx_nofilter_trans
 * and S16_opaque_D32_nofilter_DX
 */
extern void S16_D32_decal_nofilter_scale_t(uint16_t *srcAddr, int fx, int dx, SkPMColor dstC[], int count);
void ClampX_S16_D32_nofilter_trans_t(const     SkBitmapProcState& s,
                                     int       x,
                                     int       y,
                                     SkPMColor dstC[],
                                     int       count)
{
    if (count == 0) {
        return ;    
    }
    
    uint16_t* SK_RESTRICT srcAddr = (uint16_t *)s.fBitmap->getPixels();
    int xpos;

    SkPoint pt;
    s.fInvProc(s.fInvMatrix, SkIntToScalar(x) + SK_ScalarHalf,
               SkIntToScalar(y) + SK_ScalarHalf, &pt);
    int tmp =s.fIntTileProcY(SkScalarToFixed(pt.fY) >> 16,
                             s.fBitmap->height());
    srcAddr = (uint16_t *)((char*)srcAddr + tmp * s.fBitmap->rowBytes());

    xpos = SkScalarToFixed(pt.fX) >> 16;

    const int width = s.fBitmap->width();    
    if (1 == width) {
        // all of the following X values must be 0
        sk_memset32(dstC, SkPixel16ToPixel32(srcAddr[0]), count);
        return;
    }
    
    int n;

    // fill before 0 as needed
    if (xpos < 0) {
        n = -xpos;
        if (n > count) {
            n = count;
        }
        sk_memset32(dstC, SkPixel16ToPixel32(srcAddr[0]), n);
        count -= n;
        if (0 == count) {
            return;
        }
        dstC += n;
        xpos = 0;
    }

    // fill in 0..width-1 if needed
    if (xpos < width) {
        n = width - xpos;
        if (n > count) {
            n = count;
        }
        S16_D32_decal_nofilter_scale_t(srcAddr, xpos << 16, 0x10000, dstC, n);
        count -= n;
        if (0 == count) {
            return;
        }
        dstC += n;
    }

    // fill the remaining with the max value
    sk_memset32(dstC, SkPixel16ToPixel32(srcAddr[width - 1]), count);
}

// S32_opaque_D32_nofilter_DX + clampx_nofilter_trans
void ClampX_S32_D32_nofilter_trans_t(const     SkBitmapProcState& s,
                                     int       x,
                                     int       y,
                                     SkPMColor dstC[],
                                     int       count)
{
    if (count == 0) {
        return ;    
    }
    
    uint32_t* SK_RESTRICT srcAddr = (uint32_t *)s.fBitmap->getPixels();
    int xpos;

    SkPoint pt;
    s.fInvProc(s.fInvMatrix, SkIntToScalar(x) + SK_ScalarHalf,
               SkIntToScalar(y) + SK_ScalarHalf, &pt);
    uint32_t tmp =s.fIntTileProcY(SkScalarToFixed(pt.fY) >> 16,
                                  s.fBitmap->height());
    srcAddr = (uint32_t *)((char*)srcAddr + tmp * s.fBitmap->rowBytes());

    xpos = SkScalarToFixed(pt.fX) >> 16;

    const int width = s.fBitmap->width();    
    if (1 == width) {
        // all of the following X values must be 0
        sk_memset32(dstC, srcAddr[0], count);
        return;
    }
    
    int n;

    // fill before 0 as needed
    if (xpos < 0) {
        n = -xpos;
        if (n > count) {
            n = count;
        }
        sk_memset32(dstC, srcAddr[0], n);
        count -= n;
        if (0 == count) {
            return;
        }
        dstC += n;
        xpos = 0;
    }

    // fill in 0..width-1 if needed
    if (xpos < width) {
        n = width - xpos;
        if (n > count) {
            n = count;
        }
        memcpy(dstC, srcAddr + xpos, n * sizeof(uint32_t));
        count -= n;
        if (0 == count) {
            return;
        }
        dstC += n;
    }

    // fill the remaining with the max value
    sk_memset32(dstC, srcAddr[width - 1], count);
}

static void repeatx_nofilter_trans(const SkBitmapProcState& s,
                                   uint32_t xy[], int count, int x, int y) {
    SkASSERT((s.fInvType & ~SkMatrix::kTranslate_Mask) == 0);

    int xpos = nofilter_trans_preamble(s, &xy, x, y);
    const int width = s.fBitmap->width();
    if (1 == width) {
        // all of the following X values must be 0
        memset(xy, 0, count * sizeof(uint16_t));
        return;
    }

    uint16_t* xptr = reinterpret_cast<uint16_t*>(xy);
    int start = sk_int_mod(xpos, width);
    int n = width - start;
    if (n > count) {
        n = count;
    }
    fill_sequential(xptr, start, n);
    xptr += n;
    count -= n;

    while (count >= width) {
        fill_sequential(xptr, 0, width);
        xptr += width;
        count -= width;
    }

    if (count > 0) {
        fill_sequential(xptr, 0, count);
    }
}

// tao.zeng, add
// S32_opaque_D32_nofilter_DX + repeatx_nofilter_trans
void Repeatx_S32_D32_nofilter_trans_t(const     SkBitmapProcState& s,
                                      int       x,
                                      int       y,
                                      SkPMColor dstC[],
                                      int       count)
{
    SkASSERT((s.fInvType & ~SkMatrix::kTranslate_Mask) == 0);

    int xpos;  
    
    const uint32_t * SK_RESTRICT srcAddr = (const uint32_t*)s.fBitmap->getPixels();
    SkPoint pt;
    s.fInvProc(s.fInvMatrix, SkIntToScalar(x) + SK_ScalarHalf,
               SkIntToScalar(y) + SK_ScalarHalf, &pt);
    uint32_t tmp = s.fIntTileProcY(SkScalarToFixed(pt.fY) >> 16,
                                   s.fBitmap->height());
    xpos = SkScalarToFixed(pt.fX) >> 16;

    srcAddr = (uint32_t *)((char*)srcAddr + tmp * s.fBitmap->rowBytes());
    const int width = s.fBitmap->width();    
    if (1 == width) {
        // all of the following X values must be 0
        sk_memset32(dstC, srcAddr[0], count);
        return;
    }

    int start = sk_int_mod(xpos, width);
    int n = width - start;
    if (n > count) {
        n = count;
    }
    memcpy(dstC, srcAddr + start, n * sizeof(uint32_t));
    dstC += n;
    count -= n;

    while (count >= width) {
        memcpy(dstC, srcAddr, width * sizeof(uint32_t));
        dstC += width;
        count -= width;
    }

    if (count > 0) {
        memcpy(dstC, srcAddr, count * sizeof(uint32_t));
    }
}

static void fill_backwards(uint16_t xptr[], int pos, int count) {
    for (int i = 0; i < count; i++) {
        SkASSERT(pos >= 0);
        xptr[i] = pos--;
    }
}

static void mirrorx_nofilter_trans(const SkBitmapProcState& s,
                                   uint32_t xy[], int count, int x, int y) {
    SkASSERT((s.fInvType & ~SkMatrix::kTranslate_Mask) == 0);

    int xpos = nofilter_trans_preamble(s, &xy, x, y);
    const int width = s.fBitmap->width();
    if (1 == width) {
        // all of the following X values must be 0
        memset(xy, 0, count * sizeof(uint16_t));
        return;
    }

    uint16_t* xptr = reinterpret_cast<uint16_t*>(xy);
    // need to know our start, and our initial phase (forward or backward)
    bool forward;
    int n;
    int start = sk_int_mod(xpos, 2 * width);
    if (start >= width) {
        start = width + ~(start - width);
        forward = false;
        n = start + 1;  // [start .. 0]
    } else {
        forward = true;
        n = width - start;  // [start .. width)
    }
    if (n > count) {
        n = count;
    }
    if (forward) {
        fill_sequential(xptr, start, n);
    } else {
        fill_backwards(xptr, start, n);
    }
    forward = !forward;
    xptr += n;
    count -= n;

    while (count >= width) {
        if (forward) {
            fill_sequential(xptr, 0, width);
        } else {
            fill_backwards(xptr, width - 1, width);
        }
        forward = !forward;
        xptr += width;
        count -= width;
    }

    if (count > 0) {
        if (forward) {
            fill_sequential(xptr, 0, count);
        } else {
            fill_backwards(xptr, width - 1, count);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

#include <cutils/log.h>                                                 // tao.zeng, for debug

#define OPEN_LOG                                1
#define LOGD_IF     ALOGD_IF

#define CLAMPX_NOFILTER_TRANS                   0
#define REPEATX_NOFILTER_TRANS                  1
#define MIRRORX_NOFILTER_TRANS                  2
#define CLAMPX_CLAMPY_NOFILTER_SCALE_NEON       3
#define CLAMPX_CLAMPY_FILTER_SCALE_NEON         4
#define CLAMPX_CLAMPY_NOFILTER_AFFINE_NEON      5
#define CLAMPX_CLAMPY_FILTER_AFFINE_NEON        6
#define CLAMPX_CLAMPY_NOFILTER_PERSP_NEON       7
#define CLAMPX_CLAMPY_FILTER_PERSP_NEON         8
#define REPEATX_REPEATY_NOFILTER_SCALE_NEON     9 
#define REPEATX_REPEATY_FILTER_SCALE_NEON       10
#define REPEATX_REPEATY_NOFILTER_AFFINE_NEON    11
#define REPEATX_REPEATY_FILTER_AFFINE_NEON      12
#define REPEATX_REPEATY_NOFILTER_PERSP_NEON     13
#define REPEATX_REPEATY_FILTER_PERSP_NEON       14
#define GENERALXY_NOFILTER_SCALE_NEON           15 
#define GENERALXY_FILTER_SCALE_NEON             16
#define GENERALXY_NOFILTER_AFFINE_NEON          17
#define GENERALXY_FILTER_AFFINE_NEON            18
#define GENERALXY_NOFILTER_PERSP_NEON           19
#define GENERALXY_FILTER_PERSP_NEON             20 

int getMatrixProcName(void *pfunc)
{
    int i = -1;

    if (pfunc == clampx_nofilter_trans) {
        return CLAMPX_NOFILTER_TRANS;
    }
    if (pfunc == repeatx_nofilter_trans) {
        return REPEATX_NOFILTER_TRANS;
    }
    if (pfunc == mirrorx_nofilter_trans) {
        return MIRRORX_NOFILTER_TRANS;
    }
    // ClampX_ClampY ## xxx
    if (pfunc == ClampX_ClampY_nofilter_scale_neon) {
        return CLAMPX_CLAMPY_NOFILTER_SCALE_NEON;
    }
    if (pfunc == ClampX_ClampY_filter_scale_neon) {
        return CLAMPX_CLAMPY_FILTER_SCALE_NEON;
    }
    if (pfunc == ClampX_ClampY_nofilter_affine_neon) {
        return CLAMPX_CLAMPY_NOFILTER_AFFINE_NEON;
    }
    if (pfunc == ClampX_ClampY_filter_affine_neon) {
        return CLAMPX_CLAMPY_FILTER_AFFINE_NEON;
    }
    if (pfunc == ClampX_ClampY_nofilter_persp_neon) {
        return CLAMPX_CLAMPY_NOFILTER_PERSP_NEON;
    }
    if (pfunc == ClampX_ClampY_filter_persp_neon) {
        return CLAMPX_CLAMPY_FILTER_PERSP_NEON;
    }
    // RepeatX_RepeatY ## xxx
    if (pfunc == RepeatX_RepeatY_nofilter_scale_neon) {
        return REPEATX_REPEATY_NOFILTER_SCALE_NEON;
    }
    if (pfunc == RepeatX_RepeatY_filter_scale) {
        return REPEATX_REPEATY_FILTER_SCALE_NEON;
    }
    if (pfunc == RepeatX_RepeatY_nofilter_affine_neon) {
        return REPEATX_REPEATY_NOFILTER_AFFINE_NEON;
    }
    if (pfunc == RepeatX_RepeatY_filter_affine) {
        return REPEATX_REPEATY_FILTER_AFFINE_NEON;
    }
    if (pfunc == RepeatX_RepeatY_nofilter_persp_neon) {
        return REPEATX_REPEATY_NOFILTER_PERSP_NEON;
    }
    if (pfunc == RepeatX_RepeatY_filter_persp) {
        return REPEATX_REPEATY_FILTER_PERSP_NEON;
    }
    // GeneralXY ## xxx
    if (pfunc == GeneralXY_nofilter_scale) {
        return GENERALXY_NOFILTER_SCALE_NEON;
    }
    if (pfunc == GeneralXY_filter_scale) {
        return GENERALXY_FILTER_SCALE_NEON;
    }
    if (pfunc == GeneralXY_nofilter_affine) {
        return GENERALXY_NOFILTER_AFFINE_NEON;
    }
    if (pfunc == GeneralXY_filter_affine) {
        return GENERALXY_FILTER_AFFINE_NEON;
    }
    if (pfunc == GeneralXY_nofilter_persp) {
        return GENERALXY_NOFILTER_PERSP_NEON;
    }
    if (pfunc == GeneralXY_filter_persp) {
        return GENERALXY_FILTER_PERSP_NEON;
    }
    return i;
}

void disp_matrix_proc_name(void *pfunc)
{

    if (pfunc == clampx_nofilter_trans) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "clampx_nofilter_trans");
        return ;
    } else if (pfunc == repeatx_nofilter_trans) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "repeatx_nofilter_trans");
        return ;
    } else if (pfunc == mirrorx_nofilter_trans) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "mirrorx_nofilter_trans");
        return ;
    } else
    // ClampX_ClampY ## xxx
    if (pfunc == ClampX_ClampY_nofilter_scale_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "ClampX_ClampY_nofilter_scale_neon");
        return ;
    } else if (pfunc == ClampX_ClampY_filter_scale_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "ClampX_ClampY_filter_scale_neon");
        return ;
    } else if (pfunc == ClampX_ClampY_nofilter_affine_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "ClampX_ClampY_nofilter_affine_neon");
        return ;
    } else if (pfunc == ClampX_ClampY_filter_affine_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "ClampX_ClampY_filter_affine_neon");
        return ;
    } else if (pfunc == ClampX_ClampY_nofilter_persp_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "ClampX_ClampY_nofilter_persp_neon");
        return ;
    } else if (pfunc == ClampX_ClampY_filter_persp_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "ClampX_ClampY_filter_persp_neon");
        return ;
    } else
    // RepeatX_RepeatY ## xxx
    if (pfunc == RepeatX_RepeatY_nofilter_scale_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "RepeatX_RepeatY_nofilter_scale_neon");
        return ;
    } else if (pfunc == RepeatX_RepeatY_filter_scale) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "RepeatX_RepeatY_filter_scale");
        return ;
    } else if (pfunc == RepeatX_RepeatY_nofilter_affine_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "RepeatX_RepeatY_nofilter_affine_neon");
        return ;
    } else if (pfunc == RepeatX_RepeatY_filter_affine) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "RepeatX_RepeatY_filter_affine");
        return ;
    } else if (pfunc == RepeatX_RepeatY_nofilter_persp_neon) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "RepeatX_RepeatY_nofilter_persp_neon");
        return ; 
    } else if (pfunc == RepeatX_RepeatY_filter_persp) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "RepeatX_RepeatY_filter_persp");
        return ;
    } else
    // GeneralXY ## xxx
    if (pfunc == GeneralXY_nofilter_scale) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "GeneralXY_nofilter_scale");
        return ;
    } else if (pfunc == GeneralXY_filter_scale) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "GeneralXY_filter_scale");
        return ;
    } else if (pfunc == GeneralXY_nofilter_affine) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "GeneralXY_nofilter_affine");
        return ;
    } else if (pfunc == GeneralXY_filter_affine) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "GeneralXY_filter_affine");
        return ;
    } else if (pfunc == GeneralXY_nofilter_persp) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "GeneralXY_nofilter_persp");
        return ;
    } else if (pfunc == GeneralXY_filter_persp) {
        LOGD_IF(OPEN_LOG, "MatrixProc:%s", "GeneralXY_filter_persp");
        return ;
    }
}

SkBitmapProcState::MatrixProc
SkBitmapProcState::chooseMatrixProc(bool trivial_matrix) {
//    test_int_tileprocs();
    // check for our special case when there is no scale/affine/perspective
    if (trivial_matrix) {
        SkASSERT(SkPaint::kNone_FilterLevel == fFilterLevel);
        fIntTileProcY = choose_int_tile_proc(fTileModeY);
        switch (fTileModeX) {
            case SkShader::kClamp_TileMode:
                return clampx_nofilter_trans;
            case SkShader::kRepeat_TileMode:
                return repeatx_nofilter_trans;
            case SkShader::kMirror_TileMode:
                return mirrorx_nofilter_trans;
        }
    }

    int index = 0;
    if (fFilterLevel != SkPaint::kNone_FilterLevel) {
        index = 1;
    }
    if (fInvType & SkMatrix::kPerspective_Mask) {
        index += 4;
    } else if (fInvType & SkMatrix::kAffine_Mask) {
        index += 2;
    }

    if (SkShader::kClamp_TileMode == fTileModeX &&
        SkShader::kClamp_TileMode == fTileModeY)
    {
        // clamp gets special version of filterOne
        fFilterOneX = SK_Fixed1;
        fFilterOneY = SK_Fixed1;
        return SK_ARM_NEON_WRAP(ClampX_ClampY_Procs)[index];
    }

    // all remaining procs use this form for filterOne
    fFilterOneX = SK_Fixed1 / fBitmap->width();
    fFilterOneY = SK_Fixed1 / fBitmap->height();

    if (SkShader::kRepeat_TileMode == fTileModeX &&
        SkShader::kRepeat_TileMode == fTileModeY)
    {
        return SK_ARM_NEON_WRAP(RepeatX_RepeatY_Procs)[index];
    }

    fTileProcX = choose_tile_proc(fTileModeX);
    fTileProcY = choose_tile_proc(fTileModeY);
    fTileLowBitsProcX = choose_tile_lowbits_proc(fTileModeX);
    fTileLowBitsProcY = choose_tile_lowbits_proc(fTileModeY);
    return GeneralXY_Procs[index];
}
