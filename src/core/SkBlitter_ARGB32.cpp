/*
 * Copyright 2006 The Android Open Source Project
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 */

#include "SkCoreBlitters.h"
#include "SkColorPriv.h"
#include "SkShader.h"
#include "SkUtils.h"
#include "SkXfermode.h"
#include "SkBlitMask.h"
#define LOG_TAG "SKIA"
#include <utils/Log.h>
#if defined(FIMG2D_ENABLED)
#include "SkBitmap.h"
#include "SkBitmapProcShader.h"
#if defined(FIMG2D4X)
#include "SkFimgApi4x.h"
extern Fimg fimg;
extern SkMutex gG2DMutex;
#endif
#endif

///////////////////////////////////////////////////////////////////////////////
static inline int upscale31To32(int value) {
    SkASSERT((unsigned)value <= 31);
    return value + (value >> 4);
}

static inline int blend32(int src, int dst, int scale) {
    SkASSERT((unsigned)src <= 0xFF);
    SkASSERT((unsigned)dst <= 0xFF);
    SkASSERT((unsigned)scale <= 32);
    return dst + ((src - dst) * scale >> 5);
}

static void blit_lcd16_opaque(SkPMColor dst[], const uint16_t src[],
                              SkPMColor color, int width) {
    int srcR = SkGetPackedR32(color);
    int srcG = SkGetPackedG32(color);
    int srcB = SkGetPackedB32(color);

    for (int i = 0; i < width; i++) {
        uint16_t mask = src[i];
        if (0 == mask) {
            continue;
        }

        SkPMColor d = dst[i];
        
        /*  We want all of these in 5bits, hence the shifts in case one of them
         *  (green) is 6bits.
         */
        int maskR = SkGetPackedR16(mask) >> (SK_R16_BITS - 5);
        int maskG = SkGetPackedG16(mask) >> (SK_G16_BITS - 5);
        int maskB = SkGetPackedB16(mask) >> (SK_B16_BITS - 5);

        // Now upscale them to 0..256, so we can use SkAlphaBlend
        maskR = upscale31To32(maskR);
        maskG = upscale31To32(maskG);
        maskB = upscale31To32(maskB);

        int maskA = SkMax32(SkMax32(maskR, maskG), maskB);

        int dstA = SkGetPackedA32(d);
        int dstR = SkGetPackedR32(d);
        int dstG = SkGetPackedG32(d);
        int dstB = SkGetPackedB32(d);

        dst[i] = SkPackARGB32(blend32(0xFF, dstA, maskA),
                              blend32(srcR, dstR, maskR),
                              blend32(srcG, dstG, maskG),
                              blend32(srcB, dstB, maskB));
    }
}

static void blitmask_lcd16(const SkBitmap& device, const SkMask& mask,
                           const SkIRect& clip, SkPMColor srcColor) {
    int x = clip.fLeft;
    int y = clip.fTop;
    int width = clip.width();
    int height = clip.height();

    SkPMColor*		dstRow = device.getAddr32(x, y);
    const uint16_t* srcRow = mask.getAddrLCD16(x, y);

    do {
        blit_lcd16_opaque(dstRow, srcRow, srcColor, width);
        dstRow = (SkPMColor*)((char*)dstRow + device.rowBytes());
        srcRow = (const uint16_t*)((const char*)srcRow + mask.fRowBytes);
    } while (--height != 0);
}


static void SkARGB32_Blit32(const SkBitmap& device, const SkMask& mask,
							const SkIRect& clip, SkPMColor srcColor) {
	U8CPU alpha = SkGetPackedA32(srcColor);
	unsigned flags = SkBlitRow::kSrcPixelAlpha_Flag32;
	if (alpha != 255) {
		flags |= SkBlitRow::kGlobalAlpha_Flag32;
	}
	SkBlitRow::Proc32 proc = SkBlitRow::Factory32(flags);

    int x = clip.fLeft;
    int y = clip.fTop;
    int width = clip.width();
    int height = clip.height();

    SkPMColor*		 dstRow = device.getAddr32(x, y);
    const SkPMColor* srcRow = reinterpret_cast<const SkPMColor*>(mask.getAddr8(x, y));

    do {
		proc(dstRow, srcRow, width, alpha);
        dstRow = (SkPMColor*)((char*)dstRow + device.rowBytes());
        srcRow = (const SkPMColor*)((const char*)srcRow + mask.fRowBytes);
    } while (--height != 0);
}

//////////////////////////////////////////////////////////////////////////////////////

SkARGB32_Blitter::SkARGB32_Blitter(const SkBitmap& device, const SkPaint& paint)
        : INHERITED(device) {
    SkColor color = paint.getColor();
    fColor = color;

    fSrcA = SkColorGetA(color);
    unsigned scale = SkAlpha255To256(fSrcA);
    fSrcR = SkAlphaMul(SkColorGetR(color), scale);
    fSrcG = SkAlphaMul(SkColorGetG(color), scale);
    fSrcB = SkAlphaMul(SkColorGetB(color), scale);

    fPMColor = SkPackARGB32(fSrcA, fSrcR, fSrcG, fSrcB);
    fColor32Proc = SkBlitRow::ColorProcFactory();
    fBlitMaskProc = SkBlitMask::ColorFactory(SkBitmap::kARGB_8888_Config, SkMask::kARGB32_Format, color);
}

const SkBitmap* SkARGB32_Blitter::justAnOpaqueColor(uint32_t* value) {
    if (255 == fSrcA) {
        *value = fPMColor;
        return &fDevice;
    }
    return NULL;
}

#if defined _WIN32 && _MSC_VER >= 1300  // disable warning : local variable used without having been initialized
#pragma warning ( push )
#pragma warning ( disable : 4701 )
#endif

void SkARGB32_Blitter::blitH(int x, int y, int width) {
    SkASSERT(x >= 0 && y >= 0 && x + width <= fDevice.width());

    uint32_t*   device = fDevice.getAddr32(x, y);
    fColor32Proc(device, device, width, fPMColor);
}

void SkARGB32_Blitter::blitAntiH(int x, int y, const SkAlpha antialias[],
                                 const int16_t runs[]) {
    if (fSrcA == 0) {
        return;
    }

    uint32_t    color = fPMColor;
    uint32_t*   device = fDevice.getAddr32(x, y);
    unsigned    opaqueMask = fSrcA; // if fSrcA is 0xFF, then we will catch the fast opaque case

    for (;;) {
        int count = runs[0];
        SkASSERT(count >= 0);
        if (count <= 0) {
            return;
        }
        unsigned aa = antialias[0];
        if (aa) {
            if ((opaqueMask & aa) == 255) {
                sk_memset32(device, color, count);
            } else {
                uint32_t sc = SkAlphaMulQ(color, SkAlpha255To256(aa));
                fColor32Proc(device, device, count, sc);
            }
        }
        runs += count;
        antialias += count;
        device += count;
    }
}

//////////////////////////////////////////////////////////////////////////////////////

#define solid_8_pixels(mask, dst, color)    \
    do {                                    \
        if (mask & 0x80) dst[0] = color;    \
        if (mask & 0x40) dst[1] = color;    \
        if (mask & 0x20) dst[2] = color;    \
        if (mask & 0x10) dst[3] = color;    \
        if (mask & 0x08) dst[4] = color;    \
        if (mask & 0x04) dst[5] = color;    \
        if (mask & 0x02) dst[6] = color;    \
        if (mask & 0x01) dst[7] = color;    \
    } while (0)

#define SK_BLITBWMASK_NAME                  SkARGB32_BlitBW
#define SK_BLITBWMASK_ARGS                  , SkPMColor color
#define SK_BLITBWMASK_BLIT8(mask, dst)      solid_8_pixels(mask, dst, color)
#define SK_BLITBWMASK_GETADDR               getAddr32
#define SK_BLITBWMASK_DEVTYPE               uint32_t
#include "SkBlitBWMaskTemplate.h"

#define blend_8_pixels(mask, dst, sc, dst_scale)                            \
    do {                                                                    \
        if (mask & 0x80) { dst[0] = sc + SkAlphaMulQ(dst[0], dst_scale); }  \
        if (mask & 0x40) { dst[1] = sc + SkAlphaMulQ(dst[1], dst_scale); }  \
        if (mask & 0x20) { dst[2] = sc + SkAlphaMulQ(dst[2], dst_scale); }  \
        if (mask & 0x10) { dst[3] = sc + SkAlphaMulQ(dst[3], dst_scale); }  \
        if (mask & 0x08) { dst[4] = sc + SkAlphaMulQ(dst[4], dst_scale); }  \
        if (mask & 0x04) { dst[5] = sc + SkAlphaMulQ(dst[5], dst_scale); }  \
        if (mask & 0x02) { dst[6] = sc + SkAlphaMulQ(dst[6], dst_scale); }  \
        if (mask & 0x01) { dst[7] = sc + SkAlphaMulQ(dst[7], dst_scale); }  \
    } while (0)

#define SK_BLITBWMASK_NAME                  SkARGB32_BlendBW
#define SK_BLITBWMASK_ARGS                  , uint32_t sc, unsigned dst_scale
#define SK_BLITBWMASK_BLIT8(mask, dst)      blend_8_pixels(mask, dst, sc, dst_scale)
#define SK_BLITBWMASK_GETADDR               getAddr32
#define SK_BLITBWMASK_DEVTYPE               uint32_t
#include "SkBlitBWMaskTemplate.h"

void SkARGB32_Blitter::blitMask(const SkMask& mask, const SkIRect& clip) {
    SkASSERT(mask.fBounds.contains(clip));
    SkASSERT(fSrcA != 0xFF);

    if (fSrcA == 0) {
        return;
    }

    if (SkBlitMask::BlitColor(fDevice, mask, clip, fColor)) {
        return;
    }

    if (mask.fFormat == SkMask::kBW_Format) {
#if defined(FIMG2D4X)
#if (FIMGAPI_G2D_BLOKING == TRUE)
        gG2DMutex.acquire();
#else
        if (gG2DMutex.tryacquire() != 0) {
            SkARGB32_BlendBW(fDevice, mask, clip, fPMColor, SkAlpha255To256(255 - fSrcA));
            return;
        }
#endif
        int retFimg = FimgARGB32_BW(fimg, fDevice, mask, clip, fPMColor, SkAlpha255To256(255 - fSrcA));
        if (retFimg != FIMGAPI_FINISHED)
            SkARGB32_BlendBW(fDevice, mask, clip, fPMColor, SkAlpha255To256(255 - fSrcA));

        gG2DMutex.release();
#else
        SkARGB32_BlendBW(fDevice, mask, clip, fPMColor, SkAlpha255To256(255 - fSrcA));
#endif
    } else if (SkMask::kARGB32_Format == mask.fFormat) {
		SkARGB32_Blit32(fDevice, mask, clip, fPMColor);
    } else if (SkMask::kLCD16_Format == mask.fFormat) {
        blitmask_lcd16(fDevice, mask, clip, fPMColor);
        return;
    }
}

void SkARGB32_Opaque_Blitter::blitMask(const SkMask& mask,
                                       const SkIRect& clip) {
    SkASSERT(mask.fBounds.contains(clip));

    if (SkBlitMask::BlitColor(fDevice, mask, clip, fColor)) {
        return;
    }
    
    if (mask.fFormat == SkMask::kBW_Format) {
#if defined(FIMG2D4X)
#if (FIMGAPI_G2D_BLOKING == TRUE)
        gG2DMutex.acquire();
#else
        if (gG2DMutex.tryacquire() != 0) {
            SkARGB32_BlitBW(fDevice, mask, clip, fPMColor);
            return;
        }
#endif
        int retFimg = FimgARGB32_BW(fimg, fDevice, mask, clip, fPMColor, SkAlpha255To256(255));
        if (retFimg != FIMGAPI_FINISHED)
            SkARGB32_BlitBW(fDevice, mask, clip, fPMColor);

        gG2DMutex.release();
#else
        SkARGB32_BlitBW(fDevice, mask, clip, fPMColor);
#endif
    } else if (SkMask::kARGB32_Format == mask.fFormat) {
		SkARGB32_Blit32(fDevice, mask, clip, fPMColor);
	}

    int x = clip.fLeft;
    int y = clip.fTop;
    int width = clip.width();
    int height = clip.height();


#if defined(FIMG2D4X)

#if (FIMGAPI_G2D_BLOKING == TRUE)
    gG2DMutex.acquire();
#else
    if (gG2DMutex.tryacquire() != 0) {
        fBlitMaskProc(fDevice.getAddr32(x, y), fDevice.rowBytes(),
                      SkBitmap::kARGB_8888_Config,
                      mask.getAddr(x, y), mask.fRowBytes, fColor, width, height);
        return;
    }
#endif

    fimg.srcAddr        = (unsigned char *)NULL;

    fimg.isSolidFill    = false;
    fimg.fillcolor      = toARGB32(SkPreMultiplyColor(fColor));

    fimg.srcColorFormat = fDevice.config();

    fimg.dstX           = x;
    fimg.dstY           = y;
    fimg.dstW           = width;
    fimg.dstH           = height;
    fimg.dstFWStride    = fDevice.rowBytes();
    fimg.dstFH          = fDevice.height();
    fimg.dstBPP         = fDevice.bytesPerPixel();
    fimg.dstColorFormat = fDevice.config();
    fimg.dstAddr        = (unsigned char *)fDevice.getAddr(0,0);

    fimg.clipT = y;
    fimg.clipB = y + height;
    fimg.clipL = x;
    fimg.clipR = x + width;

    fimg.mskX           = 0;
    fimg.mskY           = 0;
    fimg.mskW           = width;
    fimg.mskH           = height;
    fimg.mskFWStride    = mask.fRowBytes;
    fimg.mskFH          = height;
    fimg.mskBPP         = 1;
    fimg.mskColorFormat = SkBitmap::kA8_Config;
    fimg.mskAddr        = (unsigned char *)mask.getAddr(x, y);

    fimg.rotate         = 0;

    fimg.xfermode = SkXfermode::kSrcOver_Mode;
    fimg.isDither = false;
    fimg.colorFilter = NULL;
    fimg.matrixType = 0;
    fimg.matrixSx = 0;
    fimg.matrixSy = 0;

    fimg.alpha = 0xFF;

    fimg.canusehybrid = 0;
    int retFimg = FimgApiStretch(&fimg, __func__);
    if (retFimg != FIMGAPI_FINISHED)
        fBlitMaskProc(fDevice.getAddr32(x, y), fDevice.rowBytes(),
//                      SkBitmap::kARGB_8888_Config,
                      mask.getAddr(x, y), mask.fRowBytes, fColor, width, height);

    gG2DMutex.release();

#else
    fBlitMaskProc(fDevice.getAddr32(x, y), fDevice.rowBytes(),
//                  SkBitmap::kARGB_8888_Config,
                  mask.getAddr(x, y), mask.fRowBytes, fColor, width, height);
#endif
}

///////////////////////////////////////////////////////////////////////////////

void SkARGB32_Blitter::blitV(int x, int y, int height, SkAlpha alpha) {
    if (alpha == 0 || fSrcA == 0) {
        return;
    }

    uint32_t* device = fDevice.getAddr32(x, y);
    uint32_t  color = fPMColor;

    if (alpha != 255) {
        color = SkAlphaMulQ(color, SkAlpha255To256(alpha));
    }

    unsigned dst_scale = 255 - SkGetPackedA32(color);
    uint32_t rowBytes = fDevice.rowBytes();
    while (--height >= 0) {
        device[0] = color + SkAlphaMulQ(device[0], dst_scale);
        device = (uint32_t*)((char*)device + rowBytes);
    }
}

void SkARGB32_Blitter::blitRect(int x, int y, int width, int height) {
    SkASSERT(x >= 0 && y >= 0 && x + width <= fDevice.width() && y + height <= fDevice.height());

    if (fSrcA == 0) {
        return;
    }

    uint32_t*   device = fDevice.getAddr32(x, y);
    uint32_t    color = fPMColor;
    size_t      rowBytes = fDevice.rowBytes();

    while (--height >= 0) {
        fColor32Proc(device, device, width, color);
        device = (uint32_t*)((char*)device + rowBytes);
    }
}

#if defined _WIN32 && _MSC_VER >= 1300
#pragma warning ( pop )
#endif

///////////////////////////////////////////////////////////////////////
void SkARGB32_Black_Blitter::blitMask(const SkMask& mask, const SkIRect& clip) {
    SkASSERT(mask.fBounds.contains(clip));

    if (mask.fFormat == SkMask::kBW_Format) {
        SkPMColor black = (SkPMColor)(SK_A32_MASK << SK_A32_SHIFT);

#if defined(FIMG2D4X)
#if (FIMGAPI_G2D_BLOKING == TRUE)
        gG2DMutex.acquire();
#else
        if (gG2DMutex.tryacquire() != 0) {
            SkARGB32_BlitBW(fDevice, mask, clip, black);
            return;
        }
#endif
        int retFimg = FimgARGB32_BW(fimg, fDevice, mask, clip, fPMColor, SkAlpha255To256(255));
        if (retFimg != FIMGAPI_FINISHED)
            SkARGB32_BlitBW(fDevice, mask, clip, black);

        gG2DMutex.release();
#else
        SkARGB32_BlitBW(fDevice, mask, clip, black);
#endif
    } else if (SkMask::kARGB32_Format == mask.fFormat) {
		SkARGB32_Blit32(fDevice, mask, clip, fPMColor);
    } else if (SkMask::kLCD16_Format == mask.fFormat) {
        blitmask_lcd16(fDevice, mask, clip, fPMColor);
    } else {
#if defined(SK_SUPPORT_LCDTEXT)
        const bool      lcdMode = mask.fFormat == SkMask::kHorizontalLCD_Format;
        const bool      verticalLCDMode = mask.fFormat == SkMask::kVerticalLCD_Format;
#endif

        // In LCD mode the masks have either an extra couple of rows or columns on the edges.
        unsigned        width = clip.width();
        unsigned        height = clip.height();

        SkASSERT((int)height > 0);
        SkASSERT((int)width > 0);

#if defined(SK_SUPPORT_LCDTEXT)
        if (lcdMode || verticalLCDMode) {
            int widthAdjustment, heightAdjustment;
            const uint32_t* alpha32;
            uint32_t* device = adjustForSubpixelClip(mask, clip, fDevice, &widthAdjustment, &heightAdjustment, &alpha32);

            width += widthAdjustment;
            height += heightAdjustment;

            unsigned deviceRB = fDevice.rowBytes() - (width << 2);
            unsigned alphaExtraRowWords = mask.rowWordsLCD() - width;

            do {
                unsigned w = width;
                do {
                    const uint32_t alphaPixel = *alpha32++;
                    const uint32_t originalPixel = *device;
                    *device++ = BlendLCDPixelWithBlack(alphaPixel, originalPixel);
                } while (--w != 0);
                device = (uint32_t*)((char*)device + deviceRB);
                alpha32 += alphaExtraRowWords;
            } while (--height != 0);

            return;
        }
#endif

        uint32_t*      device = fDevice.getAddr32(clip.fLeft, clip.fTop);
        unsigned       maskRB = mask.fRowBytes - width;
        unsigned       deviceRB = fDevice.rowBytes() - (width << 2);
        const uint8_t* alpha = (uint8_t*) mask.getAddr(clip.fLeft, clip.fTop);
        do {
            unsigned w = width;
            do {
                unsigned aa = *alpha++;
                *device = (aa << SK_A32_SHIFT) + SkAlphaMulQ(*device, SkAlpha255To256(255 - aa));
                device += 1;
            } while (--w != 0);
            device = (uint32_t*)((char*)device + deviceRB);
            alpha += maskRB;
        } while (--height != 0);
    }
}

void SkARGB32_Black_Blitter::blitAntiH(int x, int y, const SkAlpha antialias[],
                                       const int16_t runs[]) {
    uint32_t*   device = fDevice.getAddr32(x, y);
    SkPMColor   black = (SkPMColor)(SK_A32_MASK << SK_A32_SHIFT);

    for (;;) {
        int count = runs[0];
        SkASSERT(count >= 0);
        if (count <= 0) {
            return;
        }
        unsigned aa = antialias[0];
        if (aa) {
            if (aa == 255) {
                sk_memset32(device, black, count);
            } else {
                SkPMColor src = aa << SK_A32_SHIFT;
                unsigned dst_scale = 256 - aa;
                int n = count;
                do {
                    --n;
                    device[n] = src + SkAlphaMulQ(device[n], dst_scale);
                } while (n > 0);
            }
        }
        runs += count;
        antialias += count;
        device += count;
    }
}

///////////////////////////////////////////////////////////////////////////////

SkARGB32_Shader_Blitter::SkARGB32_Shader_Blitter(const SkBitmap& device,
                            const SkPaint& paint) : INHERITED(device, paint) {
    fBuffer = (SkPMColor*)sk_malloc_throw(device.width() * (sizeof(SkPMColor)));

    fXfermode = paint.getXfermode();
    SkSafeRef(fXfermode);

    int flags = 0;
    if (!(fShader->getFlags() & SkShader::kOpaqueAlpha_Flag)) {
        flags |= SkBlitRow::kSrcPixelAlpha_Flag32;
    }
    // we call this on the output from the shader
    fProc32 = SkBlitRow::Factory32(flags);
    // we call this on the output from the shader + alpha from the aa buffer
    fProc32Blend = SkBlitRow::Factory32(flags | SkBlitRow::kGlobalAlpha_Flag32);
}

SkARGB32_Shader_Blitter::~SkARGB32_Shader_Blitter() {
    SkSafeUnref(fXfermode);
    sk_free(fBuffer);
}

void SkARGB32_Shader_Blitter::blitH(int x, int y, int width) {
    SkASSERT(x >= 0 && y >= 0 && x + width <= fDevice.width());

    uint32_t*   device = fDevice.getAddr32(x, y);

    if (fXfermode == NULL && (fShader->getFlags() & SkShader::kOpaqueAlpha_Flag)) {
        fShader->shadeSpan(x, y, device, width);
    } else {
        SkPMColor*  span = fBuffer;
        fShader->shadeSpan(x, y, span, width);
        if (fXfermode) {
            fXfermode->xfer32(device, span, width, NULL);
        } else {
            fProc32(device, span, width, 255);
        }
    }
}

void SkARGB32_Shader_Blitter::blitAntiH(int x, int y, const SkAlpha antialias[],
                                        const int16_t runs[]) {
    SkPMColor*  span = fBuffer;
    uint32_t*   device = fDevice.getAddr32(x, y);
    SkShader*   shader = fShader;

    if (fXfermode) {
        for (;;) {
            SkXfermode* xfer = fXfermode;

            int count = *runs;
            if (count <= 0)
                break;
            int aa = *antialias;
            if (aa) {
                shader->shadeSpan(x, y, span, count);
                if (aa == 255) {
                    xfer->xfer32(device, span, count, NULL);
                } else {
                    // count is almost always 1
                    for (int i = count - 1; i >= 0; --i) {
                        xfer->xfer32(&device[i], &span[i], 1, antialias);
                    }
                }
            }
            device += count;
            runs += count;
            antialias += count;
            x += count;
        }
    } else if (fShader->getFlags() & SkShader::kOpaqueAlpha_Flag) {
        for (;;) {
            int count = *runs;
            if (count <= 0) {
                break;
            }
            int aa = *antialias;
            if (aa) {
                if (aa == 255) {
                    // cool, have the shader draw right into the device
                    shader->shadeSpan(x, y, device, count);
                } else {
                    shader->shadeSpan(x, y, span, count);
                    fProc32Blend(device, span, count, aa);
                }
            }
            device += count;
            runs += count;
            antialias += count;
            x += count;
        }
    } else {    // no xfermode but the shader not opaque
        for (;;) {
            int count = *runs;
            if (count <= 0) {
                break;
            }
            int aa = *antialias;
            if (aa) {
                fShader->shadeSpan(x, y, span, count);
                if (aa == 255) {
                    fProc32(device, span, count, 255);
                } else {
                    fProc32Blend(device, span, count, aa);
                }
            }
            device += count;
            runs += count;
            antialias += count;
            x += count;
        }
    }
}

void SkARGB32_Shader_Blitter::blitMask(const SkMask& mask, const SkIRect& clip) {
    // we only handle kA8 with an xfermode
    if (fXfermode && (SkMask::kA8_Format != mask.fFormat)) {
        this->INHERITED::blitMask(mask, clip);
        return;
    }

    SkASSERT(mask.fBounds.contains(clip));

    SkBlitMask::RowProc proc = NULL;
    if (!fXfermode) {
        unsigned flags = 0;
        if (fShader->getFlags() & SkShader::kOpaqueAlpha_Flag) {
            flags |= SkBlitMask::kSrcIsOpaque_RowFlag;
        }
        proc = SkBlitMask::RowFactory(SkBitmap::kARGB_8888_Config, mask.fFormat,
                                      (SkBlitMask::RowFlags)flags);
        if (NULL == proc) {
            this->INHERITED::blitMask(mask, clip);
            return;
        }
    }

    const int x = clip.fLeft;
    const int width = clip.width();
    int y = clip.fTop;
    int height = clip.height();

    char* dstRow = (char*)fDevice.getAddr32(x, y);
    const size_t dstRB = fDevice.rowBytes();
    const uint8_t* maskRow = (const uint8_t*)mask.getAddr(x, y);
    const size_t maskRB = mask.fRowBytes;

    SkShader* shader = fShader;
    SkPMColor* span = fBuffer;

    if (fXfermode) {
        SkASSERT(SkMask::kA8_Format == mask.fFormat);
        SkXfermode* xfer = fXfermode;
        do {
            shader->shadeSpan(x, y, span, width);
            xfer->xfer32((SkPMColor*)dstRow, span, width, maskRow);
            dstRow += dstRB;
            maskRow += maskRB;
            y += 1;
        } while (--height > 0);
    } else {
        do {
            shader->shadeSpan(x, y, span, width);
            proc(dstRow, maskRow, span, width);
            dstRow += dstRB;
            maskRow += maskRB;
            y += 1;
        } while (--height > 0);
    }
}

