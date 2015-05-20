/*
 * Copyright 2011 Google Inc.
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 */

#include "SkBlitRow.h"
#include "SkBlitMask.h"
#include "SkColorPriv.h"
#include "SkUtils.h"

#define UNROLL

SkBlitRow::ColorRectProc PlatformColorRectProcFactory();

static void S32_Opaque_BlitRow32(SkPMColor* SK_RESTRICT dst,
                                 const SkPMColor* SK_RESTRICT src,
                                 int count, U8CPU alpha) {
#if defined(__aarch64__)
    /*
     * TODO: optimize for AARCH64
     */
    SkASSERT(255 == alpha);
    sk_memcpy32(dst, src, count);
#else
    /*
     * tao.zeng@amlogic.com, src and dst are algined by 4
     */
    asm volatile (
        "cmp        %[count], #0                    \n"
        "it         eq                              \n"
        "bxeq       lr                              \n"                 // if count == 0, return
        "pld        [%[src], #32]                   \n"
        "tst        %[src], #0x04                   \n"                 // make aligne to 8 bytes
        "ittt       ne                              \n"
        "ldrne      r12, [%[src]], #4               \n"
        "strne      r12, [%[dst]], #4               \n"
        "subne      %[count], %[count], #1          \n"
        "cmp        %[count], #16                   \n"                 //
        "blt        S32_Opaque_BlitRow32_less16     \n"                 //
        "pld        [%[src], #64]                   \n"
        "sub        %[count], #16                   \n"
    "S32_Opaque_BlitRow32_loop16:                   \n"
        "vldmia     %[src]!, {q0, q1, q2, q3}       \n"
        "pld        [%[src], #64]                   \n"
        "pld        [%[src], #96]                   \n"
        "subs       %[count], %[count], #16         \n"
        "vstmia     %[dst]!, {q0, q1, q2, q3}       \n"
        "bge        S32_Opaque_BlitRow32_loop16     \n"
        "adds       %[count], %[count], #16         \n"
        "cmp        %[count], #0                    \n"
        "it         eq                              \n"
        "bxeq       lr                              \n"

    "S32_Opaque_BlitRow32_less16:                   \n"
        "cmp        %[count], #8                    \n"
        "blt        S32_Opaque_BlitRow32_less8      \n"
        "vldmia     %[src]!, {q0, q1}               \n"
        "subs       %[count], %[count], #8          \n"
        "vstmia     %[dst]!, {q0, q1}               \n"
        "it         eq                              \n"
        "bxeq       lr                              \n"
    "S32_Opaque_BlitRow32_less8:                    \n"
        "cmp        %[count], #4                    \n"
        "blt        S32_Opaque_BlitRow32_less4      \n"
        "vldmia     %[src]!, {d0, d1}               \n"
        "subs       %[count], %[count], #4          \n"
        "vstmia     %[dst]!, {d0, d1}               \n"
        "it         eq                              \n"
        "bxeq       lr                              \n"
    "S32_Opaque_BlitRow32_less4:                    \n"
        "cmp        %[count], #2                    \n"
        "itt        ge                              \n"
        "ldmge      %[src]!, {%[alpha], r12}        \n"
        "subges     %[count], #2                    \n"
        "it         ge                              \n"
        "stmge      %[dst]!, {%[alpha], r12}        \n"
        "it         eq                              \n"
        "bxeq       lr                              \n"
        "ldr        r12, [%[src]], #4               \n"
        "str        r12, [%[dst]], #4               \n"
        "bx         lr                              \n"
        :
        :[src] "r" (src), [dst] "r" (dst), [count] "r" (count), [alpha] "r" (alpha)
        :"cc", "memory", "r12"
    );
#endif
}

static void S32_Blend_BlitRow32(SkPMColor* SK_RESTRICT dst,
                                const SkPMColor* SK_RESTRICT src,
                                int count, U8CPU alpha) {
    SkASSERT(alpha <= 255);
    if (count > 0) {
        unsigned src_scale = SkAlpha255To256(alpha);
        unsigned dst_scale = 256 - src_scale;

#ifdef UNROLL
        if (count & 1) {
            *dst = SkAlphaMulQ(*(src++), src_scale) + SkAlphaMulQ(*dst, dst_scale);
            dst += 1;
            count -= 1;
        }

        const SkPMColor* SK_RESTRICT srcEnd = src + count;
        while (src != srcEnd) {
            *dst = SkAlphaMulQ(*(src++), src_scale) + SkAlphaMulQ(*dst, dst_scale);
            dst += 1;
            *dst = SkAlphaMulQ(*(src++), src_scale) + SkAlphaMulQ(*dst, dst_scale);
            dst += 1;
        }
#else
        do {
            *dst = SkAlphaMulQ(*src, src_scale) + SkAlphaMulQ(*dst, dst_scale);
            src += 1;
            dst += 1;
        } while (--count > 0);
#endif
    }
}

static void S32A_Opaque_BlitRow32(SkPMColor* SK_RESTRICT dst,
                                  const SkPMColor* SK_RESTRICT src,
                                  int count, U8CPU alpha) {
    SkASSERT(255 == alpha);
    if (count > 0) {
#ifdef UNROLL
        if (count & 1) {
            *dst = SkPMSrcOver(*(src++), *dst);
            dst += 1;
            count -= 1;
        }

        const SkPMColor* SK_RESTRICT srcEnd = src + count;
        while (src != srcEnd) {
            *dst = SkPMSrcOver(*(src++), *dst);
            dst += 1;
            *dst = SkPMSrcOver(*(src++), *dst);
            dst += 1;
        }
#else
        do {
            *dst = SkPMSrcOver(*src, *dst);
            src += 1;
            dst += 1;
        } while (--count > 0);
#endif
    }
}

static void S32A_Blend_BlitRow32(SkPMColor* SK_RESTRICT dst,
                                 const SkPMColor* SK_RESTRICT src,
                                 int count, U8CPU alpha) {
    SkASSERT(alpha <= 255);
    if (count > 0) {
#ifdef UNROLL
        if (count & 1) {
            *dst = SkBlendARGB32(*(src++), *dst, alpha);
            dst += 1;
            count -= 1;
        }

        const SkPMColor* SK_RESTRICT srcEnd = src + count;
        while (src != srcEnd) {
            *dst = SkBlendARGB32(*(src++), *dst, alpha);
            dst += 1;
            *dst = SkBlendARGB32(*(src++), *dst, alpha);
            dst += 1;
        }
#else
        do {
            *dst = SkBlendARGB32(*src, *dst, alpha);
            src += 1;
            dst += 1;
        } while (--count > 0);
#endif
    }
}

///////////////////////////////////////////////////////////////////////////////

static const SkBlitRow::Proc32 gDefault_Procs32[] = {
    S32_Opaque_BlitRow32,
    S32_Blend_BlitRow32,
    S32A_Opaque_BlitRow32,
    S32A_Blend_BlitRow32
};

SkBlitRow::Proc32 SkBlitRow::Factory32(unsigned flags) {
    SkASSERT(flags < SK_ARRAY_COUNT(gDefault_Procs32));
    // just so we don't crash
    flags &= kFlags32_Mask;

    SkBlitRow::Proc32 proc = PlatformProcs32(flags);
    if (NULL == proc) {
        proc = gDefault_Procs32[flags];
    }
    SkASSERT(proc);
    return proc;
}

SkBlitRow::Proc32 SkBlitRow::ColorProcFactory() {
    SkBlitRow::ColorProc proc = PlatformColorProc();
    if (NULL == proc) {
        proc = Color32;
    }
    SkASSERT(proc);
    return proc;
}

void SkBlitRow::Color32(SkPMColor* SK_RESTRICT dst,
                        const SkPMColor* SK_RESTRICT src,
                        int count, SkPMColor color) {
    if (count > 0) {
        if (0 == color) {
            if (src != dst) {
                memcpy(dst, src, count * sizeof(SkPMColor));
            }
            return;
        }
        unsigned colorA = SkGetPackedA32(color);
        if (255 == colorA) {
            sk_memset32(dst, color, count);
        } else {
            unsigned scale = 256 - SkAlpha255To256(colorA);
        #if defined(__ARM_HAVE_NEON) && !defined(__aarch64__)
            /*
             * tao.zeng@amlogic.com, use NEON to optimize these funciton
             */
            // tao.zeng, alpha will not be 256, because SkAlpha255To256(colorA);
            asm volatile (
                /*
                 * regigster allocate
                 * q0   --> scale
                 * q1   --> wide color
                 * q2 - q3 --> src0 -- src8
                 */
                "vdup.32    d2, %[color]                        \n"
                "vdup.8     q0, %[scale]                        \n"     // q0 = [scale].8
                "pld        [%[src], #128]                      \n"     // preload data
                "cmp        %[count], #8                        \n"
                "vshll.u8   q1, d2, #8                          \n"     // q1 = [color 00].16
                "blt        SkBlitRow_Color32_less_8            \n"
                // main loop
            "SkBlitRow_Color32_loop_8:                          \n"
                "vldmia     %[src]!, {d4, d5, d6, d7}           \n"     // load 8 src data
                "pld        [%[src], #128]                      \n"
                "sub        %[count], %[count], #8              \n"
                "cmp        %[count], #8                        \n"
                "vmull.u8   q8, d4, d0                          \n"     // multiple scale
                "vmull.u8   q9, d5, d0                          \n"
                "vmull.u8   q10, d6, d0                         \n"
                "vmull.u8   q11, d7, d0                         \n"
                "vadd.u16   q8, q8, q1                          \n"     // add src color
                "vadd.u16   q9, q9, q1                          \n"
                "vadd.u16   q10, q10, q1                        \n"
                "vadd.u16   q11, q11, q1                        \n"
                "vshrn.i16  d4, q8, #8                          \n"     // narrow result
                "vshrn.i16  d5, q9, #8                          \n"
                "vshrn.i16  d6, q10, #8                         \n"
                "vshrn.i16  d7, q11, #8                         \n"
                "vstmia     %[dst]!, {d4, d5, d6, d7}           \n"
                "bge        SkBlitRow_Color32_loop_8            \n"
                "cmp        %[count], #0                        \n"
                "beq        out                                 \n"

            "SkBlitRow_Color32_less_8:                          \n"
                "cmp        %[count], #4                        \n"
                "blt        SkBlitRow_Color32_less_4            \n"

                "vldmia     %[src]!, {d4, d5}                   \n"
                "subs       %[count], %[count], #4              \n"
                "vmull.u8   q8, d4, d0                          \n"     // multiple scale
                "vmull.u8   q9, d5, d0                          \n"
                "vadd.u16   q8, q8, q1                          \n"     // add src color
                "vadd.u16   q9, q9, q1                          \n"
                "vshrn.i16  d4, q8, #8                          \n"     // narrow result
                "vshrn.i16  d5, q9, #8                          \n"
                "vstmia     %[dst]!, {d4, d5}                   \n"
                "beq        out                                 \n"

            "SkBlitRow_Color32_less_4:                          \n"
                "cmp        %[count], #2                        \n"
                "blt        SkBlitRow_Color32_less_2            \n"
                "vldmia     %[src]!, {d4}                       \n"
                "vmull.u8   q8, d4, d0                          \n"     // multiple scale
                "subs       %[count], %[count], #2              \n"
                "vadd.u16   q8, q8, q1                          \n"     // add src color
                "vshrn.i16  d4, q8, #8                          \n"     // narrow result
                "vstmia     %[dst]!, {d4}                       \n"
                "beq        out                                 \n"

            "SkBlitRow_Color32_less_2:                          \n"
                "vld1.32    {d4[0]}, [%[src]]!                  \n"
                "vmull.u8   q8, d4, d0                          \n"     // multiple scale
                "vadd.u16   q8, q8, q1                          \n"     // add src color
                "vshrn.i16  d4, q8, #8                          \n"     // narrow result
                "vst1.32    {d4[0]}, [%[dst]]!                  \n"

            "out:                                               \n"
                :
                : [color] "r"(color), [dst] "r" (dst), [src] "r" (src), [scale] "r" (scale), [count] "r" (count)
                : "memory"
            );
        #else
            /*
             * TODO: optimize for AARCH64
             */
            do {
                *dst = color + SkAlphaMulQ(*src, scale);
                src += 1;
                dst += 1;
            } while (--count);
        #endif
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

static void D32_Mask_Color(void* dst, size_t dstRB, SkBitmap::Config,
                           const uint8_t* mask, size_t maskRB, SkColor color,
                           int width, int height) {
    SkPMColor pmc = SkPreMultiplyColor(color);
    size_t dstOffset = dstRB - (width << 2);
    size_t maskOffset = maskRB - width;
    SkPMColor *device = (SkPMColor *)dst;

    do {
        int w = width;
        do {
            unsigned aa = *mask++;
            *device = SkBlendARGB32(pmc, *device, aa);
            device += 1;
        } while (--w != 0);
        device = (uint32_t*)((char*)device + dstOffset);
        mask += maskOffset;
    } while (--height != 0);
}

static void D32_Mask_Opaque(void* dst, size_t dstRB, SkBitmap::Config,
                            const uint8_t* mask, size_t maskRB, SkColor color,
                            int width, int height) {
    SkPMColor pmc = SkPreMultiplyColor(color);
    uint32_t* device = (uint32_t*)dst;

    maskRB -= width;
    dstRB -= (width << 2);

#if !defined(__aarch64__)
    asm volatile (
        "vmov       s0,  %[pmc]                     \n"     // set ARGB of pmc
        "vdup.32    d31, %[pmc]                     \n"     // d31 = pmc.32
        "vmov.i16   q13, #0x100                     \n"
        "vmovl.u8   q0,  d0                         \n"
        "vmovl.u8   q15, d31                        \n"
        "vmov.i8    d1,  #0x01                      \n"
        :
        : [pmc] "r" (pmc)
        : "cc", "memory"
    );
#endif

    do {
    #if defined(__aarch64__)
        /*
         * TODO: optimize for AARCH64
         */
        int w = width;
        do {
            unsigned aa = *mask++;
            *device = SkAlphaMulQ(pmc, SkAlpha255To256(aa)) + SkAlphaMulQ(*device, SkAlpha255To256(255 - aa));
            device += 1;
        } while (--w != 0);
    #else
        asm volatile (
            "mov        r12, %[width]                           \n"
            "cmp        r12, #8                                 \n"     // count > 8?
            "blt        D32_Mask_Opaque_less_8                  \n"
        "D32_Mask_Opaque_loop8:                                 \n"
            "vld1.8     {d7}, [%[mask]]!                        \n"     // load 8 mask
            "vld4.8     {d2, d3, d4, d5}, [%[device]]           \n"     // load 8 dst data
            "sub        r12, r12, #8                            \n"
            "vaddl.u8   q14, d7,  d1                            \n"     // 255to256(aa)
            "vsubw.u8   q8,  q13, d7                            \n"     // 255to256(255 - aa)
            "vmovl.u8   q9,  d2                                 \n"     // expend R to 16 bits
            "vmovl.u8   q10, d3                                 \n"     // expend G to 16 bits
            "vmovl.u8   q11, d4                                 \n"     // expend B to 16 bits
            "vmovl.u8   q12, d5                                 \n"     // expend A to 16 bits

            "vmul.i16   q9,  q9,  q8                            \n"     // device.R * 255to256(255 -aa)
            "vmul.i16   q10, q10, q8                            \n"     // device.G * 255to256(255 -aa)
            "vmul.i16   q11, q11, q8                            \n"     // device.B * 255to256(255 -aa)
            "vmul.i16   q12, q12, q8                            \n"     // device.A * 255to256(255 -aa)
            "vmul.i16   q1,  q14, d0[0]                         \n"     // pmc.R * 255to256(aa)
            "vmul.i16   q2,  q14, d0[1]                         \n"     // pmc.G * 255to256(aa)
            "vmul.i16   q8,  q14, d0[2]                         \n"     // pmc.B * 255to256(aa)
            "vmul.i16   q14, q14, d0[3]                         \n"     // pmc.A * 255to256(aa)
            "vshrn.i16  d18, q9,  #8                            \n"     //
            "vshrn.i16  d19, q10, #8                            \n"
            "vshrn.i16  d20, q11, #8                            \n"
            "vshrn.i16  d21, q12, #8                            \n"
            "vshrn.i16  d2,  q1,  #8                            \n"
            "vshrn.i16  d3,  q2,  #8                            \n"
            "vshrn.i16  d4,  q8,  #8                            \n"
            "vshrn.i16  d5,  q14, #8                            \n"
            "vadd.i8    q1,  q1,  q9                            \n"
            "vadd.i8    q2,  q2,  q10                           \n"
            "cmp        r12, #8                                 \n"
            "vst4.8     {d2, d3, d4, d5}, [%[device]]!          \n"
            "bge        D32_Mask_Opaque_loop8                   \n"
            "cmp        r12, #0                                 \n"
            "beq        D32_Mask_Opaque_out                     \n"

        "D32_Mask_Opaque_less_8:                                \n"
            "ldr        r11, =D32_Mask_Opaque_table             \n"
            "cmp        r12, #4                                 \n"
            "vldmia     r11, {q14}                              \n"
            "blt        D32_Mask_Opaque_less_4                  \n"
            "vld1.32    {d7[]}, [%[mask]]!                      \n"     // d7  = [aa 0..3 0..3]
            "vldmia     %[device], {d2, d3}                     \n"
            "subs       r12, r12, #4                            \n"
            "vtbl.8     d4,  {d7}, d28                          \n"     // d4  = [1 1 1 1 0 0 0 0] of aa
            "vtbl.8     d5,  {d7}, d29                          \n"     // d5  = [3 3 3 3 2 2 2 2] of aa
            "vmovl.u8   q9,  d2                                 \n"     // exprand to 16 bits
            "vmovl.u8   q10, d3                                 \n"
            "vaddl.u8   q11, d4,  d1                            \n"     // 255to256(aa)
            "vaddl.u8   q12, d5,  d1                            \n"     // 255to256(aa)
            "vsubw.u8   q3,  q13, d4                            \n"     // 255to256(255 - aa)
            "vsubw.u8   q8,  q13, d5                            \n"     // 255to256(255 - aa)
            "vmul.i16   q9,  q9,  q3                            \n"     // device * 255to256(255 -aa)
            "vmul.i16   q10, q10, q8                            \n"     // device * 255to256(255 -aa)
            "vmul.i16   q11, q15, q11                           \n"     // pmc * 255to256(aa)
            "vmul.i16   q12, q15, q12                           \n"     // pmc * 255to256(aa)
            "vshrn.i16  d2,  q9,  #8                            \n"
            "vshrn.i16  d3,  q10, #8                            \n"
            "vshrn.i16  d4,  q11, #8                            \n"
            "vshrn.i16  d5,  q12, #8                            \n"
            "vadd.i8    q1,  q1,  q2                            \n"
            "vstmia     %[device]!, {d2, d3}                    \n"
            "beq        D32_Mask_Opaque_out                     \n"

        "D32_Mask_Opaque_less_4:                                \n"
            "cmp        r12, #2                                 \n"
            "blt        D32_Mask_Opaque_less_2                  \n"
            "vld1.16    {d7[]}, [%[mask]]!                      \n"
            "vldmia     %[device], {d2}                         \n"
            "subs       r12, r12, #2                            \n"
            "vtbl.8     d4,  {d7}, d28                          \n"     // d4  = [1 1 1 1 0 0 0 0] of aa
            "vmovl.u8   q9,  d2                                 \n"     // exprand to 16 bits
            "vaddl.u8   q11, d4,  d1                            \n"     // 255to256(aa)
            "vsubw.u8   q3,  q13, d4                            \n"     // 255to256(255 - aa)
            "vmul.i16   q9,  q9,  q3                            \n"     // device * 255to256(255 -aa)
            "vmul.i16   q10, q15, q11                           \n"     // pmc * 255to256(aa)
            "vshrn.i16  d2,  q9,  #8                            \n"
            "vshrn.i16  d3,  q10, #8                            \n"
            "vadd.i8    d2,  d3,  d2                            \n"
            "vstmia     %[device]!, {d2}                        \n"
            "beq        D32_Mask_Opaque_out                     \n"

        "D32_Mask_Opaque_less_2:                                \n"
            "vld1.8     {d7[]}, [%[mask]]!                      \n"
            "vld1.32    {d2[0]}, [%[device]]                    \n"
            "vtbl.8     d4,  {d7}, d28                          \n"     // d4  = [1 1 1 1 0 0 0 0] of aa
            "vaddl.u8   q11, d4,  d1                            \n"     // 255to256(aa)
            "vsubw.u8   q3,  q13, d4                            \n"     // 255to256(255 - aa)
            "vmovl.u8   q9,  d2                                 \n"     // exprand to 16 bits
            "vmul.i16   q9,  q9,  q3                            \n"     // device * 255to256(255 -aa)
            "vmul.i16   q10, q15, q11                           \n"     // pmc * 255to256(aa)
            "vshrn.i16  d2,  q9,  #8                            \n"
            "vshrn.i16  d3,  q10, #8                            \n"
            "vadd.i8    d2,  d3,  d2                            \n"
            "vst1.32    {d2[0]}, [%[device]]!                   \n"

        "D32_Mask_Opaque_out:                                   \n"
            :
            : [device] "r" (device) , [mask] "r" (mask), [width] "r" (width)
            : "cc", "memory", "r12", "r11"
        );
    #endif
        device = (uint32_t*)((char*)device + dstRB);
        mask += maskRB;
    } while (--height != 0);

#if !defined(__aarch64__)
    asm volatile (
    "D32_Mask_Opaque_table:                         \n"
        ".word      0x00000000                      \n"
        ".word      0x01010101                      \n"
        ".word      0x02020202                      \n"
        ".word      0x03030303                      \n"
    );
#endif
}

static void D32_Mask_Black(void* dst, size_t dstRB, SkBitmap::Config,
                           const uint8_t* mask, size_t maskRB, SkColor,
                           int width, int height) {
    uint32_t* device = (uint32_t*)dst;

    maskRB -= width;
    dstRB -= (width << 2);

    do {
        int w = width;
        do {
            unsigned aa = *mask++;
            *device = (aa << SK_A32_SHIFT) + SkAlphaMulQ(*device, SkAlpha255To256(255 - aa));
            device += 1;
        } while (--w != 0);
        device = (uint32_t*)((char*)device + dstRB);
        mask += maskRB;
    } while (--height != 0);
}


template <size_t N> void assignLoop(SkPMColor* dst, SkPMColor color) {
    for (size_t i = 0; i < N; ++i) {
        *dst++ = color;
    }
}

static inline void assignLoop(SkPMColor dst[], SkPMColor color, int count) {
    while (count >= 4) {
        *dst++ = color;
        *dst++ = color;
        *dst++ = color;
        *dst++ = color;
        count -= 4;
    }
    if (count >= 2) {
        *dst++ = color;
        *dst++ = color;
        count -= 2;
    }
    if (count > 0) {
        *dst++ = color;
    }
}

void SkBlitRow::ColorRect32(SkPMColor* dst, int width, int height,
                            size_t rowBytes, SkPMColor color) {
    if (width <= 0 || height <= 0 || 0 == color) {
        return;
    }

    // Just made up this value, since I saw it once in a SSE2 file.
    // We should consider writing some tests to find the optimimal break-point
    // (or query the Platform proc?)
    static const int MIN_WIDTH_FOR_SCANLINE_PROC = 32;

    bool isOpaque = (0xFF == SkGetPackedA32(color));

    if (!isOpaque || width >= MIN_WIDTH_FOR_SCANLINE_PROC) {
        SkBlitRow::ColorProc proc = SkBlitRow::ColorProcFactory();
        while (--height >= 0) {
            (*proc)(dst, dst, width, color);
            dst = (SkPMColor*) ((char*)dst + rowBytes);
        }
    } else {
        switch (width) {
            case 1:
                while (--height >= 0) {
                    assignLoop<1>(dst, color);
                    dst = (SkPMColor*) ((char*)dst + rowBytes);
                }
                break;
            case 2:
                while (--height >= 0) {
                    assignLoop<2>(dst, color);
                    dst = (SkPMColor*) ((char*)dst + rowBytes);
                }
                break;
            case 3:
                while (--height >= 0) {
                    assignLoop<3>(dst, color);
                    dst = (SkPMColor*) ((char*)dst + rowBytes);
                }
                break;
            default:
                while (--height >= 0) {
                    assignLoop(dst, color, width);
                    dst = (SkPMColor*) ((char*)dst + rowBytes);
                }
                break;
        }
    }
}

SkBlitRow::ColorRectProc SkBlitRow::ColorRectProcFactory() {
    SkBlitRow::ColorRectProc proc = PlatformColorRectProcFactory();
    if (NULL == proc) {
        proc = ColorRect32;
    }
    SkASSERT(proc);
    return proc;
}
