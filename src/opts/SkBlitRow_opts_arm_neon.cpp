/*
 * Copyright 2012 The Android Open Source Project
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 */

#include "SkBlitRow_opts_arm_neon.h"

#include "SkBlitMask.h"
#include "SkBlitRow.h"
#include "SkColorPriv.h"
#include "SkDither.h"
#include "SkMathPriv.h"
#include "SkUtils.h"

#include "SkCachePreload_arm.h"
#include "SkColor_opts_neon.h"
#include <arm_neon.h>

void S32_D565_Opaque_neon(uint16_t* SK_RESTRICT dst,
                           const SkPMColor* SK_RESTRICT src, int count,
                           U8CPU alpha, int /*x*/, int /*y*/) {
    SkASSERT(255 == alpha);

    while (count >= 8) {
        uint8x8x4_t vsrc;
        uint16x8_t vdst;

        // Load
        vsrc = vld4_u8((uint8_t*)src);

        // Convert src to 565
        vdst = SkPixel32ToPixel16_neon8(vsrc);

        // Store
        vst1q_u16(dst, vdst);

        // Prepare next iteration
        dst += 8;
        src += 8;
        count -= 8;
    };

    // Leftovers
    while (count > 0) {
        SkPMColor c = *src++;
        SkPMColorAssert(c);
        *dst = SkPixel32ToPixel16_ToU16(c);
        dst++;
        count--;
    };
}

void S32A_D565_Opaque_neon(uint16_t* SK_RESTRICT dst,
                           const SkPMColor* SK_RESTRICT src, int count,
                           U8CPU alpha, int /*x*/, int /*y*/) {
    SkASSERT(255 == alpha);

    if (count >= 8) {
        uint16_t* SK_RESTRICT keep_dst = 0;

        asm volatile (
                      "ands       ip, %[count], #7            \n\t"
                      "vmov.u8    d31, #1<<7                  \n\t"
                      "vld1.16    {q12}, [%[dst]]             \n\t"
                      "vld4.8     {d0-d3}, [%[src]]           \n\t"
                      // Thumb does not support the standard ARM conditional
                      // instructions but instead requires the 'it' instruction
                      // to signal conditional execution
                      "it eq                                  \n\t"
                      "moveq      ip, #8                      \n\t"
                      "mov        %[keep_dst], %[dst]         \n\t"

                      "add        %[src], %[src], ip, LSL#2   \n\t"
                      "add        %[dst], %[dst], ip, LSL#1   \n\t"
                      "subs       %[count], %[count], ip      \n\t"
                      "b          9f                          \n\t"
                      // LOOP
                      "2:                                         \n\t"

                      "vld1.16    {q12}, [%[dst]]!            \n\t"
                      "vld4.8     {d0-d3}, [%[src]]!          \n\t"
                      "vst1.16    {q10}, [%[keep_dst]]        \n\t"
                      "sub        %[keep_dst], %[dst], #8*2   \n\t"
                      "subs       %[count], %[count], #8      \n\t"
                      "9:                                         \n\t"
                      "pld        [%[dst],#32]                \n\t"
                      // expand 0565 q12 to 8888 {d4-d7}
                      "vmovn.u16  d4, q12                     \n\t"
                      "vshr.u16   q11, q12, #5                \n\t"
                      "vshr.u16   q10, q12, #6+5              \n\t"
                      "vmovn.u16  d5, q11                     \n\t"
                      "vmovn.u16  d6, q10                     \n\t"
                      "vshl.u8    d4, d4, #3                  \n\t"
                      "vshl.u8    d5, d5, #2                  \n\t"
                      "vshl.u8    d6, d6, #3                  \n\t"

                      "vmovl.u8   q14, d31                    \n\t"
                      "vmovl.u8   q13, d31                    \n\t"
                      "vmovl.u8   q12, d31                    \n\t"

                      // duplicate in 4/2/1 & 8pix vsns
                      "vmvn.8     d30, d3                     \n\t"
                      "vmlal.u8   q14, d30, d6                \n\t"
                      "vmlal.u8   q13, d30, d5                \n\t"
                      "vmlal.u8   q12, d30, d4                \n\t"
                      "vshr.u16   q8, q14, #5                 \n\t"
                      "vshr.u16   q9, q13, #6                 \n\t"
                      "vaddhn.u16 d6, q14, q8                 \n\t"
                      "vshr.u16   q8, q12, #5                 \n\t"
                      "vaddhn.u16 d5, q13, q9                 \n\t"
                      "vqadd.u8   d6, d6, d0                  \n\t"  // moved up
                      "vaddhn.u16 d4, q12, q8                 \n\t"
                      // intentionally don't calculate alpha
                      // result in d4-d6

                      "vqadd.u8   d5, d5, d1                  \n\t"
                      "vqadd.u8   d4, d4, d2                  \n\t"

                      // pack 8888 {d4-d6} to 0565 q10
                      "vshll.u8   q10, d6, #8                 \n\t"
                      "vshll.u8   q3, d5, #8                  \n\t"
                      "vshll.u8   q2, d4, #8                  \n\t"
                      "vsri.u16   q10, q3, #5                 \n\t"
                      "vsri.u16   q10, q2, #11                \n\t"

                      "bne        2b                          \n\t"

                      "1:                                         \n\t"
                      "vst1.16      {q10}, [%[keep_dst]]      \n\t"
                      : [count] "+r" (count)
                      : [dst] "r" (dst), [keep_dst] "r" (keep_dst), [src] "r" (src)
                      : "ip", "cc", "memory", "d0","d1","d2","d3","d4","d5","d6","d7",
                      "d16","d17","d18","d19","d20","d21","d22","d23","d24","d25","d26","d27","d28","d29",
                      "d30","d31"
                      );
    }
    else
    {   // handle count < 8
        uint16_t* SK_RESTRICT keep_dst = 0;

        asm volatile (
                      "vmov.u8    d31, #1<<7                  \n\t"
                      "mov        %[keep_dst], %[dst]         \n\t"

                      "tst        %[count], #4                \n\t"
                      "beq        14f                         \n\t"
                      "vld1.16    {d25}, [%[dst]]!            \n\t"
                      "vld1.32    {q1}, [%[src]]!             \n\t"

                      "14:                                        \n\t"
                      "tst        %[count], #2                \n\t"
                      "beq        12f                         \n\t"
                      "vld1.32    {d24[1]}, [%[dst]]!         \n\t"
                      "vld1.32    {d1}, [%[src]]!             \n\t"

                      "12:                                        \n\t"
                      "tst        %[count], #1                \n\t"
                      "beq        11f                         \n\t"
                      "vld1.16    {d24[1]}, [%[dst]]!         \n\t"
                      "vld1.32    {d0[1]}, [%[src]]!          \n\t"

                      "11:                                        \n\t"
                      // unzips achieve the same as a vld4 operation
                      "vuzpq.u16  q0, q1                      \n\t"
                      "vuzp.u8    d0, d1                      \n\t"
                      "vuzp.u8    d2, d3                      \n\t"
                      // expand 0565 q12 to 8888 {d4-d7}
                      "vmovn.u16  d4, q12                     \n\t"
                      "vshr.u16   q11, q12, #5                \n\t"
                      "vshr.u16   q10, q12, #6+5              \n\t"
                      "vmovn.u16  d5, q11                     \n\t"
                      "vmovn.u16  d6, q10                     \n\t"
                      "vshl.u8    d4, d4, #3                  \n\t"
                      "vshl.u8    d5, d5, #2                  \n\t"
                      "vshl.u8    d6, d6, #3                  \n\t"

                      "vmovl.u8   q14, d31                    \n\t"
                      "vmovl.u8   q13, d31                    \n\t"
                      "vmovl.u8   q12, d31                    \n\t"

                      // duplicate in 4/2/1 & 8pix vsns
                      "vmvn.8     d30, d3                     \n\t"
                      "vmlal.u8   q14, d30, d6                \n\t"
                      "vmlal.u8   q13, d30, d5                \n\t"
                      "vmlal.u8   q12, d30, d4                \n\t"
                      "vshr.u16   q8, q14, #5                 \n\t"
                      "vshr.u16   q9, q13, #6                 \n\t"
                      "vaddhn.u16 d6, q14, q8                 \n\t"
                      "vshr.u16   q8, q12, #5                 \n\t"
                      "vaddhn.u16 d5, q13, q9                 \n\t"
                      "vqadd.u8   d6, d6, d0                  \n\t"  // moved up
                      "vaddhn.u16 d4, q12, q8                 \n\t"
                      // intentionally don't calculate alpha
                      // result in d4-d6

                      "vqadd.u8   d5, d5, d1                  \n\t"
                      "vqadd.u8   d4, d4, d2                  \n\t"

                      // pack 8888 {d4-d6} to 0565 q10
                      "vshll.u8   q10, d6, #8                 \n\t"
                      "vshll.u8   q3, d5, #8                  \n\t"
                      "vshll.u8   q2, d4, #8                  \n\t"
                      "vsri.u16   q10, q3, #5                 \n\t"
                      "vsri.u16   q10, q2, #11                \n\t"

                      // store
                      "tst        %[count], #4                \n\t"
                      "beq        24f                         \n\t"
                      "vst1.16    {d21}, [%[keep_dst]]!       \n\t"

                      "24:                                        \n\t"
                      "tst        %[count], #2                \n\t"
                      "beq        22f                         \n\t"
                      "vst1.32    {d20[1]}, [%[keep_dst]]!    \n\t"

                      "22:                                        \n\t"
                      "tst        %[count], #1                \n\t"
                      "beq        21f                         \n\t"
                      "vst1.16    {d20[1]}, [%[keep_dst]]!    \n\t"

                      "21:                                        \n\t"
                      : [count] "+r" (count)
                      : [dst] "r" (dst), [keep_dst] "r" (keep_dst), [src] "r" (src)
                      : "ip", "cc", "memory", "d0","d1","d2","d3","d4","d5","d6","d7",
                      "d16","d17","d18","d19","d20","d21","d22","d23","d24","d25","d26","d27","d28","d29",
                      "d30","d31"
                      );
    }
}

void S32A_D565_Blend_neon(uint16_t* SK_RESTRICT dst,
                          const SkPMColor* SK_RESTRICT src, int count,
                          U8CPU alpha, int /*x*/, int /*y*/) {

    U8CPU alpha_for_asm = alpha;

    asm volatile (
    /* This code implements a Neon version of S32A_D565_Blend. The output differs from
     * the original in two respects:
     *  1. The results have a few mismatches compared to the original code. These mismatches
     *     never exceed 1. It's possible to improve accuracy vs. a floating point
     *     implementation by introducing rounding right shifts (vrshr) for the final stage.
     *     Rounding is not present in the code below, because although results would be closer
     *     to a floating point implementation, the number of mismatches compared to the
     *     original code would be far greater.
     *  2. On certain inputs, the original code can overflow, causing colour channels to
     *     mix. Although the Neon code can also overflow, it doesn't allow one colour channel
     *     to affect another.
     */

#if 1
        /* reflects SkAlpha255To256()'s change from a+a>>7 to a+1 */
                  "add        %[alpha], %[alpha], #1         \n\t"   // adjust range of alpha 0-256
#else
                  "add        %[alpha], %[alpha], %[alpha], lsr #7    \n\t"   // adjust range of alpha 0-256
#endif
                  "vmov.u16   q3, #255                        \n\t"   // set up constant
                  "movs       r4, %[count], lsr #3            \n\t"   // calc. count>>3
                  "vmov.u16   d2[0], %[alpha]                 \n\t"   // move alpha to Neon
                  "beq        2f                              \n\t"   // if count8 == 0, exit
                  "vmov.u16   q15, #0x1f                      \n\t"   // set up blue mask

                  "1:                                             \n\t"
                  "vld1.u16   {d0, d1}, [%[dst]]              \n\t"   // load eight dst RGB565 pixels
                  "subs       r4, r4, #1                      \n\t"   // decrement loop counter
                  "vld4.u8    {d24, d25, d26, d27}, [%[src]]! \n\t"   // load eight src ABGR32 pixels
                  //  and deinterleave

                  "vshl.u16   q9, q0, #5                      \n\t"   // shift green to top of lanes
                  "vand       q10, q0, q15                    \n\t"   // extract blue
                  "vshr.u16   q8, q0, #11                     \n\t"   // extract red
                  "vshr.u16   q9, q9, #10                     \n\t"   // extract green
                  // dstrgb = {q8, q9, q10}

                  "vshr.u8    d24, d24, #3                    \n\t"   // shift red to 565 range
                  "vshr.u8    d25, d25, #2                    \n\t"   // shift green to 565 range
                  "vshr.u8    d26, d26, #3                    \n\t"   // shift blue to 565 range

                  "vmovl.u8   q11, d24                        \n\t"   // widen red to 16 bits
                  "vmovl.u8   q12, d25                        \n\t"   // widen green to 16 bits
                  "vmovl.u8   q14, d27                        \n\t"   // widen alpha to 16 bits
                  "vmovl.u8   q13, d26                        \n\t"   // widen blue to 16 bits
                  // srcrgba = {q11, q12, q13, q14}

                  "vmul.u16   q2, q14, d2[0]                  \n\t"   // sa * src_scale
                  "vmul.u16   q11, q11, d2[0]                 \n\t"   // red result = src_red * src_scale
                  "vmul.u16   q12, q12, d2[0]                 \n\t"   // grn result = src_grn * src_scale
                  "vmul.u16   q13, q13, d2[0]                 \n\t"   // blu result = src_blu * src_scale

                  "vshr.u16   q2, q2, #8                      \n\t"   // sa * src_scale >> 8
                  "vsub.u16   q2, q3, q2                      \n\t"   // 255 - (sa * src_scale >> 8)
                  // dst_scale = q2

                  "vmla.u16   q11, q8, q2                     \n\t"   // red result += dst_red * dst_scale
                  "vmla.u16   q12, q9, q2                     \n\t"   // grn result += dst_grn * dst_scale
                  "vmla.u16   q13, q10, q2                    \n\t"   // blu result += dst_blu * dst_scale

#if 1
    // trying for a better match with SkDiv255Round(a)
    // C alg is:  a+=128; (a+a>>8)>>8
    // we'll use just a rounding shift [q2 is available for scratch]
                  "vrshr.u16   q11, q11, #8                    \n\t"   // shift down red
                  "vrshr.u16   q12, q12, #8                    \n\t"   // shift down green
                  "vrshr.u16   q13, q13, #8                    \n\t"   // shift down blue
#else
    // arm's original "truncating divide by 256"
                  "vshr.u16   q11, q11, #8                    \n\t"   // shift down red
                  "vshr.u16   q12, q12, #8                    \n\t"   // shift down green
                  "vshr.u16   q13, q13, #8                    \n\t"   // shift down blue
#endif

                  "vsli.u16   q13, q12, #5                    \n\t"   // insert green into blue
                  "vsli.u16   q13, q11, #11                   \n\t"   // insert red into green/blue
                  "vst1.16    {d26, d27}, [%[dst]]!           \n\t"   // write pixel back to dst, update ptr

                  "bne        1b                              \n\t"   // if counter != 0, loop
                  "2:                                             \n\t"   // exit

                  : [src] "+r" (src), [dst] "+r" (dst), [count] "+r" (count), [alpha] "+r" (alpha_for_asm)
                  :
                  : "cc", "memory", "r4", "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d16", "d17", "d18", "d19", "d20", "d21", "d22", "d23", "d24", "d25", "d26", "d27", "d28", "d29", "d30", "d31"
                  );

    count &= 7;
    if (count > 0) {
        do {
            SkPMColor sc = *src++;
            if (sc) {
                uint16_t dc = *dst;
                unsigned dst_scale = 255 - SkMulDiv255Round(SkGetPackedA32(sc), alpha);
                unsigned dr = SkMulS16(SkPacked32ToR16(sc), alpha) + SkMulS16(SkGetPackedR16(dc), dst_scale);
                unsigned dg = SkMulS16(SkPacked32ToG16(sc), alpha) + SkMulS16(SkGetPackedG16(dc), dst_scale);
                unsigned db = SkMulS16(SkPacked32ToB16(sc), alpha) + SkMulS16(SkGetPackedB16(dc), dst_scale);
                *dst = SkPackRGB16(SkDiv255Round(dr), SkDiv255Round(dg), SkDiv255Round(db));
            }
            dst += 1;
        } while (--count != 0);
    }
}

/* dither matrix for Neon, derived from gDitherMatrix_3Bit_16.
 * each dither value is spaced out into byte lanes, and repeated
 * to allow an 8-byte load from offsets 0, 1, 2 or 3 from the
 * start of each row.
 */
static const uint8_t gDitherMatrix_Neon[48] = {
    0, 4, 1, 5, 0, 4, 1, 5, 0, 4, 1, 5,
    6, 2, 7, 3, 6, 2, 7, 3, 6, 2, 7, 3,
    1, 5, 0, 4, 1, 5, 0, 4, 1, 5, 0, 4,
    7, 3, 6, 2, 7, 3, 6, 2, 7, 3, 6, 2,

};

void S32_D565_Blend_Dither_neon(uint16_t *dst, const SkPMColor *src,
                                int count, U8CPU alpha, int x, int y)
{

    SkASSERT(255 > alpha);

    // rescale alpha to range 1 - 256
    int scale = SkAlpha255To256(alpha);

    if (count >= 8) {
        /* select row and offset for dither array */
        const uint8_t *dstart = &gDitherMatrix_Neon[(y&3)*12 + (x&3)];

        uint8x8_t vdither = vld1_u8(dstart);         // load dither values
        uint8x8_t vdither_g = vshr_n_u8(vdither, 1); // calc. green dither values

        int16x8_t vscale = vdupq_n_s16(scale);        // duplicate scale into neon reg
        uint16x8_t vmask_b = vdupq_n_u16(0x1F);         // set up blue mask

        do {

            uint8x8_t vsrc_r, vsrc_g, vsrc_b;
            uint8x8_t vsrc565_r, vsrc565_g, vsrc565_b;
            uint16x8_t vsrc_dit_r, vsrc_dit_g, vsrc_dit_b;
            uint16x8_t vsrc_res_r, vsrc_res_g, vsrc_res_b;
            uint16x8_t vdst;
            uint16x8_t vdst_r, vdst_g, vdst_b;
            int16x8_t vres_r, vres_g, vres_b;
            int8x8_t vres8_r, vres8_g, vres8_b;

            // Load source and add dither
            {
            register uint8x8_t d0 asm("d0");
            register uint8x8_t d1 asm("d1");
            register uint8x8_t d2 asm("d2");
            register uint8x8_t d3 asm("d3");

            asm (
                "vld4.8    {d0-d3},[%[src]]!  /* r=%P0 g=%P1 b=%P2 a=%P3 */"
                : "=w" (d0), "=w" (d1), "=w" (d2), "=w" (d3), [src] "+&r" (src)
                :
            );
            vsrc_g = d1;
#if SK_PMCOLOR_BYTE_ORDER(B,G,R,A)
            vsrc_r = d2; vsrc_b = d0;
#elif SK_PMCOLOR_BYTE_ORDER(R,G,B,A)
            vsrc_r = d0; vsrc_b = d2;
#endif
            }

            vsrc565_g = vshr_n_u8(vsrc_g, 6); // calc. green >> 6
            vsrc565_r = vshr_n_u8(vsrc_r, 5); // calc. red >> 5
            vsrc565_b = vshr_n_u8(vsrc_b, 5); // calc. blue >> 5

            vsrc_dit_g = vaddl_u8(vsrc_g, vdither_g); // add in dither to green and widen
            vsrc_dit_r = vaddl_u8(vsrc_r, vdither);   // add in dither to red and widen
            vsrc_dit_b = vaddl_u8(vsrc_b, vdither);   // add in dither to blue and widen

            vsrc_dit_r = vsubw_u8(vsrc_dit_r, vsrc565_r);  // sub shifted red from result
            vsrc_dit_g = vsubw_u8(vsrc_dit_g, vsrc565_g);  // sub shifted green from result
            vsrc_dit_b = vsubw_u8(vsrc_dit_b, vsrc565_b);  // sub shifted blue from result

            vsrc_res_r = vshrq_n_u16(vsrc_dit_r, 3);
            vsrc_res_g = vshrq_n_u16(vsrc_dit_g, 2);
            vsrc_res_b = vshrq_n_u16(vsrc_dit_b, 3);

            // Load dst and unpack
            vdst = vld1q_u16(dst);
            vdst_g = vshrq_n_u16(vdst, 5);                   // shift down to get green
            vdst_r = vshrq_n_u16(vshlq_n_u16(vdst, 5), 5+5); // double shift to extract red
            vdst_b = vandq_u16(vdst, vmask_b);               // mask to get blue

            // subtract dst from src and widen
            vres_r = vsubq_s16(vreinterpretq_s16_u16(vsrc_res_r), vreinterpretq_s16_u16(vdst_r));
            vres_g = vsubq_s16(vreinterpretq_s16_u16(vsrc_res_g), vreinterpretq_s16_u16(vdst_g));
            vres_b = vsubq_s16(vreinterpretq_s16_u16(vsrc_res_b), vreinterpretq_s16_u16(vdst_b));

            // multiply diffs by scale and shift
            vres_r = vmulq_s16(vres_r, vscale);
            vres_g = vmulq_s16(vres_g, vscale);
            vres_b = vmulq_s16(vres_b, vscale);

            vres8_r = vshrn_n_s16(vres_r, 8);
            vres8_g = vshrn_n_s16(vres_g, 8);
            vres8_b = vshrn_n_s16(vres_b, 8);

            // add dst to result
            vres_r = vaddw_s8(vreinterpretq_s16_u16(vdst_r), vres8_r);
            vres_g = vaddw_s8(vreinterpretq_s16_u16(vdst_g), vres8_g);
            vres_b = vaddw_s8(vreinterpretq_s16_u16(vdst_b), vres8_b);

            // put result into 565 format
            vres_b = vsliq_n_s16(vres_b, vres_g, 5);   // shift up green and insert into blue
            vres_b = vsliq_n_s16(vres_b, vres_r, 6+5); // shift up red and insert into blue

            // Store result
            vst1q_u16(dst, vreinterpretq_u16_s16(vres_b));

            // Next iteration
            dst += 8;
            count -= 8;

        } while (count >= 8);
    }

    // Leftovers
    if (count > 0) {
        int scale = SkAlpha255To256(alpha);
        DITHER_565_SCAN(y);
        do {
            SkPMColor c = *src++;
            SkPMColorAssert(c);

            int dither = DITHER_VALUE(x);
            int sr = SkGetPackedR32(c);
            int sg = SkGetPackedG32(c);
            int sb = SkGetPackedB32(c);
            sr = SkDITHER_R32To565(sr, dither);
            sg = SkDITHER_G32To565(sg, dither);
            sb = SkDITHER_B32To565(sb, dither);

            uint16_t d = *dst;
            *dst++ = SkPackRGB16(SkAlphaBlend(sr, SkGetPackedR16(d), scale),
                                 SkAlphaBlend(sg, SkGetPackedG16(d), scale),
                                 SkAlphaBlend(sb, SkGetPackedB16(d), scale));
            DITHER_INC_X(x);
        } while (--count != 0);
    }
}

void S32A_Opaque_BlitRow32_neon(SkPMColor* SK_RESTRICT dst,
                                const SkPMColor* SK_RESTRICT src,
                                int count, U8CPU alpha) {

    SkASSERT(255 == alpha);
    if (count > 0) {


    uint8x8_t alpha_mask;

    static const uint8_t alpha_mask_setup[] = {3,3,3,3,7,7,7,7};
    alpha_mask = vld1_u8(alpha_mask_setup);

    /* do the NEON unrolled code */
#define    UNROLL    4
    while (count >= UNROLL) {
        uint8x8_t src_raw, dst_raw, dst_final;
        uint8x8_t src_raw_2, dst_raw_2, dst_final_2;

        /* The two prefetches below may make the code slighlty
         * slower for small values of count but are worth having
         * in the general case.
         */
        __builtin_prefetch(src+32);
        __builtin_prefetch(dst+32);

        /* get the source */
        src_raw = vreinterpret_u8_u32(vld1_u32(src));
#if    UNROLL > 2
        src_raw_2 = vreinterpret_u8_u32(vld1_u32(src+2));
#endif

        /* get and hold the dst too */
        dst_raw = vreinterpret_u8_u32(vld1_u32(dst));
#if    UNROLL > 2
        dst_raw_2 = vreinterpret_u8_u32(vld1_u32(dst+2));
#endif

    /* 1st and 2nd bits of the unrolling */
    {
        uint8x8_t dst_cooked;
        uint16x8_t dst_wide;
        uint8x8_t alpha_narrow;
        uint16x8_t alpha_wide;

        /* get the alphas spread out properly */
        alpha_narrow = vtbl1_u8(src_raw, alpha_mask);
        alpha_wide = vsubw_u8(vdupq_n_u16(256), alpha_narrow);

        /* spread the dest */
        dst_wide = vmovl_u8(dst_raw);

        /* alpha mul the dest */
        dst_wide = vmulq_u16 (dst_wide, alpha_wide);
        dst_cooked = vshrn_n_u16(dst_wide, 8);

        /* sum -- ignoring any byte lane overflows */
        dst_final = vadd_u8(src_raw, dst_cooked);
    }

#if    UNROLL > 2
    /* the 3rd and 4th bits of our unrolling */
    {
        uint8x8_t dst_cooked;
        uint16x8_t dst_wide;
        uint8x8_t alpha_narrow;
        uint16x8_t alpha_wide;

        alpha_narrow = vtbl1_u8(src_raw_2, alpha_mask);
        alpha_wide = vsubw_u8(vdupq_n_u16(256), alpha_narrow);

        /* spread the dest */
        dst_wide = vmovl_u8(dst_raw_2);

        /* alpha mul the dest */
        dst_wide = vmulq_u16 (dst_wide, alpha_wide);
        dst_cooked = vshrn_n_u16(dst_wide, 8);

        /* sum -- ignoring any byte lane overflows */
        dst_final_2 = vadd_u8(src_raw_2, dst_cooked);
    }
#endif

        vst1_u32(dst, vreinterpret_u32_u8(dst_final));
#if    UNROLL > 2
        vst1_u32(dst+2, vreinterpret_u32_u8(dst_final_2));
#endif

        src += UNROLL;
        dst += UNROLL;
        count -= UNROLL;
    }
#undef    UNROLL

    /* do any residual iterations */
        while (--count >= 0) {
            *dst = SkPMSrcOver(*src, *dst);
            src += 1;
            dst += 1;
        }
    }
}

extern "C" {
    void S32A_Opaque_BlitRow32_neon_o(SkPMColor* SK_RESTRICT dst,
                                const SkPMColor* SK_RESTRICT src,
                                int count, U8CPU alpha);
}
void S32A_Opaque_BlitRow32_neon_src_alpha(SkPMColor* SK_RESTRICT dst,
                                const SkPMColor* SK_RESTRICT src,
                                int count, U8CPU alpha) {
    SkASSERT(255 == alpha);

    if (count <= 0)
    return;

    /* Use these to check if src is transparent or opaque */
    const unsigned int ALPHA_OPAQ  = 0xFF000000;
    const unsigned int ALPHA_TRANS = 0x00FFFFFF;

#define UNROLL  4
    const SkPMColor* SK_RESTRICT src_end = src + count - (UNROLL + 1);
    const SkPMColor* SK_RESTRICT src_temp = src;

    /* set up the NEON variables */
    uint8x8_t alpha_mask;
    static const uint8_t alpha_mask_setup[] = {3,3,3,3,7,7,7,7};
    alpha_mask = vld1_u8(alpha_mask_setup);

    uint8x8_t src_raw, dst_raw, dst_final;
    uint8x8_t src_raw_2, dst_raw_2, dst_final_2;
    uint8x8_t dst_cooked;
    uint16x8_t dst_wide;
    uint8x8_t alpha_narrow;
    uint16x8_t alpha_wide;

    /* choose the first processing type */
    if( src >= src_end)
        goto TAIL;
    if(*src <= ALPHA_TRANS)
        goto ALPHA_0;
    if(*src >= ALPHA_OPAQ)
        goto ALPHA_255;
    /* fall-thru */

ALPHA_1_TO_254:
    do {

        /* get the source */
        src_raw = vreinterpret_u8_u32(vld1_u32(src));
        src_raw_2 = vreinterpret_u8_u32(vld1_u32(src+2));

        /* get and hold the dst too */
        dst_raw = vreinterpret_u8_u32(vld1_u32(dst));
        dst_raw_2 = vreinterpret_u8_u32(vld1_u32(dst+2));


        /* get the alphas spread out properly */
        alpha_narrow = vtbl1_u8(src_raw, alpha_mask);
        /* reflect SkAlpha255To256() semantics a+1 vs a+a>>7 */
        /* we collapsed (255-a)+1 ... */
        alpha_wide = vsubw_u8(vdupq_n_u16(256), alpha_narrow);

        /* spread the dest */
        dst_wide = vmovl_u8(dst_raw);

        /* alpha mul the dest */
        dst_wide = vmulq_u16 (dst_wide, alpha_wide);
        dst_cooked = vshrn_n_u16(dst_wide, 8);

        /* sum -- ignoring any byte lane overflows */
        dst_final = vadd_u8(src_raw, dst_cooked);

        alpha_narrow = vtbl1_u8(src_raw_2, alpha_mask);
        /* reflect SkAlpha255To256() semantics a+1 vs a+a>>7 */
        /* we collapsed (255-a)+1 ... */
        alpha_wide = vsubw_u8(vdupq_n_u16(256), alpha_narrow);

        /* spread the dest */
        dst_wide = vmovl_u8(dst_raw_2);

        /* alpha mul the dest */
        dst_wide = vmulq_u16 (dst_wide, alpha_wide);
        dst_cooked = vshrn_n_u16(dst_wide, 8);

        /* sum -- ignoring any byte lane overflows */
        dst_final_2 = vadd_u8(src_raw_2, dst_cooked);

        vst1_u32(dst, vreinterpret_u32_u8(dst_final));
        vst1_u32(dst+2, vreinterpret_u32_u8(dst_final_2));

        src += UNROLL;
        dst += UNROLL;

        /* if 2 of the next pixels aren't between 1 and 254
        it might make sense to go to the optimized loops */
        if((src[0] <= ALPHA_TRANS && src[1] <= ALPHA_TRANS) || (src[0] >= ALPHA_OPAQ && src[1] >= ALPHA_OPAQ))
            break;

    } while(src < src_end);

    if (src >= src_end)
        goto TAIL;

    if(src[0] >= ALPHA_OPAQ && src[1] >= ALPHA_OPAQ)
        goto ALPHA_255;

    /*fall-thru*/

ALPHA_0:

    /*In this state, we know the current alpha is 0 and
     we optimize for the next alpha also being zero. */
    src_temp = src;  //so we don't have to increment dst every time
    do {
        if(*(++src) > ALPHA_TRANS)
            break;
        if(*(++src) > ALPHA_TRANS)
            break;
        if(*(++src) > ALPHA_TRANS)
            break;
        if(*(++src) > ALPHA_TRANS)
            break;
    } while(src < src_end);

    dst += (src - src_temp);

    /* no longer alpha 0, so determine where to go next. */
    if( src >= src_end)
        goto TAIL;
    if(*src >= ALPHA_OPAQ)
        goto ALPHA_255;
    else
        goto ALPHA_1_TO_254;

ALPHA_255:
    while((src[0] & src[1] & src[2] & src[3]) >= ALPHA_OPAQ) {
        dst[0]=src[0];
        dst[1]=src[1];
        dst[2]=src[2];
        dst[3]=src[3];
        src+=UNROLL;
        dst+=UNROLL;
        if(src >= src_end)
            goto TAIL;
    }

    //Handle remainder.
    if(*src >= ALPHA_OPAQ) { *dst++ = *src++;
        if(*src >= ALPHA_OPAQ) { *dst++ = *src++;
            if(*src >= ALPHA_OPAQ) { *dst++ = *src++; }
        }
    }

    if( src >= src_end)
        goto TAIL;
    if(*src <= ALPHA_TRANS)
        goto ALPHA_0;
    else
        goto ALPHA_1_TO_254;

TAIL:
    /* do any residual iterations */
    src_end += UNROLL + 1;  //goto the real end
    while(src != src_end) {
        if( *src != 0 ) {
            if( *src >= ALPHA_OPAQ ) {
                *dst = *src;
            }
            else {
                *dst = SkPMSrcOver(*src, *dst);
            }
        }
        src++;
        dst++;
    }

#undef    UNROLL
    return;
}

/* Neon version of S32_Blend_BlitRow32()
 * portable version is in src/core/SkBlitRow_D32.cpp
 */
void S32_Blend_BlitRow32_neon(SkPMColor* SK_RESTRICT dst,
                              const SkPMColor* SK_RESTRICT src,
                              int count, U8CPU alpha) {
    SkASSERT(alpha <= 255);
    if (count > 0) {
        uint16_t src_scale = SkAlpha255To256(alpha);
        uint16_t dst_scale = 256 - src_scale;

    /* run them N at a time through the NEON unit */
    /* note that each 1 is 4 bytes, each treated exactly the same,
     * so we can work under that guise. We *do* know that the src&dst
     * will be 32-bit aligned quantities, so we can specify that on
     * the load/store ops and do a neon 'reinterpret' to get us to
     * byte-sized (pun intended) pieces that we widen/multiply/shift
     * we're limited at 128 bits in the wide ops, which is 8x16bits
     * or a pair of 32 bit src/dsts.
     */
    /* we *could* manually unroll this loop so that we load 128 bits
     * (as a pair of 64s) from each of src and dst, processing them
     * in pieces. This might give us a little better management of
     * the memory latency, but my initial attempts here did not
     * produce an instruction stream that looked all that nice.
     */
#define    UNROLL    2
    while (count >= UNROLL) {
        uint8x8_t  src_raw, dst_raw, dst_final;
        uint16x8_t  src_wide, dst_wide;

        /* get 64 bits of src, widen it, multiply by src_scale */
        src_raw = vreinterpret_u8_u32(vld1_u32(src));
        src_wide = vmovl_u8(src_raw);
        /* gcc hoists vdupq_n_u16(), better than using vmulq_n_u16() */
        src_wide = vmulq_u16 (src_wide, vdupq_n_u16(src_scale));

        /* ditto with dst */
        dst_raw = vreinterpret_u8_u32(vld1_u32(dst));
        dst_wide = vmovl_u8(dst_raw);

        /* combine add with dst multiply into mul-accumulate */
        dst_wide = vmlaq_u16(src_wide, dst_wide, vdupq_n_u16(dst_scale));

        dst_final = vshrn_n_u16(dst_wide, 8);
        vst1_u32(dst, vreinterpret_u32_u8(dst_final));

        src += UNROLL;
        dst += UNROLL;
        count -= UNROLL;
    }
    /* RBE: well, i don't like how gcc manages src/dst across the above
     * loop it's constantly calculating src+bias, dst+bias and it only
     * adjusts the real ones when we leave the loop. Not sure why
     * it's "hoisting down" (hoisting implies above in my lexicon ;))
     * the adjustments to src/dst/count, but it does...
     * (might be SSA-style internal logic...
     */

#if    UNROLL == 2
    if (count == 1) {
            *dst = SkAlphaMulQ(*src, src_scale) + SkAlphaMulQ(*dst, dst_scale);
    }
#else
    if (count > 0) {
            do {
                *dst = SkAlphaMulQ(*src, src_scale) + SkAlphaMulQ(*dst, dst_scale);
                src += 1;
                dst += 1;
            } while (--count > 0);
    }
#endif

#undef    UNROLL
    }
}

void S32A_Blend_BlitRow32_neon(SkPMColor* SK_RESTRICT dst,
                         const SkPMColor* SK_RESTRICT src,
                         int count, U8CPU alpha) {

    SkASSERT(255 >= alpha);

    if (count <= 0) {
        return;
    }

    unsigned alpha256 = SkAlpha255To256(alpha);

    // First deal with odd counts
    if (count & 1) {
        uint8x8_t vsrc = vdup_n_u8(0), vdst = vdup_n_u8(0), vres;
        uint16x8_t vdst_wide, vsrc_wide;
        unsigned dst_scale;

        // Load
        vsrc = vreinterpret_u8_u32(vld1_lane_u32(src, vreinterpret_u32_u8(vsrc), 0));
        vdst = vreinterpret_u8_u32(vld1_lane_u32(dst, vreinterpret_u32_u8(vdst), 0));

        // Calc dst_scale
        dst_scale = vget_lane_u8(vsrc, 3);
        dst_scale *= alpha256;
        dst_scale >>= 8;
        dst_scale = 256 - dst_scale;

        // Process src
        vsrc_wide = vmovl_u8(vsrc);
        vsrc_wide = vmulq_n_u16(vsrc_wide, alpha256);

        // Process dst
        vdst_wide = vmovl_u8(vdst);
        vdst_wide = vmulq_n_u16(vdst_wide, dst_scale);

        // Combine
        vres = vshrn_n_u16(vdst_wide, 8) + vshrn_n_u16(vsrc_wide, 8);

        vst1_lane_u32(dst, vreinterpret_u32_u8(vres), 0);
        dst++;
        src++;
        count--;
    }

    if (count) {
        uint8x8_t alpha_mask;
        static const uint8_t alpha_mask_setup[] = {3,3,3,3,7,7,7,7};
        alpha_mask = vld1_u8(alpha_mask_setup);

        do {

            uint8x8_t vsrc, vdst, vres, vsrc_alphas;
            uint16x8_t vdst_wide, vsrc_wide, vsrc_scale, vdst_scale;

            __builtin_prefetch(src+32);
            __builtin_prefetch(dst+32);

            // Load
            vsrc = vreinterpret_u8_u32(vld1_u32(src));
            vdst = vreinterpret_u8_u32(vld1_u32(dst));

            // Prepare src_scale
            vsrc_scale = vdupq_n_u16(alpha256);

            // Calc dst_scale
            vsrc_alphas = vtbl1_u8(vsrc, alpha_mask);
            vdst_scale = vmovl_u8(vsrc_alphas);
            vdst_scale *= vsrc_scale;
            vdst_scale = vshrq_n_u16(vdst_scale, 8);
            vdst_scale = vsubq_u16(vdupq_n_u16(256), vdst_scale);

            // Process src
            vsrc_wide = vmovl_u8(vsrc);
            vsrc_wide *= vsrc_scale;

            // Process dst
            vdst_wide = vmovl_u8(vdst);
            vdst_wide *= vdst_scale;

            // Combine
            vres = vshrn_n_u16(vdst_wide, 8) + vshrn_n_u16(vsrc_wide, 8);

            vst1_u32(dst, vreinterpret_u32_u8(vres));

            src += 2;
            dst += 2;
            count -= 2;
        } while(count);
    }
}

///////////////////////////////////////////////////////////////////////////////

#undef    DEBUG_OPAQUE_DITHER

#if    defined(DEBUG_OPAQUE_DITHER)
static void showme8(char *str, void *p, int len)
{
    static char buf[256];
    char tbuf[32];
    int i;
    char *pc = (char*) p;
    sprintf(buf,"%8s:", str);
    for(i=0;i<len;i++) {
        sprintf(tbuf, "   %02x", pc[i]);
        strcat(buf, tbuf);
    }
    SkDebugf("%s\n", buf);
}
static void showme16(char *str, void *p, int len)
{
    static char buf[256];
    char tbuf[32];
    int i;
    uint16_t *pc = (uint16_t*) p;
    sprintf(buf,"%8s:", str);
    len = (len / sizeof(uint16_t));    /* passed as bytes */
    for(i=0;i<len;i++) {
        sprintf(tbuf, " %04x", pc[i]);
        strcat(buf, tbuf);
    }
    SkDebugf("%s\n", buf);
}
#endif

void S32A_D565_Opaque_Dither_neon (uint16_t * SK_RESTRICT dst,
                                   const SkPMColor* SK_RESTRICT src,
                                   int count, U8CPU alpha, int x, int y) {
    SkASSERT(255 == alpha);
    /* residuals */
#if 0
    if (count > 0) {
        DITHER_565_SCAN(y);
        do {
            SkPMColor c = *src++;
            SkPMColorAssert(c);
            if (c) {
                unsigned a = SkGetPackedA32(c);

                // dither and alpha are just temporary variables to work-around
                // an ICE in debug.
                unsigned dither = DITHER_VALUE(x);
                unsigned alpha = SkAlpha255To256(a);
                int d = SkAlphaMul(dither, alpha);

                unsigned sr = SkGetPackedR32(c);
                unsigned sg = SkGetPackedG32(c);
                unsigned sb = SkGetPackedB32(c);
                sr = SkDITHER_R32_FOR_565(sr, d);
                sg = SkDITHER_G32_FOR_565(sg, d);
                sb = SkDITHER_B32_FOR_565(sb, d);

                uint32_t src_expanded = (sg << 24) | (sr << 13) | (sb << 2);
                uint32_t dst_expanded = SkExpand_rgb_16(*dst);
                dst_expanded = dst_expanded * (SkAlpha255To256(255 - a) >> 3);
                // now src and dst expanded are in g:11 r:10 x:1 b:10
                *dst = SkCompact_rgb_16((src_expanded + dst_expanded) >> 5);
            }
            dst += 1;
            DITHER_INC_X(x);
        } while (--count != 0);
    }
#else
    const uint8_t *dstart = &gDitherMatrix_Neon[(y&3)*12 + (x&3)];
    asm volatile (
        "cmp        r2, #0                              \n"
        "bxeq       lr                                  \n"
        "cmp        %[count], #8                        \n"
        "blt        .Lless_than8                        \n"
        "vld1.8     {d31}, [%[dstart]]                  \n" // d31 = dither[x]
        "vmov.i8    d29, #0xff                          \n" 
        "vmovl.u8   q15, d31                            \n" // wide of dither 
        "vmov.i8    d28, #1                             \n" 
        "vmov.i16   q7, #0x00f8                         \n"
        "vmov.i16   q12, #0x00fc                        \n"
    ".LMain_loop8:                                      \n" 
        "vld4.8     {d0, d1, d2, d3}, [%[src]]!         \n" // d0 = [R], d1 = [G], d2 = [B], d3 = [A]
        "sub        %[count], %[count], #8              \n"
        "vmov       r9, r10, s6, s7                     \n" // compare flag
        "cmp        r9, #0xffffffff                     \n"
        "cmpeq      r10, #0xffffffff                    \n" // alpha is eq to 255
        "vaddl.u8   q13, d28, d3                        \n" // q13 = SkAlpha255To256(a)
        "vshr.u8    d8,  d0, #5                         \n" // d8 = r >> 5
        "vshr.u8    d9,  d1, #6                         \n" // d9 = g >> 6
        "vmul.i16   q13, q13, q15                       \n" // q13 = SkAlphaMul(dither, alpha)
        "vshr.u8    d10, d2, #5                         \n" // d10 = b >> 5
        "vsub.i8    d4,  d0, d8                         \n" // d0 = r - (r >> 5)
        "vsub.i8    d5,  d2, d10                        \n" // d2 = b - (b >> 5)
        "vshrn.u16  d27, q13, #8                        \n" // d27 = SkAlphaMul(dither, alpha);
        "vsub.i8    d6,  d1, d9                         \n" // d1 = g - (g >> 6)
        "vshr.u8    d26, d27, #1                        \n" // d27 = d >> 1
        "vaddl.u8   q4,  d4, d27                        \n" // q4 = SkDITHER_R32_FOR_565(r, d), r + d - (r >> 5)
        "vaddl.u8   q5,  d5, d27                        \n" // q5 = SkDITHER_R32_FOR_565(b, d), b + d - (b >> 5)
        "vaddl.u8   q6,  d6, d26                        \n" // q6 = SkDITHER_G32_FOR_565(g, d),(g + (d >> 1) - (g >> 6))
        "beq        .Lalpha_255                         \n"

        /*
         * alpha blend
         */
        "vld1.16    {d4, d5}, [%[dst]]                  \n" // d4, d5
        "vsub.i8    d27, d29, d3                        \n" // (255 - a)
        "vaddl.u8   q13, d27, d28                       \n" // SkAlpha255To256(255 - a)
        "vshr.u16   q0,  q2, #8                         \n" // q0 = [dR]
        "vshr.u16   q1,  q2, #3                         \n" // q1 = [dG]
        "vshl.u16   q2,  q2, #3                         \n" // q2 = [dB]
        "vand.i16   q0,  q0, q7                         \n"
        "vand.i16   q1,  q1, q12                        \n"
        "vand.i16   q2,  q2, q7                         \n"
        "vmul.i16   q0,  q0, q13                        \n" // q0 = [dr] * SkAlpha255To256(255 - a)
        "vmul.i16   q1,  q1, q13                        \n"
        "vmul.i16   q2,  q2, q13                        \n"
        "vsra.u16   q4,  q0, #8                         \n" // q4 = sr + dr
        "vsra.u16   q5,  q2, #8                         \n" // q5 = sb + db 
        "vsra.u16   q6,  q1, #8                         \n" // q6 = sg + dg
    ".Lalpha_255:                                       \n"
        "vshr.u16   q6,  q6, #2                         \n"
        "vshr.u16   q4,  q4, #3                         \n"
        "vshr.u16   q2,  q5, #3                         \n" // b
        "vsli.16    q2,  q6, #5                         \n"
        "vsli.16    q2,  q4, #11                        \n"
        "cmp        %[count], #8                        \n"
        "vst1.16    {d4, d5}, [%[dst]]!                 \n"
        "bge        .LMain_loop8                        \n"

    ".Lless_than8:                                      \n"
        :[dst] "+r" (dst), [src] "+r" (src), [count] "+r" (count), [dstart] "+r" (dstart)
        :
        : "memory", "cc", "r9", "r10"
    );
    /* residuals */
    if (count > 0) {
        DITHER_565_SCAN(y);
        do {
            SkPMColor c = *src++;
            SkPMColorAssert(c);
            if (c) {
                unsigned a = SkGetPackedA32(c);

                // dither and alpha are just temporary variables to work-around
                // an ICE in debug.
                unsigned dither = DITHER_VALUE(x);
                unsigned alpha = SkAlpha255To256(a);
                int d = SkAlphaMul(dither, alpha);

                unsigned sr = SkGetPackedR32(c);
                unsigned sg = SkGetPackedG32(c);
                unsigned sb = SkGetPackedB32(c);
                sr = SkDITHER_R32_FOR_565(sr, d);
                sg = SkDITHER_G32_FOR_565(sg, d);
                sb = SkDITHER_B32_FOR_565(sb, d);

                uint32_t src_expanded = (sg << 24) | (sr << 13) | (sb << 2);
                uint32_t dst_expanded = SkExpand_rgb_16(*dst);
                dst_expanded = dst_expanded * (SkAlpha255To256(255 - a) >> 3);
                // now src and dst expanded are in g:11 r:10 x:1 b:10
                *dst = SkCompact_rgb_16((src_expanded + dst_expanded) >> 5);
            }
            dst += 1;
            DITHER_INC_X(x);
        } while (--count != 0);
    }
#endif
}

///////////////////////////////////////////////////////////////////////////////

#undef    DEBUG_S32_OPAQUE_DITHER

void S32_D565_Opaque_Dither_neon(uint16_t* SK_RESTRICT dst,
                                 const SkPMColor* SK_RESTRICT src,
                                 int count, U8CPU alpha, int x, int y) {
    SkASSERT(255 == alpha);

#define    UNROLL    8
    if (count >= UNROLL) {
    uint8x8_t d;
    const uint8_t *dstart = &gDitherMatrix_Neon[(y&3)*12 + (x&3)];
    d = vld1_u8(dstart);

    while (count >= UNROLL) {
        uint8x8_t sr, sg, sb;
        uint16x8_t dr, dg, db;
        uint16x8_t dst8;

        {
        register uint8x8_t d0 asm("d0");
        register uint8x8_t d1 asm("d1");
        register uint8x8_t d2 asm("d2");
        register uint8x8_t d3 asm("d3");

        asm (
            "vld4.8    {d0-d3},[%[src]]!  /* r=%P0 g=%P1 b=%P2 a=%P3 */"
            : "=w" (d0), "=w" (d1), "=w" (d2), "=w" (d3), [src] "+&r" (src)
            :
        );
        sg = d1;
#if SK_PMCOLOR_BYTE_ORDER(B,G,R,A)
        sr = d2; sb = d0;
#elif SK_PMCOLOR_BYTE_ORDER(R,G,B,A)
        sr = d0; sb = d2;
#endif
        }
        /* XXX: if we want to prefetch, hide it in the above asm()
         * using the gcc __builtin_prefetch(), the prefetch will
         * fall to the bottom of the loop -- it won't stick up
         * at the top of the loop, just after the vld4.
         */

        // sr = sr - (sr>>5) + d
        sr = vsub_u8(sr, vshr_n_u8(sr, 5));
        dr = vaddl_u8(sr, d);

        // sb = sb - (sb>>5) + d
        sb = vsub_u8(sb, vshr_n_u8(sb, 5));
        db = vaddl_u8(sb, d);

        // sg = sg - (sg>>6) + d>>1; similar logic for overflows
        sg = vsub_u8(sg, vshr_n_u8(sg, 6));
        dg = vaddl_u8(sg, vshr_n_u8(d, 1));

        // pack high bits of each into 565 format  (rgb, b is lsb)
        dst8 = vshrq_n_u16(db, 3);
        dst8 = vsliq_n_u16(dst8, vshrq_n_u16(dg, 2), 5);
        dst8 = vsliq_n_u16(dst8, vshrq_n_u16(dr, 3), 11);

        // store it
        vst1q_u16(dst, dst8);

#if    defined(DEBUG_S32_OPAQUE_DITHER)
        // always good to know if we generated good results
        {
        int i, myx = x, myy = y;
        DITHER_565_SCAN(myy);
        for (i=0;i<UNROLL;i++) {
            // the '!' in the asm block above post-incremented src by the 8 pixels it reads.
            SkPMColor c = src[i-8];
            unsigned dither = DITHER_VALUE(myx);
            uint16_t val = SkDitherRGB32To565(c, dither);
            if (val != dst[i]) {
            SkDebugf("RBE: src %08x dither %02x, want %04x got %04x dbas[i] %02x\n",
                c, dither, val, dst[i], dstart[i]);
            }
            DITHER_INC_X(myx);
        }
        }
#endif

        dst += UNROLL;
        // we don't need to increment src as the asm above has already done it
        count -= UNROLL;
        x += UNROLL;        // probably superfluous
    }
    }
#undef    UNROLL

    // residuals
    if (count > 0) {
        DITHER_565_SCAN(y);
        do {
            SkPMColor c = *src++;
            SkPMColorAssert(c);
            SkASSERT(SkGetPackedA32(c) == 255);

            unsigned dither = DITHER_VALUE(x);
            *dst++ = SkDitherRGB32To565(c, dither);
            DITHER_INC_X(x);
        } while (--count != 0);
    }
}

void Color32_arm_neon(SkPMColor* dst, const SkPMColor* src, int count,
                      SkPMColor color) {
    if (count <= 0) {
        return;
    }

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

        if (count >= 8) {
            // at the end of this assembly, count will have been decremented
            // to a negative value. That is, if count mod 8 = x, it will be
            // -8 +x coming out.
            asm volatile (
                PLD128(src, 0)

                "vdup.32    q0, %[color]                \n\t"

                PLD128(src, 128)

                // scale numerical interval [0-255], so load as 8 bits
                "vdup.8     d2, %[scale]                \n\t"

                PLD128(src, 256)

                "subs       %[count], %[count], #8      \n\t"

                PLD128(src, 384)

                "Loop_Color32:                          \n\t"

                // load src color, 8 pixels, 4 64 bit registers
                // (and increment src).
                "vld1.32    {d4-d7}, [%[src]]!          \n\t"

                PLD128(src, 384)

                // multiply long by scale, 64 bits at a time,
                // destination into a 128 bit register.
                "vmull.u8   q4, d4, d2                  \n\t"
                "vmull.u8   q5, d5, d2                  \n\t"
                "vmull.u8   q6, d6, d2                  \n\t"
                "vmull.u8   q7, d7, d2                  \n\t"

                // shift the 128 bit registers, containing the 16
                // bit scaled values back to 8 bits, narrowing the
                // results to 64 bit registers.
                "vshrn.i16  d8, q4, #8                  \n\t"
                "vshrn.i16  d9, q5, #8                  \n\t"
                "vshrn.i16  d10, q6, #8                 \n\t"
                "vshrn.i16  d11, q7, #8                 \n\t"

                // adding back the color, using 128 bit registers.
                "vadd.i8    q6, q4, q0                  \n\t"
                "vadd.i8    q7, q5, q0                  \n\t"

                // store back the 8 calculated pixels (2 128 bit
                // registers), and increment dst.
                "vst1.32    {d12-d15}, [%[dst]]!        \n\t"

                "subs       %[count], %[count], #8      \n\t"
                "bge        Loop_Color32                \n\t"
                : [src] "+r" (src), [dst] "+r" (dst), [count] "+r" (count)
                : [color] "r" (color), [scale] "r" (scale)
                : "cc", "memory",
                  "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7",
                  "d8", "d9", "d10", "d11", "d12", "d13", "d14", "d15"
                          );
            // At this point, if we went through the inline assembly, count is
            // a negative value:
            // if the value is -8, there is no pixel left to process.
            // if the value is -7, there is one pixel left to process
            // ...
            // And'ing it with 7 will give us the number of pixels
            // left to process.
            count = count & 0x7;
        }

        while (count > 0) {
            *dst = color + SkAlphaMulQ(*src, scale);
            src += 1;
            dst += 1;
            count--;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

const SkBlitRow::Proc sk_blitrow_platform_565_procs_arm_neon[] = {
    // no dither
    // NOTE: For the S32_D565_Blend function below, we don't have a special
    //       version that assumes that each source pixel is opaque. But our
    //       S32A is still faster than the default, so use it.
    S32_D565_Opaque_neon,
    S32A_D565_Blend_neon,   // really S32_D565_Blend
    S32A_D565_Opaque_neon,
    S32A_D565_Blend_neon,

    // dither
    S32_D565_Opaque_Dither_neon,
    S32_D565_Blend_Dither_neon,
    S32A_D565_Opaque_Dither_neon,
    NULL,   // S32A_D565_Blend_Dither
};

const SkBlitRow::Proc32 sk_blitrow_platform_32_procs_arm_neon[] = {
    NULL,   // S32_Opaque,
    S32_Blend_BlitRow32_neon,        // S32_Blend,
    /*
     * We have two choices for S32A_Opaque procs. The one reads the src alpha
     * value and attempts to optimize accordingly.  The optimization is
     * sensitive to the source content and is not a win in all cases. For
     * example, if there are a lot of transitions between the alpha states,
     * the performance will almost certainly be worse.  However, for many
     * common cases the performance is equivalent or better than the standard
     * case where we do not inspect the src alpha.
     */
#if SK_A32_SHIFT == 24
    // This proc assumes the alpha value occupies bits 24-32 of each SkPMColor
    //S32A_Opaque_BlitRow32_neon_src_alpha,   // S32A_Opaque,
    S32A_Opaque_BlitRow32_neon_o,
#else
    S32A_Opaque_BlitRow32_neon,     // S32A_Opaque,
#endif
    S32A_Blend_BlitRow32_neon        // S32A_Blend
};
