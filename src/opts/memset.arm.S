/*
 * Copyright 2010 The Android Open Source Project
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 */

/* Changes:
 * 2010-08-11 Steve McIntyre <steve.mcintyre@arm.com>
 *    Added small changes to the two functions to make them work on the
 *    specified number of 16- or 32-bit values rather than the original
 *    code which was specified as a count of bytes. More verbose comments
 *    to aid future maintenance.
 */

    .text
    .align

    .global arm_memset32
    .type   arm_memset32, %function
    .global arm_memset16
    .type   arm_memset16, %function

/*
 * Optimized memset functions for ARM.
 *
 * void arm_memset16(uint16_t* dst, uint16_t value, int count);
 * void arm_memset32(uint32_t* dst, uint32_t value, int count);
 *
 */
arm_memset16:
        .fnstart
        cmp         r2, #0x0
        bxle        lr
        push        {lr}

        /* expand the data to 32 bits */
        orr         r1, r1, lsl #16

        /* align to 32 bits */
        tst         r0, #2
        strneh      r1, [r0], #2
        subne       r2, r2, #1

        /* Multiply count by 2 - go from the number of 16-bit shorts
         * to the number of bytes desired. */
        movs        r2, r2, lsl #1

        popeq       {pc}
        /* Now jump into the main loop below. */
        b           .Lwork_32
        .fnend

arm_memset32:
        .fnstart
        cmp         r2, #0x0
        bxle        lr
        push        {lr}

        /* Multiply count by 4 - go from the number of 32-bit words to
         * the number of bytes desired. */
        mov         r2, r2, lsl #2

.Lwork_32:
        /* Set up registers ready for writing them out. */
        mov         ip, r1
        mov         lr, r1

        /* Try to align the destination to a cache line. Assume 32
         * byte (8 word) cache lines, it's the common case. */
        cmp         r2, #32
        mov         r3, r1
        blt         .Lmemset32_less32
        tst         r0, #4
        beq         .Laligned32
        str         r1, [r0], #4
        sub         r2, r2, #4
        cmp         r2, #32
        blt         .Lmemset32_less32

        /* Now quickly loop through the cache-aligned data. */
.Laligned32:
1:      subs        r2, r2, #32
        stmhsia     r0!, {r1,r3,ip,lr}
        stmhsia     r0!, {r1,r3,ip,lr}
        bhs         1b
        add         r2, r2, #32

        /* (Optionally) store any remaining trailing bytes.
         * (0-30 bytes, length in r2) */
.Lmemset32_less32:
        movs        r2, r2, lsl #28
        stmcsia     r0!, {r1,r3,ip,lr}
        stmmiia     r0!, {r1,lr}
        movs        r2, r2, lsl #2
        strcs       r1, [r0], #4
        strmih      lr, [r0], #2

        pop         {pc}
        .fnend
