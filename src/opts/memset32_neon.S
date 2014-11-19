/***************************************************************************
 * Copyright (c) 2009,2010, Code Aurora Forum. All rights reserved.
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 ***************************************************************************/

    .code 32
    .fpu neon
    .align 8

    .global k_bzero 
    .type   k_bzero, %function

    /*
     * tao.zeng, most calling case, the input number is large enough
     * implement sk_bzero
     */
     .func
k_bzero:
    cmp         r1, #0                                                  // check count
    bxeq        lr
    ands        r12, r0, #0x03                                          // test if aligned 4
    mov         r2, r1                                                  // backup count
    mov         r1, #0
    beq         sk_bzero_align_4

    tst         r0, #1                                                  // aligned to 2 bytes
    strneb      r1, [r0], #1
    subs        r2, r2, #1
    bxeq        lr
    cmp         r2, #2
    bge         align_to_4

    strb        r1, [r0]                                                // only remain 1 bytes
    bx          lr

align_to_4:
    tst         r0, #2
    strgeh      r1, [r0], #2
    subs        r2, r2, #2
    bxeq        lr

sk_bzero_align_4:
    mov         r12, r2
    lsr         r2,  #2
    mov         r3,  lr
    bl          memset32_neon                                           // call memset32_neon

    ands        r12, r12, #3                                            // process remain data
    bxeq        r3
    cmp         r12, #2
    strgeh      r1, [r0], #2
    subs        r12, #2
    bxeq        r3
    strb        r1, [r0]
    bx          r3
    
    .endfunc

    .global memset32_neon
    .type   memset32_neon, %function
    .func
memset32_neon:
    cmp         r2, #0
    bxeq        lr
    vdup.32     q0, r1                                                  // copy data
    vdup.32     q1, r1
    tst         r0, #0x04
    strne       r1, [r0], #4                                            // aligned to 8 bytes
    subne       r2, r2, #1

    cmp         r2, #128                                                // main loop is 64
    blt         memset32_less_128
    sub         r2, #128
memset32_loop128:
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    subs        r2, r2, #128
    bge         memset32_loop128
    add         r2, r2, #128
memset32_less_128:
    cmp         r2, #0
    bxeq        lr
    cmp         r2, #64
    blt         memset32_less_64
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    subs        r2, r2, #64
    bxeq        lr
memset32_less_64:
    cmp         r2, #32
    blt         memset32_less_32
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    subs        r2, r2, #32
    bxeq        lr
memset32_less_32:
#if 0
    cmp         r2, #16
    vstmge      r0!, {q0, q1}
    vstmge      r0!, {q0, q1}
    subs        r2, r2, #16
    bxeq        lr
    cmp         r2, #8
    vstmge      r0!, {q0, q1}
    subs        r2, r2, #8
    bxeq        lr
    cmp         r2, #4
    vstmge      r0!, {q0}
    subs        r2, r2, #4
    bxeq        lr
    cmp         r2, #2
    vstmge      r0!, {d0}
    subs        r2, #2
    bxeq        lr
    str         r1, [r0], #4
    bx          lr
#else
    cmp         r2, #16
    blt         memset32_less_16
    vstmia      r0!, {q0, q1}
    vstmia      r0!, {q0, q1}
    subs        r2, r2, #16
    bxeq        lr
memset32_less_16:
    cmp         r2, #8
    blt         memset32_less_8
    vstmia      r0!, {q0, q1}
    subs        r2, r2, #8
    bxeq        lr
memset32_less_8:
    cmp         r2, #4
    blt         memset32_less_4
    vstmia      r0!, {q0}
    subs        r2, r2, #4
    bxeq        lr
memset32_less_4:
    cmp         r2, #2
    blt         memset32_less_2
    vstmia      r0!, {d0}
    subs        r2, r2, #2
    bxeq        lr
memset32_less_2:
    str         r1, [r0], #4
    bx          lr
#endif
    .endfunc
    .end
