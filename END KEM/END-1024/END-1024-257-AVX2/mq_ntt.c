#include <string.h>
#include <stdio.h>
#include <immintrin.h>
#include "mq_ntt_param.h"

static void montgomery_reduce_257(__m256i *r, __m256i *a, __m256i *b)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(257);
    __m256i tmp_QINV = _mm256_set1_epi16(-255);

    tmp_c0 = _mm256_mullo_epi16(*a, *b);
    tmp_c1 = _mm256_mulhi_epi16(*a, *b);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *r     = _mm256_sub_epi16(tmp_c1, tmp_x);
}

void mq_poly_mul_ntt_257(int16_t *r, int16_t *a, int16_t *b, int16_t *muldata)
{
    int i;
    __m256i zeta, tmp1,tmp2;
    __m256i a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7;
    // __m256i a8, a9, a10, a11, a12, a13, a14, a15, b8, b9, b10, b11, b12, b13, b14, b15;
    __m256i r0, r1, r2, r3, r4, r5, r6, r7;
    for(i = 0;i < DIM_N;i += 256)
    {
        zeta = _mm256_loadu_si256((__m256i *)(muldata + i / 16));

        a0 = _mm256_loadu_si256((__m256i *)(a + i));
        a1 = _mm256_loadu_si256((__m256i *)(a + i + 16));
        a2 = _mm256_loadu_si256((__m256i *)(a + i + 32));
        a3 = _mm256_loadu_si256((__m256i *)(a + i + 48));
        a4 = _mm256_loadu_si256((__m256i *)(a + i + 128));
        a5 = _mm256_loadu_si256((__m256i *)(a + i + 144));
        a6 = _mm256_loadu_si256((__m256i *)(a + i + 160));
        a7 = _mm256_loadu_si256((__m256i *)(a + i + 176));

        b0 = _mm256_loadu_si256((__m256i *)(b + i));
        b1 = _mm256_loadu_si256((__m256i *)(b + i + 16));
        b2 = _mm256_loadu_si256((__m256i *)(b + i + 32));
        b3 = _mm256_loadu_si256((__m256i *)(b + i + 48));
        b4 = _mm256_loadu_si256((__m256i *)(b + i + 128));
        b5 = _mm256_loadu_si256((__m256i *)(b + i + 144));
        b6 = _mm256_loadu_si256((__m256i *)(b + i + 160));
        b7 = _mm256_loadu_si256((__m256i *)(b + i + 176));

        montgomery_reduce_257(&r0, &a0, &b0);
        montgomery_reduce_257(&tmp1, &a7, &b1);
        montgomery_reduce_257(&tmp2, &a6, &b2);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r0 = _mm256_add_epi16(r0, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i), r0);

        montgomery_reduce_257(&r1, &a1, &b0);
        montgomery_reduce_257(&tmp2, &a0, &b1);
        r1 = _mm256_add_epi16(r1, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b2);
        montgomery_reduce_257(&tmp2, &a6, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r1 = _mm256_add_epi16(r1, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 128), r1);

        montgomery_reduce_257(&r2, &a2, &b0);
        montgomery_reduce_257(&tmp2, &a1, &b1);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b2);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b3);
        montgomery_reduce_257(&tmp2, &a6, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r2 = _mm256_add_epi16(r2, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 16), r2);

        montgomery_reduce_257(&r3, &a3, &b0);
        montgomery_reduce_257(&tmp2, &a2, &b1);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b2);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b3);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b4);
        montgomery_reduce_257(&tmp2, &a6, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r3 = _mm256_add_epi16(r3, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 144), r3);

        montgomery_reduce_257(&r4, &a4, &b0);
        montgomery_reduce_257(&tmp2, &a3, &b1);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b2);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b3);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b4);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b5);
        montgomery_reduce_257(&tmp2, &a6, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r4 = _mm256_add_epi16(r4, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 32), r4);

        montgomery_reduce_257(&r5, &a5, &b0);
        montgomery_reduce_257(&tmp2, &a4, &b1);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b2);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b3);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b4);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b5);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b6);
        montgomery_reduce_257(&tmp2, &a6, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r5 = _mm256_add_epi16(r5, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 160), r5);

        montgomery_reduce_257(&r6, &a6, &b0);
        montgomery_reduce_257(&tmp2, &a5, &b1);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b2);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b3);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b4);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b5);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b6);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b7);
        montgomery_reduce_257(&tmp1, &tmp1, &zeta);
        r6 = _mm256_add_epi16(r6, tmp1);
        _mm256_storeu_si256((__m256i *)(r + i + 48), r6);

        montgomery_reduce_257(&r7, &a7, &b0);
        montgomery_reduce_257(&tmp2, &a6, &b1);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b2);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b3);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b4);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b5);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b6);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b7);
        r7 = _mm256_add_epi16(r7, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 176), r7);


        a0 = _mm256_loadu_si256((__m256i *)(a + i + 64));
        a1 = _mm256_loadu_si256((__m256i *)(a + i + 80));
        a2 = _mm256_loadu_si256((__m256i *)(a + i + 96));
        a3 = _mm256_loadu_si256((__m256i *)(a + i + 112));
        a4 = _mm256_loadu_si256((__m256i *)(a + i + 192));
        a5 = _mm256_loadu_si256((__m256i *)(a + i + 208));
        a6 = _mm256_loadu_si256((__m256i *)(a + i + 224));
        a7 = _mm256_loadu_si256((__m256i *)(a + i + 240));

        b0 = _mm256_loadu_si256((__m256i *)(b + i + 64));
        b1 = _mm256_loadu_si256((__m256i *)(b + i + 80));
        b2 = _mm256_loadu_si256((__m256i *)(b + i + 96));
        b3 = _mm256_loadu_si256((__m256i *)(b + i + 112));
        b4 = _mm256_loadu_si256((__m256i *)(b + i + 192));
        b5 = _mm256_loadu_si256((__m256i *)(b + i + 208));
        b6 = _mm256_loadu_si256((__m256i *)(b + i + 224));
        b7 = _mm256_loadu_si256((__m256i *)(b + i + 240));

        montgomery_reduce_257(&r0, &a0, &b0);
        montgomery_reduce_257(&tmp1, &a7, &b1);
        montgomery_reduce_257(&tmp2, &a6, &b2);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r0 = _mm256_sub_epi16(r0, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 64), r0);

        montgomery_reduce_257(&r1, &a1, &b0);
        montgomery_reduce_257(&tmp2, &a0, &b1);
        r1 = _mm256_add_epi16(r1, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b2);
        montgomery_reduce_257(&tmp2, &a6, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r1 = _mm256_sub_epi16(r1, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 192), r1);

        montgomery_reduce_257(&r2, &a2, &b0);
        montgomery_reduce_257(&tmp2, &a1, &b1);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b2);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b3);
        montgomery_reduce_257(&tmp2, &a6, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r2 = _mm256_sub_epi16(r2, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 80), r2);

        montgomery_reduce_257(&r3, &a3, &b0);
        montgomery_reduce_257(&tmp2, &a2, &b1);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b2);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b3);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b4);
        montgomery_reduce_257(&tmp2, &a6, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r3 = _mm256_sub_epi16(r3, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 208), r3);

        montgomery_reduce_257(&r4, &a4, &b0);
        montgomery_reduce_257(&tmp2, &a3, &b1);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b2);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b3);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b4);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b5);
        montgomery_reduce_257(&tmp2, &a6, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r4 = _mm256_sub_epi16(r4, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 96), r4);

        montgomery_reduce_257(&r5, &a5, &b0);
        montgomery_reduce_257(&tmp2, &a4, &b1);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b2);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b3);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b4);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b5);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b6);
        montgomery_reduce_257(&tmp2, &a6, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        r5 = _mm256_sub_epi16(r5, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 224), r5);

        montgomery_reduce_257(&r6, &a6, &b0);
        montgomery_reduce_257(&tmp2, &a5, &b1);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b2);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b3);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b4);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b5);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b6);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b7);
        montgomery_reduce_257(&tmp1, &tmp1, &zeta);
        r6 = _mm256_sub_epi16(r6, tmp1);
        _mm256_storeu_si256((__m256i *)(r + i + 112), r6);

        montgomery_reduce_257(&r7, &a7, &b0);
        montgomery_reduce_257(&tmp2, &a6, &b1);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b2);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b3);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b4);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b5);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b6);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b7);
        r7 = _mm256_add_epi16(r7, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 240), r7);
    }  
}

static void montgomery_reduce_769(__m256i *r, __m256i *a, __m256i *b)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(769);
    __m256i tmp_QINV = _mm256_set1_epi16(-767);

    tmp_c0 = _mm256_mullo_epi16(*a, *b);
    tmp_c1 = _mm256_mulhi_epi16(*a, *b);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *r     = _mm256_sub_epi16(tmp_c1, tmp_x);
}

void mq_poly_mul_ntt_769(int16_t *r, int16_t *a, int16_t *b, int16_t *muldata)
{
    int i;
    __m256i zeta, tmp1,tmp2;
    __m256i a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7;
    __m256i r0, r1, r2, r3, r4, r5, r6, r7;
    for(i = 0;i < DIM_N;i += 256)
    {
        zeta = _mm256_loadu_si256((__m256i *)(muldata + i / 16));

        a0 = _mm256_loadu_si256((__m256i *)(a + i));
        a1 = _mm256_loadu_si256((__m256i *)(a + i + 16));
        a2 = _mm256_loadu_si256((__m256i *)(a + i + 32));
        a3 = _mm256_loadu_si256((__m256i *)(a + i + 48));
        a4 = _mm256_loadu_si256((__m256i *)(a + i + 128));
        a5 = _mm256_loadu_si256((__m256i *)(a + i + 144));
        a6 = _mm256_loadu_si256((__m256i *)(a + i + 160));
        a7 = _mm256_loadu_si256((__m256i *)(a + i + 176));

        b0 = _mm256_loadu_si256((__m256i *)(b + i));
        b1 = _mm256_loadu_si256((__m256i *)(b + i + 16));
        b2 = _mm256_loadu_si256((__m256i *)(b + i + 32));
        b3 = _mm256_loadu_si256((__m256i *)(b + i + 48));
        b4 = _mm256_loadu_si256((__m256i *)(b + i + 128));
        b5 = _mm256_loadu_si256((__m256i *)(b + i + 144));
        b6 = _mm256_loadu_si256((__m256i *)(b + i + 160));
        b7 = _mm256_loadu_si256((__m256i *)(b + i + 176));

        montgomery_reduce_769(&r0, &a0, &b0);
        montgomery_reduce_769(&tmp1, &a7, &b1);
        montgomery_reduce_769(&tmp2, &a6, &b2);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r0 = _mm256_add_epi16(r0, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i), r0);

        montgomery_reduce_769(&r1, &a1, &b0);
        montgomery_reduce_769(&tmp2, &a0, &b1);
        r1 = _mm256_add_epi16(r1, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b2);
        montgomery_reduce_769(&tmp2, &a6, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r1 = _mm256_add_epi16(r1, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 128), r1);

        montgomery_reduce_769(&r2, &a2, &b0);
        montgomery_reduce_769(&tmp2, &a1, &b1);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b2);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b3);
        montgomery_reduce_769(&tmp2, &a6, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r2 = _mm256_add_epi16(r2, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 16), r2);

        montgomery_reduce_769(&r3, &a3, &b0);
        montgomery_reduce_769(&tmp2, &a2, &b1);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b2);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b3);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b4);
        montgomery_reduce_769(&tmp2, &a6, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r3 = _mm256_add_epi16(r3, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 144), r3);

        montgomery_reduce_769(&r4, &a4, &b0);
        montgomery_reduce_769(&tmp2, &a3, &b1);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b2);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b3);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b4);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b5);
        montgomery_reduce_769(&tmp2, &a6, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r4 = _mm256_add_epi16(r4, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 32), r4);

        montgomery_reduce_769(&r5, &a5, &b0);
        montgomery_reduce_769(&tmp2, &a4, &b1);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b2);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b3);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b4);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b5);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b6);
        montgomery_reduce_769(&tmp2, &a6, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r5 = _mm256_add_epi16(r5, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 160), r5);

        montgomery_reduce_769(&r6, &a6, &b0);
        montgomery_reduce_769(&tmp2, &a5, &b1);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b2);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b3);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b4);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b5);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b6);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b7);
        montgomery_reduce_769(&tmp1, &tmp1, &zeta);
        r6 = _mm256_add_epi16(r6, tmp1);
        _mm256_storeu_si256((__m256i *)(r + i + 48), r6);

        montgomery_reduce_769(&r7, &a7, &b0);
        montgomery_reduce_769(&tmp2, &a6, &b1);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b2);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b3);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b4);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b5);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b6);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b7);
        r7 = _mm256_add_epi16(r7, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 176), r7);


        a0 = _mm256_loadu_si256((__m256i *)(a + i + 64));
        a1 = _mm256_loadu_si256((__m256i *)(a + i + 80));
        a2 = _mm256_loadu_si256((__m256i *)(a + i + 96));
        a3 = _mm256_loadu_si256((__m256i *)(a + i + 112));
        a4 = _mm256_loadu_si256((__m256i *)(a + i + 192));
        a5 = _mm256_loadu_si256((__m256i *)(a + i + 208));
        a6 = _mm256_loadu_si256((__m256i *)(a + i + 224));
        a7 = _mm256_loadu_si256((__m256i *)(a + i + 240));

        b0 = _mm256_loadu_si256((__m256i *)(b + i + 64));
        b1 = _mm256_loadu_si256((__m256i *)(b + i + 80));
        b2 = _mm256_loadu_si256((__m256i *)(b + i + 96));
        b3 = _mm256_loadu_si256((__m256i *)(b + i + 112));
        b4 = _mm256_loadu_si256((__m256i *)(b + i + 192));
        b5 = _mm256_loadu_si256((__m256i *)(b + i + 208));
        b6 = _mm256_loadu_si256((__m256i *)(b + i + 224));
        b7 = _mm256_loadu_si256((__m256i *)(b + i + 240));

        montgomery_reduce_769(&r0, &a0, &b0);
        montgomery_reduce_769(&tmp1, &a7, &b1);
        montgomery_reduce_769(&tmp2, &a6, &b2);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r0 = _mm256_sub_epi16(r0, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 64), r0);

        montgomery_reduce_769(&r1, &a1, &b0);
        montgomery_reduce_769(&tmp2, &a0, &b1);
        r1 = _mm256_add_epi16(r1, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b2);
        montgomery_reduce_769(&tmp2, &a6, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r1 = _mm256_sub_epi16(r1, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 192), r1);

        montgomery_reduce_769(&r2, &a2, &b0);
        montgomery_reduce_769(&tmp2, &a1, &b1);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b2);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b3);
        montgomery_reduce_769(&tmp2, &a6, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r2 = _mm256_sub_epi16(r2, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 80), r2);

        montgomery_reduce_769(&r3, &a3, &b0);
        montgomery_reduce_769(&tmp2, &a2, &b1);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b2);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b3);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b4);
        montgomery_reduce_769(&tmp2, &a6, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r3 = _mm256_sub_epi16(r3, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 208), r3);

        montgomery_reduce_769(&r4, &a4, &b0);
        montgomery_reduce_769(&tmp2, &a3, &b1);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b2);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b3);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b4);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b5);
        montgomery_reduce_769(&tmp2, &a6, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r4 = _mm256_sub_epi16(r4, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 96), r4);

        montgomery_reduce_769(&r5, &a5, &b0);
        montgomery_reduce_769(&tmp2, &a4, &b1);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b2);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b3);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b4);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b5);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b6);
        montgomery_reduce_769(&tmp2, &a6, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_769(&tmp2, &tmp1, &zeta);
        r5 = _mm256_sub_epi16(r5, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 224), r5);

        montgomery_reduce_769(&r6, &a6, &b0);
        montgomery_reduce_769(&tmp2, &a5, &b1);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b2);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b3);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b4);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b5);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b6);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_769(&tmp1, &a7, &b7);
        montgomery_reduce_769(&tmp1, &tmp1, &zeta);
        r6 = _mm256_sub_epi16(r6, tmp1);
        _mm256_storeu_si256((__m256i *)(r + i + 112), r6);

        montgomery_reduce_769(&r7, &a7, &b0);
        montgomery_reduce_769(&tmp2, &a6, &b1);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a5, &b2);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a4, &b3);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a3, &b4);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a2, &b5);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a1, &b6);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_769(&tmp2, &a0, &b7);
        r7 = _mm256_add_epi16(r7, tmp2);
        _mm256_storeu_si256((__m256i *)(r + i + 240), r7);
    }  
}

static void montgomery_mul_769(__m256i *tmp_c,__m256i *tmp_a,__m256i *tmp_b,__m256i *tmp_Q,__m256i *tmp_QINV)
{
    __m256i tmp_c0, tmp_c1, tmp_x;
    tmp_c0 = _mm256_mullo_epi16(*tmp_a, *tmp_b);
    tmp_c1 = _mm256_mulhi_epi16(*tmp_a, *tmp_b);
    tmp_x = _mm256_mullo_epi16(tmp_c0, *tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, *tmp_Q);
    *tmp_c = _mm256_sub_epi16(tmp_c1, tmp_x);

}

void mq_poly_mul_ntt_257_mq(int16_t *r, int16_t *a, int16_t *b, int16_t *muldata)
{
    int i;
    __m256i zeta, tmp1, tmp2;
    __m256i a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7;
    __m256i r0, r1, r2, r3, r4, r5, r6, r7;
    __m256i tmp_Q = _mm256_set1_epi16(257);
    __m256i set1 = _mm256_set1_epi16(1);
    for(i = 0; i < DIM_N; i += 256)
    {
        zeta = _mm256_load_si256((__m256i *)(muldata + i / 16));

        a0 = _mm256_load_si256((__m256i *)(a + i));
        a1 = _mm256_load_si256((__m256i *)(a + i + 16));
        a2 = _mm256_load_si256((__m256i *)(a + i + 32));
        a3 = _mm256_load_si256((__m256i *)(a + i + 48));
        a4 = _mm256_load_si256((__m256i *)(a + i + 128));
        a5 = _mm256_load_si256((__m256i *)(a + i + 144));
        a6 = _mm256_load_si256((__m256i *)(a + i + 160));
        a7 = _mm256_load_si256((__m256i *)(a + i + 176));

        b0 = _mm256_load_si256((__m256i *)(b + i));
        b1 = _mm256_load_si256((__m256i *)(b + i + 16));
        b2 = _mm256_load_si256((__m256i *)(b + i + 32));
        b3 = _mm256_load_si256((__m256i *)(b + i + 48));
        b4 = _mm256_load_si256((__m256i *)(b + i + 128));
        b5 = _mm256_load_si256((__m256i *)(b + i + 144));
        b6 = _mm256_load_si256((__m256i *)(b + i + 160));
        b7 = _mm256_load_si256((__m256i *)(b + i + 176));

        montgomery_reduce_257(&r0, &a0, &b0);
        montgomery_reduce_257(&tmp1, &a7, &b1);
        montgomery_reduce_257(&tmp2, &a6, &b2);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_add_epi16(r0, tmp2);
        montgomery_reduce_257(&r0, &tmp1, &set1);
        r0 = _mm256_add_epi16(r0, _mm256_and_si256(_mm256_srai_epi16(r0, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i), r0);

        montgomery_reduce_257(&r1, &a1, &b0);
        montgomery_reduce_257(&tmp2, &a0, &b1);
        r1 = _mm256_add_epi16(r1, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b2);
        montgomery_reduce_257(&tmp2, &a6, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_add_epi16(r1, tmp2);
        montgomery_reduce_257(&r1, &tmp1, &set1);
        r1 = _mm256_add_epi16(r1, _mm256_and_si256(_mm256_srai_epi16(r1, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 16), r1);

        montgomery_reduce_257(&r2, &a2, &b0);
        montgomery_reduce_257(&tmp2, &a1, &b1);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b2);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b3);
        montgomery_reduce_257(&tmp2, &a6, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&r2, &tmp1, &set1);
        r2 = _mm256_add_epi16(r2, _mm256_and_si256(_mm256_srai_epi16(r2, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 32), r2);

        montgomery_reduce_257(&r3, &a3, &b0);
        montgomery_reduce_257(&tmp2, &a2, &b1);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b2);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b3);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b4);
        montgomery_reduce_257(&tmp2, &a6, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&r3, &tmp1, &set1);
        r3 = _mm256_add_epi16(r3, _mm256_and_si256(_mm256_srai_epi16(r3, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 48), r3);

        montgomery_reduce_257(&r4, &a4, &b0);
        montgomery_reduce_257(&tmp2, &a3, &b1);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b2);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b3);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b4);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b5);
        montgomery_reduce_257(&tmp2, &a6, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&r4, &tmp1, &set1);
        r4 = _mm256_add_epi16(r4, _mm256_and_si256(_mm256_srai_epi16(r4, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 128), r4);

        montgomery_reduce_257(&r5, &a5, &b0);
        montgomery_reduce_257(&tmp2, &a4, &b1);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b2);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b3);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b4);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b5);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b6);
        montgomery_reduce_257(&tmp2, &a6, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&r5, &tmp1, &set1);
        r5 = _mm256_add_epi16(r5, _mm256_and_si256(_mm256_srai_epi16(r5, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 144), r5);

        montgomery_reduce_257(&r6, &a6, &b0);
        montgomery_reduce_257(&tmp2, &a5, &b1);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b2);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b3);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b4);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b5);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b6);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b7);
        montgomery_reduce_257(&tmp1, &tmp1, &zeta);
        tmp1 = _mm256_add_epi16(r6, tmp1);
        montgomery_reduce_257(&r6, &tmp1, &set1);
        r6 = _mm256_add_epi16(r6, _mm256_and_si256(_mm256_srai_epi16(r6, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 160), r6);

        montgomery_reduce_257(&r7, &a7, &b0);
        montgomery_reduce_257(&tmp2, &a6, &b1);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b2);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b3);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b4);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b5);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b6);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b7);
        tmp1 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&r7, &tmp1, &set1);
        r7 = _mm256_add_epi16(r7, _mm256_and_si256(_mm256_srai_epi16(r7, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 176), r7);


        a0 = _mm256_load_si256((__m256i *)(a + i + 64));
        a1 = _mm256_load_si256((__m256i *)(a + i + 80));
        a2 = _mm256_load_si256((__m256i *)(a + i + 96));
        a3 = _mm256_load_si256((__m256i *)(a + i + 112));
        a4 = _mm256_load_si256((__m256i *)(a + i + 192));
        a5 = _mm256_load_si256((__m256i *)(a + i + 208));
        a6 = _mm256_load_si256((__m256i *)(a + i + 224));
        a7 = _mm256_load_si256((__m256i *)(a + i + 240));

        b0 = _mm256_load_si256((__m256i *)(b + i + 64));
        b1 = _mm256_load_si256((__m256i *)(b + i + 80));
        b2 = _mm256_load_si256((__m256i *)(b + i + 96));
        b3 = _mm256_load_si256((__m256i *)(b + i + 112));
        b4 = _mm256_load_si256((__m256i *)(b + i + 192));
        b5 = _mm256_load_si256((__m256i *)(b + i + 208));
        b6 = _mm256_load_si256((__m256i *)(b + i + 224));
        b7 = _mm256_load_si256((__m256i *)(b + i + 240));

        montgomery_reduce_257(&r0, &a0, &b0);
        montgomery_reduce_257(&tmp1, &a7, &b1);
        montgomery_reduce_257(&tmp2, &a6, &b2);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_sub_epi16(r0, tmp2);
        montgomery_reduce_257(&r0, &tmp1, &set1);
        r0 = _mm256_add_epi16(r0, _mm256_and_si256(_mm256_srai_epi16(r0, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 64), r0);

        montgomery_reduce_257(&r1, &a1, &b0);
        montgomery_reduce_257(&tmp2, &a0, &b1);
        r1 = _mm256_add_epi16(r1, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b2);
        montgomery_reduce_257(&tmp2, &a6, &b3);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_sub_epi16(r1, tmp2);
        montgomery_reduce_257(&r1, &tmp1, &set1);
        r1 = _mm256_add_epi16(r1, _mm256_and_si256(_mm256_srai_epi16(r1, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 80), r1);

        montgomery_reduce_257(&r2, &a2, &b0);
        montgomery_reduce_257(&tmp2, &a1, &b1);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b2);
        r2 = _mm256_add_epi16(r2, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b3);
        montgomery_reduce_257(&tmp2, &a6, &b4);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_sub_epi16(r2, tmp2);
        montgomery_reduce_257(&r2, &tmp1, &set1);
        r2 = _mm256_add_epi16(r2, _mm256_and_si256(_mm256_srai_epi16(r2, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 96), r2);

        montgomery_reduce_257(&r3, &a3, &b0);
        montgomery_reduce_257(&tmp2, &a2, &b1);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b2);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b3);
        r3 = _mm256_add_epi16(r3, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b4);
        montgomery_reduce_257(&tmp2, &a6, &b5);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_sub_epi16(r3, tmp2);
        montgomery_reduce_257(&r3, &tmp1, &set1);
        r3 = _mm256_add_epi16(r3, _mm256_and_si256(_mm256_srai_epi16(r3, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 112), r3);

        montgomery_reduce_257(&r4, &a4, &b0);
        montgomery_reduce_257(&tmp2, &a3, &b1);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b2);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b3);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b4);
        r4 = _mm256_add_epi16(r4, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b5);
        montgomery_reduce_257(&tmp2, &a6, &b6);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_sub_epi16(r4, tmp2);
        montgomery_reduce_257(&r4, &tmp1, &set1);
        r4 = _mm256_add_epi16(r4, _mm256_and_si256(_mm256_srai_epi16(r4, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 192), r4);

        montgomery_reduce_257(&r5, &a5, &b0);
        montgomery_reduce_257(&tmp2, &a4, &b1);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b2);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b3);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b4);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b5);
        r5 = _mm256_add_epi16(r5, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b6);
        montgomery_reduce_257(&tmp2, &a6, &b7);
        tmp1 = _mm256_add_epi16(tmp1, tmp2);
        montgomery_reduce_257(&tmp2, &tmp1, &zeta);
        tmp1 = _mm256_sub_epi16(r5, tmp2);
        montgomery_reduce_257(&r5, &tmp1, &set1);
        r5 = _mm256_add_epi16(r5, _mm256_and_si256(_mm256_srai_epi16(r5, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 208), r5);

        montgomery_reduce_257(&r6, &a6, &b0);
        montgomery_reduce_257(&tmp2, &a5, &b1);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b2);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b3);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b4);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b5);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b6);
        r6 = _mm256_add_epi16(r6, tmp2);
        montgomery_reduce_257(&tmp1, &a7, &b7);
        montgomery_reduce_257(&tmp1, &tmp1, &zeta);
        tmp1 = _mm256_sub_epi16(r6, tmp1);
        montgomery_reduce_257(&r6, &tmp1, &set1);
        r6 = _mm256_add_epi16(r6, _mm256_and_si256(_mm256_srai_epi16(r6, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 224), r6);

        montgomery_reduce_257(&r7, &a7, &b0);
        montgomery_reduce_257(&tmp2, &a6, &b1);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a5, &b2);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a4, &b3);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a3, &b4);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a2, &b5);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a1, &b6);
        r7 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&tmp2, &a0, &b7);
        tmp1 = _mm256_add_epi16(r7, tmp2);
        montgomery_reduce_257(&r7, &tmp1, &set1);
        r7 = _mm256_add_epi16(r7, _mm256_and_si256(_mm256_srai_epi16(r7, 15), tmp_Q));
        _mm256_store_si256((__m256i *)(r + i + 240), r7);
    } 
}

int mq_poly_inv_ntt_257(int16_t *r, int16_t *a)
{
    int i, bad = 0;
    __m128i v0, v1;
    __m256i zeta, t, k;
    __m256i a0, a1, a2, a3, a4, a5, a6, a7;
    __m256i b0, b1, b2, b3;
    __m256i c0, c1, e;
    __m256i f0, f1, f2, f3;
    __m256i r0, r1, r2, r3, r4, r5, r6, r7;
    __m256i tmp_Q = _mm256_set1_epi16(257);
    __m256i set0 = _mm256_set1_epi16(0);
    __m256i set2 = _mm256_set1_epi16(2);
    __m256i setn1 = _mm256_set1_epi16(-1);

    int16_t __attribute__((aligned(32))) test[DIM_N];

    for(i = 0; i < DIM_N; i += 256)
    {
        zeta = _mm256_load_si256((__m256i *)(last_level_zeta + i / 8));

        a0 = _mm256_load_si256((__m256i *)(a + i));
        a1 = _mm256_load_si256((__m256i *)(a + i + 16));
        a2 = _mm256_load_si256((__m256i *)(a + i + 32));
        a3 = _mm256_load_si256((__m256i *)(a + i + 48));
        a4 = _mm256_load_si256((__m256i *)(a + i + 128));
        a5 = _mm256_load_si256((__m256i *)(a + i + 144));
        a6 = _mm256_load_si256((__m256i *)(a + i + 160));
        a7 = _mm256_load_si256((__m256i *)(a + i + 176));

        montgomery_reduce_257(&b0, &a1, &a7);
        montgomery_reduce_257(&t, &a3, &a5);
        b0 = _mm256_add_epi16(b0, t);
        montgomery_reduce_257(&t, &a2, &a6);
        b0 = _mm256_sub_epi16(b0, t);
        montgomery_reduce_257(&b0, &b0, &set2);
        montgomery_reduce_257(&t, &a4, &a4);
        b0 = _mm256_sub_epi16(b0, t);
        montgomery_reduce_257(&b0, &b0, &zeta);
        montgomery_reduce_257(&t, &a0, &a0);
        b0 = _mm256_add_epi16(b0, t);

        montgomery_reduce_257(&b1, &a3, &a7);
        montgomery_reduce_257(&t, &a4, &a6);
        b1 = _mm256_sub_epi16(b1, t);
        montgomery_reduce_257(&b1, &b1, &set2);
        montgomery_reduce_257(&t, &a5, &a5);
        b1 = _mm256_add_epi16(b1, t);
        montgomery_reduce_257(&b1, &b1, &zeta);
        montgomery_reduce_257(&t, &a0, &a2);
        montgomery_reduce_257(&t, &t, &set2);
        b1 = _mm256_add_epi16(b1, t);
        montgomery_reduce_257(&t, &a1, &a1);
        b1 = _mm256_sub_epi16(b1, t);

        montgomery_reduce_257(&b2, &a5, &a7);
        montgomery_reduce_257(&b2, &b2, &set2);
        montgomery_reduce_257(&t, &a6, &a6);
        b2 = _mm256_sub_epi16(b2, t);
        montgomery_reduce_257(&b2, &b2, &zeta);
        montgomery_reduce_257(&t, &a0, &a4);
        montgomery_reduce_257(&k, &a1, &a3);
        t = _mm256_sub_epi16(t, k);
        montgomery_reduce_257(&t, &t, &set2);
        b2 = _mm256_add_epi16(b2, t);
        montgomery_reduce_257(&t, &a2, &a2);
        b2 = _mm256_add_epi16(b2, t);

        montgomery_reduce_257(&b3, &a7, &a7);
        montgomery_reduce_257(&b3, &b3, &zeta);
        montgomery_reduce_257(&t, &a0, &a6);
        montgomery_reduce_257(&k, &a2, &a4);
        t = _mm256_add_epi16(t, k);
        montgomery_reduce_257(&k, &a1, &a5);
        t = _mm256_sub_epi16(t, k);
        montgomery_reduce_257(&t, &t, &set2);
        b3 = _mm256_add_epi16(b3, t);
        montgomery_reduce_257(&t, &a3, &a3);
        b3 = _mm256_sub_epi16(b3, t);

        montgomery_reduce_257(&c0, &b1, &b3);
        montgomery_reduce_257(&c0, &c0, &set2);
        montgomery_reduce_257(&t, &b2, &b2);
        c0 = _mm256_sub_epi16(c0, t);
        montgomery_reduce_257(&c0, &c0, &zeta);
        montgomery_reduce_257(&t, &b0, &b0);
        c0 = _mm256_add_epi16(c0, t);

        montgomery_reduce_257(&c1, &b3, &b3);
        montgomery_reduce_257(&c1, &c1, &zeta);
        montgomery_reduce_257(&t, &b0, &b2);
        montgomery_reduce_257(&t, &t, &set2);
        c1 = _mm256_add_epi16(c1, t);
        montgomery_reduce_257(&t, &b1, &b1);
        c1 = _mm256_sub_epi16(c1, t);

        montgomery_reduce_257(&e, &c1, &c1);
        montgomery_reduce_257(&e, &e, &zeta);
        montgomery_reduce_257(&t, &c0, &c0);
        e = _mm256_add_epi16(e, t);
        e = _mm256_add_epi16(e, _mm256_and_si256(_mm256_srai_epi16(e, 15), tmp_Q));
        bad += (_mm256_movemask_epi8(_mm256_cmpeq_epi16(e, set0)) == 0xFFFFFFFF);
        v0 = _mm256_extracti128_si256(e, 0);
        t=_mm256_cvtepi16_epi32(v0);
        t=_mm256_i32gather_epi32(qinv, t, sizeof(int32_t));
        v1=_mm256_extracti128_si256(e, 1);
        k=_mm256_cvtepi16_epi32(v1);
        k=_mm256_i32gather_epi32(qinv, k, sizeof(int32_t));
        t=_mm256_packs_epi32(t, k);
        e=_mm256_permute4x64_epi64(t, 0xd8);

        montgomery_reduce_257(&c0, &e, &c0);

        montgomery_reduce_257(&c1, &e, &c1);
        montgomery_reduce_257(&c1, &setn1, &c1);

        montgomery_reduce_257(&f0, &c1, &b2);
        montgomery_reduce_257(&f0, &f0, &zeta);
        montgomery_reduce_257(&t, &c0, &b0);
        f0 = _mm256_sub_epi16(t, f0);

        montgomery_reduce_257(&f1, &c1, &b3);
        montgomery_reduce_257(&f1, &f1, &zeta);
        montgomery_reduce_257(&t, &c0, &b1);
        f1 = _mm256_sub_epi16(f1, t);

        montgomery_reduce_257(&f2, &c0, &b2);
        montgomery_reduce_257(&t, &c1, &b0);
        f2 = _mm256_add_epi16(f2, t);

        montgomery_reduce_257(&f3, &c0, &b3);
        montgomery_reduce_257(&t, &c1, &b1);
        f3 = _mm256_add_epi16(f3, t);
        montgomery_reduce_257(&f3, &setn1, &f3);

        montgomery_reduce_257(&r0, &f1, &a6);
        montgomery_reduce_257(&t, &f2, &a4);
        r0 = _mm256_add_epi16(r0, t);
        montgomery_reduce_257(&t, &f3, &a2);
        r0 = _mm256_add_epi16(r0, t);
        montgomery_reduce_257(&r0, &r0, &zeta);
        montgomery_reduce_257(&t, &f0, &a0);
        r0 = _mm256_sub_epi16(t, r0);

        montgomery_reduce_257(&r1, &f1, &a7);
        montgomery_reduce_257(&t, &f2, &a5);
        r1 = _mm256_add_epi16(r1, t);
        montgomery_reduce_257(&t, &f3, &a3);
        r1 = _mm256_add_epi16(r1, t);
        montgomery_reduce_257(&r1, &r1, &zeta);
        montgomery_reduce_257(&t, &f0, &a1);
        r1 = _mm256_sub_epi16(r1, t);

        montgomery_reduce_257(&r2, &f2, &a6);
        montgomery_reduce_257(&t, &f3, &a4);
        r2 = _mm256_add_epi16(r2, t);
        montgomery_reduce_257(&r2, &r2, &zeta);
        montgomery_reduce_257(&t, &f0, &a2);
        montgomery_reduce_257(&k, &f1, &a0);
        t = _mm256_add_epi16(t, k);
        r2 = _mm256_sub_epi16(t, r2);

        montgomery_reduce_257(&r3, &f2, &a7);
        montgomery_reduce_257(&t, &f3, &a5);
        r3 = _mm256_add_epi16(r3, t);
        montgomery_reduce_257(&r3, &r3, &zeta);
        montgomery_reduce_257(&t, &f0, &a3);
        montgomery_reduce_257(&k, &f1, &a1);
        t = _mm256_add_epi16(t, k);
        r3 = _mm256_sub_epi16(r3, t);

        montgomery_reduce_257(&r4, &f3, &a6);
        montgomery_reduce_257(&r4, &r4, &zeta);
        montgomery_reduce_257(&t, &f0, &a4);
        montgomery_reduce_257(&k, &f1, &a2);
        t = _mm256_add_epi16(t, k);
        montgomery_reduce_257(&k, &f2, &a0);
        t = _mm256_add_epi16(t, k);
        r4 = _mm256_sub_epi16(t, r4);

        montgomery_reduce_257(&r5, &f3, &a7);
        montgomery_reduce_257(&r5, &r5, &zeta);
        montgomery_reduce_257(&t, &f0, &a5);
        montgomery_reduce_257(&k, &f1, &a3);
        t = _mm256_add_epi16(t, k);
        montgomery_reduce_257(&k, &f2, &a1);
        t = _mm256_add_epi16(t, k);
        r5 = _mm256_sub_epi16(r5, t);

        montgomery_reduce_257(&r6, &f1, &a4);
        montgomery_reduce_257(&t, &f2, &a2);
        r6 = _mm256_add_epi16(r6, t);
        montgomery_reduce_257(&t, &f3, &a0);
        r6 = _mm256_add_epi16(r6, t);
        montgomery_reduce_257(&t, &f0, &a6);
        r6 = _mm256_add_epi16(r6, t);

        montgomery_reduce_257(&r7, &f0, &a7);
        montgomery_reduce_257(&t, &f1, &a5);
        r7 = _mm256_add_epi16(r7, t);
        montgomery_reduce_257(&t, &f2, &a3);
        r7 = _mm256_add_epi16(r7, t);
        montgomery_reduce_257(&t, &f3, &a1);
        r7 = _mm256_add_epi16(r7, t);
        montgomery_reduce_257(&r7, &r7, &setn1);

        _mm256_store_si256((__m256i *)(r + i), r0);
        _mm256_store_si256((__m256i *)(r + i + 16), r1);
        _mm256_store_si256((__m256i *)(r + i + 32), r2);
        _mm256_store_si256((__m256i *)(r + i + 48), r3);
        _mm256_store_si256((__m256i *)(r + i + 128), r4);
        _mm256_store_si256((__m256i *)(r + i + 144), r5);
        _mm256_store_si256((__m256i *)(r + i + 160), r6);
        _mm256_store_si256((__m256i *)(r + i + 176), r7);
        
        zeta = _mm256_load_si256((__m256i *)(last_level_zeta + i / 8 + 16));

        a0 = _mm256_load_si256((__m256i *)(a + i + 64));
        a1 = _mm256_load_si256((__m256i *)(a + i + 80));
        a2 = _mm256_load_si256((__m256i *)(a + i + 96));
        a3 = _mm256_load_si256((__m256i *)(a + i + 112));
        a4 = _mm256_load_si256((__m256i *)(a + i + 192));
        a5 = _mm256_load_si256((__m256i *)(a + i + 208));
        a6 = _mm256_load_si256((__m256i *)(a + i + 224));
        a7 = _mm256_load_si256((__m256i *)(a + i + 240));

        montgomery_reduce_257(&b0, &a1, &a7);
        montgomery_reduce_257(&t, &a3, &a5);
        b0 = _mm256_add_epi16(b0, t);
        montgomery_reduce_257(&t, &a2, &a6);
        b0 = _mm256_sub_epi16(b0, t);
        montgomery_reduce_257(&b0, &b0, &set2);
        montgomery_reduce_257(&t, &a4, &a4);
        b0 = _mm256_sub_epi16(b0, t);
        montgomery_reduce_257(&b0, &b0, &zeta);
        montgomery_reduce_257(&t, &a0, &a0);
        b0 = _mm256_add_epi16(b0, t);

        montgomery_reduce_257(&b1, &a3, &a7);
        montgomery_reduce_257(&t, &a4, &a6);
        b1 = _mm256_sub_epi16(b1, t);
        montgomery_reduce_257(&b1, &b1, &set2);
        montgomery_reduce_257(&t, &a5, &a5);
        b1 = _mm256_add_epi16(b1, t);
        montgomery_reduce_257(&b1, &b1, &zeta);
        montgomery_reduce_257(&t, &a0, &a2);
        montgomery_reduce_257(&t, &t, &set2);
        b1 = _mm256_add_epi16(b1, t);
        montgomery_reduce_257(&t, &a1, &a1);
        b1 = _mm256_sub_epi16(b1, t);

        montgomery_reduce_257(&b2, &a5, &a7);
        montgomery_reduce_257(&b2, &b2, &set2);
        montgomery_reduce_257(&t, &a6, &a6);
        b2 = _mm256_sub_epi16(b2, t);
        montgomery_reduce_257(&b2, &b2, &zeta);
        montgomery_reduce_257(&t, &a0, &a4);
        montgomery_reduce_257(&k, &a1, &a3);
        t = _mm256_sub_epi16(t, k);
        montgomery_reduce_257(&t, &t, &set2);
        b2 = _mm256_add_epi16(b2, t);
        montgomery_reduce_257(&t, &a2, &a2);
        b2 = _mm256_add_epi16(b2, t);

        montgomery_reduce_257(&b3, &a7, &a7);
        montgomery_reduce_257(&b3, &b3, &zeta);
        montgomery_reduce_257(&t, &a0, &a6);
        montgomery_reduce_257(&k, &a2, &a4);
        t = _mm256_add_epi16(t, k);
        montgomery_reduce_257(&k, &a1, &a5);
        t = _mm256_sub_epi16(t, k);
        montgomery_reduce_257(&t, &t, &set2);
        b3 = _mm256_add_epi16(b3, t);
        montgomery_reduce_257(&t, &a3, &a3);
        b3 = _mm256_sub_epi16(b3, t);

        montgomery_reduce_257(&c0, &b1, &b3);
        montgomery_reduce_257(&c0, &c0, &set2);
        montgomery_reduce_257(&t, &b2, &b2);
        c0 = _mm256_sub_epi16(c0, t);
        montgomery_reduce_257(&c0, &c0, &zeta);
        montgomery_reduce_257(&t, &b0, &b0);
        c0 = _mm256_add_epi16(c0, t);

        montgomery_reduce_257(&c1, &b3, &b3);
        montgomery_reduce_257(&c1, &c1, &zeta);
        montgomery_reduce_257(&t, &b0, &b2);
        montgomery_reduce_257(&t, &t, &set2);
        c1 = _mm256_add_epi16(c1, t);
        montgomery_reduce_257(&t, &b1, &b1);
        c1 = _mm256_sub_epi16(c1, t);

        montgomery_reduce_257(&e, &c1, &c1);
        montgomery_reduce_257(&e, &e, &zeta);
        montgomery_reduce_257(&t, &c0, &c0);
        e = _mm256_add_epi16(e, t);
        e = _mm256_add_epi16(e, _mm256_and_si256(_mm256_srai_epi16(e, 15), tmp_Q));
        v0 = _mm256_extracti128_si256(e, 0);
        t=_mm256_cvtepi16_epi32(v0);
        t=_mm256_i32gather_epi32(qinv, t, sizeof(int32_t));
        v1=_mm256_extracti128_si256(e, 1);
        k=_mm256_cvtepi16_epi32(v1);
        k=_mm256_i32gather_epi32(qinv, k, sizeof(int32_t));
        t=_mm256_packs_epi32(t, k);
        e=_mm256_permute4x64_epi64(t, 0xd8);

        montgomery_reduce_257(&c0, &e, &c0);

        montgomery_reduce_257(&c1, &e, &c1);
        montgomery_reduce_257(&c1, &setn1, &c1);

        montgomery_reduce_257(&f0, &c1, &b2);
        montgomery_reduce_257(&f0, &f0, &zeta);
        montgomery_reduce_257(&t, &c0, &b0);
        f0 = _mm256_sub_epi16(t, f0);

        montgomery_reduce_257(&f1, &c1, &b3);
        montgomery_reduce_257(&f1, &f1, &zeta);
        montgomery_reduce_257(&t, &c0, &b1);
        f1 = _mm256_sub_epi16(f1, t);

        montgomery_reduce_257(&f2, &c0, &b2);
        montgomery_reduce_257(&t, &c1, &b0);
        f2 = _mm256_add_epi16(f2, t);

        montgomery_reduce_257(&f3, &c0, &b3);
        montgomery_reduce_257(&t, &c1, &b1);
        f3 = _mm256_add_epi16(f3, t);
        montgomery_reduce_257(&f3, &setn1, &f3);

        montgomery_reduce_257(&f3, &c0, &b3);
        montgomery_reduce_257(&t, &c1, &b1);
        f3 = _mm256_add_epi16(f3, t);
        montgomery_reduce_257(&f3, &setn1, &f3);

        montgomery_reduce_257(&r0, &f1, &a6);
        montgomery_reduce_257(&t, &f2, &a4);
        r0 = _mm256_add_epi16(r0, t);
        montgomery_reduce_257(&t, &f3, &a2);
        r0 = _mm256_add_epi16(r0, t);
        montgomery_reduce_257(&r0, &r0, &zeta);
        montgomery_reduce_257(&t, &f0, &a0);
        r0 = _mm256_sub_epi16(t, r0);

        montgomery_reduce_257(&r1, &f1, &a7);
        montgomery_reduce_257(&t, &f2, &a5);
        r1 = _mm256_add_epi16(r1, t);
        montgomery_reduce_257(&t, &f3, &a3);
        r1 = _mm256_add_epi16(r1, t);
        montgomery_reduce_257(&r1, &r1, &zeta);
        montgomery_reduce_257(&t, &f0, &a1);
        r1 = _mm256_sub_epi16(r1, t);

        montgomery_reduce_257(&r2, &f2, &a6);
        montgomery_reduce_257(&t, &f3, &a4);
        r2 = _mm256_add_epi16(r2, t);
        montgomery_reduce_257(&r2, &r2, &zeta);
        montgomery_reduce_257(&t, &f0, &a2);
        montgomery_reduce_257(&k, &f1, &a0);
        t = _mm256_add_epi16(t, k);
        r2 = _mm256_sub_epi16(t, r2);

        montgomery_reduce_257(&r3, &f2, &a7);
        montgomery_reduce_257(&t, &f3, &a5);
        r3 = _mm256_add_epi16(r3, t);
        montgomery_reduce_257(&r3, &r3, &zeta);
        montgomery_reduce_257(&t, &f0, &a3);
        montgomery_reduce_257(&k, &f1, &a1);
        t = _mm256_add_epi16(t, k);
        r3 = _mm256_sub_epi16(r3, t);

        montgomery_reduce_257(&r4, &f3, &a6);
        montgomery_reduce_257(&r4, &r4, &zeta);
        montgomery_reduce_257(&t, &f0, &a4);
        montgomery_reduce_257(&k, &f1, &a2);
        t = _mm256_add_epi16(t, k);
        montgomery_reduce_257(&k, &f2, &a0);
        t = _mm256_add_epi16(t, k);
        r4 = _mm256_sub_epi16(t, r4);

        montgomery_reduce_257(&r5, &f3, &a7);
        montgomery_reduce_257(&r5, &r5, &zeta);
        montgomery_reduce_257(&t, &f0, &a5);
        montgomery_reduce_257(&k, &f1, &a3);
        t = _mm256_add_epi16(t, k);
        montgomery_reduce_257(&k, &f2, &a1);
        t = _mm256_add_epi16(t, k);
        r5 = _mm256_sub_epi16(r5, t);

        montgomery_reduce_257(&r6, &f1, &a4);
        montgomery_reduce_257(&t, &f2, &a2);
        r6 = _mm256_add_epi16(r6, t);
        montgomery_reduce_257(&t, &f3, &a0);
        r6 = _mm256_add_epi16(r6, t);
        montgomery_reduce_257(&t, &f0, &a6);
        r6 = _mm256_add_epi16(r6, t);

        montgomery_reduce_257(&r7, &f0, &a7);
        montgomery_reduce_257(&t, &f1, &a5);
        r7 = _mm256_add_epi16(r7, t);
        montgomery_reduce_257(&t, &f2, &a3);
        r7 = _mm256_add_epi16(r7, t);
        montgomery_reduce_257(&t, &f3, &a1);
        r7 = _mm256_add_epi16(r7, t);
        montgomery_reduce_257(&r7, &r7, &setn1);

        _mm256_store_si256((__m256i *)(r + i + 64), r0);
        _mm256_store_si256((__m256i *)(r + i + 80), r1);
        _mm256_store_si256((__m256i *)(r + i + 96), r2);
        _mm256_store_si256((__m256i *)(r + i + 112), r3);
        _mm256_store_si256((__m256i *)(r + i + 192), r4);
        _mm256_store_si256((__m256i *)(r + i + 208), r5);
        _mm256_store_si256((__m256i *)(r + i + 224), r6);
        _mm256_store_si256((__m256i *)(r + i + 240), r7);

    }
    return bad;
}


static void montgomery_reduce_pre_769(__m256i *r, __m256i *a, __m256i *b, __m256i *b_qinv)
{
    __m256i tmp_c0, tmp_c1;

    __m256i tmp_Q = _mm256_set1_epi16(769);

    tmp_c0 = _mm256_mullo_epi16(*a, *b_qinv);
    tmp_c1 = _mm256_mulhi_epi16(*a, *b);

    tmp_c0  = _mm256_mulhi_epi16(tmp_c0, tmp_Q);
    *r     = _mm256_sub_epi16(tmp_c1, tmp_c0);
}

static void montgomery_reduce_pre_257(__m256i *r, __m256i *a, __m256i *b, __m256i *b_qinv)
{
    __m256i tmp_c0, tmp_c1;

    __m256i tmp_Q = _mm256_set1_epi16(257);

    tmp_c0 = _mm256_mullo_epi16(*a, *b_qinv);
    tmp_c1 = _mm256_mulhi_epi16(*a, *b);

    tmp_c0  = _mm256_mulhi_epi16(tmp_c0, tmp_Q);
    *r     = _mm256_sub_epi16(tmp_c1, tmp_c0);
}


void mq_poly_adjoint_ntt_769(int16_t *r, int16_t *a)
{
    int i;
    __m256i tmp0, tmp1, tmp2, tmp3;
    __m256i tmp_zeta1, tmp_zeta_qinv1,tmp_zeta2, tmp_zeta_qinv2;;
    for(i=0;i<DIM_N;i+=256)
    {
        tmp_zeta_qinv1 = _mm256_loadu_si256((__m256i *)(adjoint_769 + i/4));
        tmp_zeta1 = _mm256_loadu_si256((__m256i *)(adjoint_769 + i/4 + 16));

        tmp_zeta_qinv2 = _mm256_loadu_si256((__m256i *)(adjoint_769 + i/4 + 32));
        tmp_zeta2 = _mm256_loadu_si256((__m256i *)(adjoint_769 + i/4 + 48));

        tmp0 = _mm256_loadu_si256((__m256i *)(a + i));
        // montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 192), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 176 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 176), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 160 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 160), tmp3);
        
        tmp0 = _mm256_loadu_si256((__m256i *)(a + 144 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 144), tmp3);

        

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 128 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 64), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 48 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 48), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 32 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 32), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 16 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 16), tmp3);





        tmp0 = _mm256_loadu_si256((__m256i *)(a + 64 + i));
        // montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 256), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 240 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 240), tmp3);
        
        tmp0 = _mm256_loadu_si256((__m256i *)(a + 224 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 224), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 208 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 208), tmp3);



        tmp0 = _mm256_loadu_si256((__m256i *)(a + 192 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 128), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 112 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 112), tmp3);
        
        tmp0 = _mm256_loadu_si256((__m256i *)(a + 96 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 96), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 80 + i));
        montgomery_reduce_pre_769(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 80), tmp3);
    }
    
}

void mq_poly_adjoint_ntt_257(int16_t *r, int16_t *a)
{
    int i;
    __m256i tmp0, tmp1, tmp2, tmp3;
    __m256i tmp_zeta1, tmp_zeta_qinv1,tmp_zeta2, tmp_zeta_qinv2;;
    for(i=0;i<DIM_N;i+=256)
    {
        tmp_zeta_qinv1 = _mm256_loadu_si256((__m256i *)(adjoint_257 + i/4));
        tmp_zeta1 = _mm256_loadu_si256((__m256i *)(adjoint_257 + i/4 + 16));

        tmp_zeta_qinv2 = _mm256_loadu_si256((__m256i *)(adjoint_257 + i/4 + 32));
        tmp_zeta2 = _mm256_loadu_si256((__m256i *)(adjoint_257 + i/4 + 48));

        tmp0 = _mm256_loadu_si256((__m256i *)(a + i));
        // montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 192), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 176 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 176), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 160 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 160), tmp3);
        
        tmp0 = _mm256_loadu_si256((__m256i *)(a + 144 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 144), tmp3);

        

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 128 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 64), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 48 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 48), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 32 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 32), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 16 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta1, &tmp_zeta_qinv1);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 16), tmp3);





        tmp0 = _mm256_loadu_si256((__m256i *)(a + 64 + i));
        // montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 256), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 240 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 240), tmp3);
        
        tmp0 = _mm256_loadu_si256((__m256i *)(a + 224 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 224), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 208 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 208), tmp3);



        tmp0 = _mm256_loadu_si256((__m256i *)(a + 192 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 128), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 112 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 112), tmp3);
        
        tmp0 = _mm256_loadu_si256((__m256i *)(a + 96 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 96), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a + 80 + i));
        montgomery_reduce_pre_257(&tmp0, &tmp0, &tmp_zeta2, &tmp_zeta_qinv2);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 80), tmp3);
    }
    
}

int mq_poly_mul_ntt_12289(int16_t *c, int16_t *a, int16_t *b)
{
    int i;

    __m256i tmp_a, tmp_b, tmp_c1, tmp_c0, tmp_x;
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);
    __m256i tmp_Q = _mm256_set1_epi16(12289);

    for (i = 0; i < DIM_N; i += 16)
    {
        //load
        tmp_a = _mm256_load_si256((__m256i *)(a + i));
        tmp_b = _mm256_load_si256((__m256i *)(b + i));
        //c=a*b
        tmp_c0 = _mm256_mullo_epi16(tmp_a, tmp_b);
        tmp_c1 = _mm256_mulhi_epi16(tmp_a, tmp_b);
        //montgomery mod 12289
        tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
        tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
        tmp_c1 = _mm256_sub_epi16(tmp_c1, tmp_x);
        //if c1<0 then +12289
        tmp_c1 = _mm256_add_epi16(tmp_c1, _mm256_and_si256(_mm256_srai_epi16(tmp_c1, 15), tmp_Q));
        //store
        _mm256_store_si256((__m256i *)(c + i), tmp_c1);
    }

    return 0;
}

void mq_poly_adjoint_ntt_12289(int16_t *c, int16_t *a)
{
    int i;
    __m256i tmp0, tmp1, tmp2, tmp3;
    for(i=0;i<DIM_N;i+=16)
    {
        tmp0 = _mm256_loadu_si256((__m256i *)(a + i));
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(c + DIM_N - i - 16), tmp3);
    }
    
}


// V=V*S%12289, A1=(U+V)%12289, A2=(U-V)%12289
static int ntt_core(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287); //12289^{-1} mod 2^16
    //V*S
    tmp_c0 = _mm256_mullo_epi16(*V, *S);
    tmp_c1 = _mm256_mulhi_epi16(*V, *S);
    //V*S montgomery mod 12289
    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *V     = _mm256_sub_epi16(tmp_c1, tmp_x);
    //if c1<0 then +12289
    *V     = _mm256_add_epi16(*V, _mm256_and_si256(_mm256_srai_epi16(*V, 15), tmp_Q));
    //U%12289 
    *U     = _mm256_add_epi16(*U, _mm256_and_si256(_mm256_srai_epi16(*U, 15), tmp_Q));
    // (U mod 12289) +V -12289
    *A1    = _mm256_add_epi16(*U, *V);
    *A1    = _mm256_sub_epi16(*A1, tmp_Q);
    //(U mod 12289)-V
    *A2   = _mm256_sub_epi16(*U, *V);

    return 0;
}

static int ntt_core_lazymod(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;
    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287); //12289^{-1} mod 2^16
    //V*S
    tmp_c0 = _mm256_mullo_epi16(*V, *S);
    tmp_c1 = _mm256_mulhi_epi16(*V, *S);
    //V*S montgomery mod 12289
    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *V     = _mm256_sub_epi16(tmp_c1, tmp_x);
    //U+V
    *A1    = _mm256_add_epi16(*U, *V);
    //U-V
    *A2   = _mm256_sub_epi16(*U, *V);

    return 0;
}

//  A1=(U+V)%12289, A2=(U-V)*S%12289
static int intt_core(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287); //12289^{-1} mod 2^16
    //A1=(U+V)%12289
    *A1 = _mm256_add_epi16(*U, *V);
    *A1 = _mm256_sub_epi16(*A1, tmp_Q);
    *A1 = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));
    //A2=V=(U-V)%12289 
    *A2 = _mm256_sub_epi16(*U, *V);
    //*A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));
    //V*S
    tmp_c0 = _mm256_mullo_epi16(*A2, *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2, *S);
    //mod 12289
    //montgomery mod 12289
    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2    = _mm256_sub_epi16(tmp_c1, tmp_x);
    //if c1<0 then +12289
    *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}

//level n-1 of intt, lazymod
static int intt_core_lazymod(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287); //12289^{-1} mod 2^16
    //A1=(U+V)%12289
    *A1 = _mm256_add_epi16(*U, *V);
    *A1 = _mm256_sub_epi16(*A1, tmp_Q);
    //V=(U-V)
    *A2 = _mm256_sub_epi16(*U, *V);
    //V*S
    tmp_c0 = _mm256_mullo_epi16(*A2, *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2, *S);
    //mod 12289
    //montgomery mod 12289
    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2    = _mm256_sub_epi16(tmp_c1, tmp_x);

    return 0;
}

static int intt_core_final(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287); //12289^{-1} mod 2^16
    __m256i tmp_b = _mm256_set1_epi16(3755);//N^-1*2^16*2*16 mod 12289
     //A1=(U+V)%12289
    *A1 = _mm256_add_epi16(*U, *V);
    //c=A1*b
    tmp_c0 = _mm256_mullo_epi16(*A1, tmp_b);
    tmp_c1 = _mm256_mulhi_epi16(*A1, tmp_b);
    //montgomery mod 12289
    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A1   = _mm256_sub_epi16(tmp_c1, tmp_x);
    //if c1<0 then +12289
    // *A1 = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));
    //A2=(U-V)%12289 
    *A2 = _mm256_sub_epi16(*U, *V);
    //*A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));
     //c=A2*S
    tmp_c0 = _mm256_mullo_epi16(*A2,  *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2,  *S);
    //montgomery mod 12289
    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2   = _mm256_sub_epi16(tmp_c1, tmp_x);
    //if A2<0 then +12289
    // *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}

static int intt_core_final_n(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287); //12289^{-1} mod 2^16
    __m256i tmp_b = _mm256_set1_epi16(3755);//N^-1*2^16*2*16 mod 12289
     //A1=(U+V)%12289
    *A1 = _mm256_add_epi16(*U, *V);
    //c=A1*b
    tmp_c0 = _mm256_mullo_epi16(*A1, tmp_b);
    tmp_c1 = _mm256_mulhi_epi16(*A1, tmp_b);
    //montgomery mod 12289
    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A1   = _mm256_sub_epi16(tmp_c1, tmp_x);
    //if c1<0 then +12289
    *A1 = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));
    //A2=(U-V)%12289 
    *A2 = _mm256_sub_epi16(*U, *V);
    *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));
     //c=A2*S
    tmp_c0 = _mm256_mullo_epi16(*A2,  *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2,  *S);
    //montgomery mod 12289
    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2   = _mm256_sub_epi16(tmp_c1, tmp_x);
    //if A2<0 then +12289
    *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}

int mq_poly_ntt_12289(int16_t *a)
{
    int i, j, s;

    __m256i tmp_U, tmp_V, tmp_A1, tmp_A2, tmp_S;

    __m256i tmp_shuffle16_pre = _mm256_set_epi8(0x0f, 0x0e, 0x0b, 0x0a, 0x07, 0x06, 0x03, 0x02, 0x0d, 0x0c, 0x09, 0x08, 0x05, 0x04, 0x01, 0x00,
                                        0x0f, 0x0e, 0x0b, 0x0a, 0x07, 0x06, 0x03, 0x02, 0x0d, 0x0c, 0x09, 0x08, 0x05, 0x04, 0x01, 0x00);

    //level 1
    tmp_S = _mm256_set1_epi16(f_12289[1]);
    for(j=0;j<512;j+=16)
    {
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a +j));
        tmp_V = _mm256_load_si256((__m256i *)(a  + 512+j));
        // V*S%12289, (U+V)%12289, (U-V)%12289
        ntt_core_lazymod(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a +j), tmp_A1);
        _mm256_store_si256((__m256i *)(a  + 512+j), tmp_A2);
    }
    
    //level 2
    s=0;
    for (i = 0; i < DIM_N; i += 512)
    {
        //load root
        tmp_S = _mm256_set1_epi16(f_12289[2 +(s++)]);
        for(j=0;j<256;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 256+j));
            // V*S%12289, (U+V)%12289, (U-V)%12289
            ntt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 256+j), tmp_A2);
        }
    }

    //level 3
    s=0;
    for (i = 0; i < DIM_N; i += 256)
    {
        //load root
        tmp_S = _mm256_set1_epi16(f_12289[4 +(s++)]);
        for(j=0;j<128;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 128+j));
            // V*S%12289, (U+V)%12289, (U-V)%12289
            ntt_core_lazymod(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 128+j), tmp_A2);
        }
    }

    //level 4
    s=0;
    for (i = 0; i < DIM_N; i += 128)
    {
        //load root
        tmp_S = _mm256_set1_epi16(f_12289[8 +(s++)]);
        for(j=0;j<64;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 64+j));
            // V*S%12289, (U+V)%12289, (U-V)%12289
            ntt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 64+j), tmp_A2);
        }
    }

    //level 5
    s=0;
    for (i = 0; i < DIM_N; i += 64)
    {
        //load root
        tmp_S = _mm256_set1_epi16(f_12289[16 +(s++)]);
        //load block 1
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 32));
        // V*S%12289, (U+V)%12289, (U-V)%12289
        ntt_core_lazymod(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 32), tmp_A2);

         //load block 2
        tmp_U = _mm256_load_si256((__m256i *)(a + i+16));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 48));
        // V*S%12289, (U+V)%12289, (U-V)%12289
        ntt_core_lazymod(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i+16), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 48), tmp_A2);
    }

    //level 6
    s=0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //load
        tmp_S = _mm256_set1_epi16(f_12289[32 +(s++)]);
        //j=0
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
        // V*S%12289, (U+V)%12289, (U-V)%12289
        ntt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }
    
    
    //level 7
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f_12289[m+i];
        tmp_S = _mm256_set_epi16(f_12289[64 + s + 1], f_12289[64 + s + 1], f_12289[64 + s + 1], f_12289[64 + s + 1], f_12289[64 + s + 1], f_12289[64 + s + 1], f_12289[64 + s + 1], f_12289[64 + s + 1],
                                 f_12289[64 + s], f_12289[64 + s], f_12289[64 + s], f_12289[64 + s], f_12289[64 + s], f_12289[64 + s], f_12289[64 + s], f_12289[64 + s]);
        s += 2;
        //load
        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));
        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        //V*S%12289, (U+V)%12289, (U-V)%12289
        ntt_core_lazymod(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }


    //level 8
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f_12289[m+i];
        tmp_S = _mm256_set_epi16(f_12289[128 + s + 3], f_12289[128 + s + 3], f_12289[128 + s + 3], f_12289[128 + s + 3], f_12289[128 + s + 1], f_12289[128 + s + 1], f_12289[128 + s + 1], f_12289[128 + s + 1],
                                 f_12289[128 + s + 2], f_12289[128 + s + 2], f_12289[128 + s + 2], f_12289[128 + s + 2], f_12289[128 + s], f_12289[128 + s], f_12289[128 + s], f_12289[128 + s]);
        s += 4;
        //load
        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));
        //permute
        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        //V*S%12289, (U+V)%12289, (U-V)%12289
        ntt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    //level 9
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f_12289[m+i];
        tmp_S = _mm256_set_epi16(f_12289[256 + s + 7], f_12289[256 + s + 7], f_12289[256 + s + 3], f_12289[256 + s + 3], f_12289[256 + s + 5], f_12289[256 + s + 5], f_12289[256 + s + 1], f_12289[256 + s + 1],
                                 f_12289[256 + s + 6], f_12289[256 + s + 6], f_12289[256 + s + 2], f_12289[256 + s + 2], f_12289[256 + s + 4], f_12289[256 + s + 4], f_12289[256 + s], f_12289[256 + s]);
        s += 8;
        //load
        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));
        //permute
        tmp_A1 = _mm256_shuffle_epi32(tmp_A1, 0xd8);
        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);

        tmp_A2 = _mm256_shuffle_epi32(tmp_A2, 0xd8);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        //V*S
        //V*S%12289, (U+V)%12289, (U-V)%12289
        ntt_core_lazymod(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    //level 10
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f_12289[m+i];
        tmp_S = _mm256_set_epi16(f_12289[512 + s + 15], f_12289[512 + s + 7], f_12289[512 + s + 11], f_12289[512 + s + 3], f_12289[512 + s + 13], f_12289[512 + s + 5], f_12289[512 + s + 9], f_12289[512 + s + 1],
                                 f_12289[512 + s + 14], f_12289[512 + s + 6], f_12289[512 + s + 10], f_12289[512 + s + 2], f_12289[512 + s + 12], f_12289[512 + s + 4], f_12289[512 + s + 8], f_12289[512 + s]);
        s += 16;
        //load
        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));
        //permute
        tmp_A1 = _mm256_shuffle_epi8(tmp_A1, tmp_shuffle16_pre);
        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);

        tmp_A2 = _mm256_shuffle_epi8(tmp_A2, tmp_shuffle16_pre);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        //V*S%12289, (U+V)%12289, (U-V)%12289
        ntt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    return 0;
}

int mq_poly_intt_12289(int16_t *a)
{
    int i, j, s;

    __m256i tmp_U, tmp_V, tmp_A1, tmp_A2, tmp_S;
    __m256i tmp_shuffle16_pos = _mm256_set_epi8(0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00,
                                        0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00);

    //level 10
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[512 + s + 15], fn_12289[512 + s + 7], fn_12289[512 + s + 11], fn_12289[512 + s + 3], fn_12289[512 + s + 13], fn_12289[512 + s + 5], fn_12289[512 + s + 9], fn_12289[512 + s + 1],
                                 fn_12289[512 + s + 14], fn_12289[512 + s + 6], fn_12289[512 + s + 10], fn_12289[512 + s + 2], fn_12289[512 + s + 12], fn_12289[512 + s + 4], fn_12289[512 + s + 8], fn_12289[512 + s]);
        s += 16;
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
        //  A1=(U+V-12289)%12289, A2=MontFull((U-V)*S)
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi8(tmp_U, tmp_shuffle16_pos);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi8(tmp_V, tmp_shuffle16_pos);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    //level 9
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[256 + s + 7], fn_12289[256 + s + 7], fn_12289[256 + s + 3], fn_12289[256 + s + 3], fn_12289[256 + s + 5], fn_12289[256 + s + 5], fn_12289[256 + s + 1], fn_12289[256 + s + 1],
                                 fn_12289[256 + s + 6], fn_12289[256 + s + 6], fn_12289[256 + s + 2], fn_12289[256 + s + 2], fn_12289[256 + s + 4], fn_12289[256 + s + 4], fn_12289[256 + s], fn_12289[256 + s]);
        s += 8;
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi32(tmp_U, 0xd8);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi32(tmp_V, 0xd8);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    //level 8
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[128 + s + 3], fn_12289[128 + s + 3], fn_12289[128 + s + 3], fn_12289[128 + s + 3], fn_12289[128 + s + 1], fn_12289[128 + s + 1], fn_12289[128 + s + 1], fn_12289[128 + s + 1],
                                 fn_12289[128 + s + 2], fn_12289[128 + s + 2], fn_12289[128 + s + 2], fn_12289[128 + s + 2], fn_12289[128 + s], fn_12289[128 + s], fn_12289[128 + s], fn_12289[128 + s]);
        s += 4;
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    //level 7
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1],
                                 fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s]);
        s += 2;
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }
    //level 6
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set1_epi16(fn_12289[32 + s++]);
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
       
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    //level 5
    s=0;
    for (i = 0; i < DIM_N; i += 64)
    {
        //load root
        tmp_S = _mm256_set1_epi16(fn_12289[16 +(s++)]);
        //load block 1
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 32));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 32), tmp_A2);

         //load block 2
        tmp_U = _mm256_load_si256((__m256i *)(a + i+16));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 48));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i+16), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 48), tmp_A2);
    }

    //level 4
    s=0;
    for (i = 0; i < DIM_N; i += 128)
    {
        //load root
        tmp_S = _mm256_set1_epi16(fn_12289[8 +(s++)]);
        for(j=0;j<64;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 64+j));
            //V*S%12289, (U+V)%12289, (U-V)%12289
            intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 64+j), tmp_A2);
        }
    }

     //level 3  m=4
    s=0;
    for (i = 0; i < DIM_N; i += 256)
    {
        //load root
        tmp_S = _mm256_set1_epi16(fn_12289[4 +(s++)]);
        for(j=0;j<128;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 128+j));
            //V*S%12289, (U+V)%12289, (U-V)%12289
            intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 128+j), tmp_A2);
        }
    }

    //level 2
    s=0;
    for (i = 0; i < DIM_N; i += 512)
    {
        //load root
        tmp_S = _mm256_set1_epi16(fn_12289[2 +(s++)]);
        for(j=0;j<256;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 256+j));
            //V*S%12289, (U+V)%12289, (U-V)%12289
            intt_core_lazymod(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 256+j), tmp_A2);
        }
    }

    //level 1  m=1
    tmp_S = _mm256_set1_epi16(983);
    for(j=0;j<512;j+=16)
    {
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a +j));
        tmp_V = _mm256_load_si256((__m256i *)(a  + 512+j));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core_final(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a +j), tmp_A1);
        _mm256_store_si256((__m256i *)(a  + 512+j), tmp_A2);
    }


    return 0;
}

int mq_poly_intt_12289_n(int16_t *a)
{
    int i, j, s;

    __m256i tmp_U, tmp_V, tmp_A1, tmp_A2, tmp_S;
    __m256i tmp_shuffle16_pos = _mm256_set_epi8(0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00,
                                        0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00);

    //level 10
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[512 + s + 15], fn_12289[512 + s + 7], fn_12289[512 + s + 11], fn_12289[512 + s + 3], fn_12289[512 + s + 13], fn_12289[512 + s + 5], fn_12289[512 + s + 9], fn_12289[512 + s + 1],
                                 fn_12289[512 + s + 14], fn_12289[512 + s + 6], fn_12289[512 + s + 10], fn_12289[512 + s + 2], fn_12289[512 + s + 12], fn_12289[512 + s + 4], fn_12289[512 + s + 8], fn_12289[512 + s]);
        s += 16;
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
        //  A1=(U+V-12289)%12289, A2=MontFull((U-V)*S)
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi8(tmp_U, tmp_shuffle16_pos);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi8(tmp_V, tmp_shuffle16_pos);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    //level 9
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[256 + s + 7], fn_12289[256 + s + 7], fn_12289[256 + s + 3], fn_12289[256 + s + 3], fn_12289[256 + s + 5], fn_12289[256 + s + 5], fn_12289[256 + s + 1], fn_12289[256 + s + 1],
                                 fn_12289[256 + s + 6], fn_12289[256 + s + 6], fn_12289[256 + s + 2], fn_12289[256 + s + 2], fn_12289[256 + s + 4], fn_12289[256 + s + 4], fn_12289[256 + s], fn_12289[256 + s]);
        s += 8;
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi32(tmp_U, 0xd8);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi32(tmp_V, 0xd8);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    //level 8
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[128 + s + 3], fn_12289[128 + s + 3], fn_12289[128 + s + 3], fn_12289[128 + s + 3], fn_12289[128 + s + 1], fn_12289[128 + s + 1], fn_12289[128 + s + 1], fn_12289[128 + s + 1],
                                 fn_12289[128 + s + 2], fn_12289[128 + s + 2], fn_12289[128 + s + 2], fn_12289[128 + s + 2], fn_12289[128 + s], fn_12289[128 + s], fn_12289[128 + s], fn_12289[128 + s]);
        s += 4;
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    //level 7
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1],
                                 fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s]);
        s += 2;
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //permute
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }
    //level 6
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set1_epi16(fn_12289[32 + s++]);
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
       
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    //level 5
    s=0;
    for (i = 0; i < DIM_N; i += 64)
    {
        //load root
        tmp_S = _mm256_set1_epi16(fn_12289[16 +(s++)]);
        //load block 1
        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 32));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 32), tmp_A2);

         //load block 2
        tmp_U = _mm256_load_si256((__m256i *)(a + i+16));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 48));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i+16), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 48), tmp_A2);
    }

    //level 4
    s=0;
    for (i = 0; i < DIM_N; i += 128)
    {
        //load root
        tmp_S = _mm256_set1_epi16(fn_12289[8 +(s++)]);
        for(j=0;j<64;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 64+j));
            //V*S%12289, (U+V)%12289, (U-V)%12289
            intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 64+j), tmp_A2);
        }
    }

     //level 3  m=4
    s=0;
    for (i = 0; i < DIM_N; i += 256)
    {
        //load root
        tmp_S = _mm256_set1_epi16(fn_12289[4 +(s++)]);
        for(j=0;j<128;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 128+j));
            //V*S%12289, (U+V)%12289, (U-V)%12289
            intt_core(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 128+j), tmp_A2);
        }
    }

    //level 2
    s=0;
    for (i = 0; i < DIM_N; i += 512)
    {
        //load root
        tmp_S = _mm256_set1_epi16(fn_12289[2 +(s++)]);
        for(j=0;j<256;j+=16)
        {
            //load
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 256+j));
            //V*S%12289, (U+V)%12289, (U-V)%12289
            intt_core_lazymod(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
            //store
            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 256+j), tmp_A2);
        }
    }

    //level 1  m=1
    tmp_S = _mm256_set1_epi16(983);
    for(j=0;j<512;j+=16)
    {
        //load
        tmp_U = _mm256_load_si256((__m256i *)(a +j));
        tmp_V = _mm256_load_si256((__m256i *)(a  + 512+j));
        //V*S%12289, (U+V)%12289, (U-V)%12289
        intt_core_final_n(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a +j), tmp_A1);
        _mm256_store_si256((__m256i *)(a  + 512+j), tmp_A2);
    }


    return 0;
}