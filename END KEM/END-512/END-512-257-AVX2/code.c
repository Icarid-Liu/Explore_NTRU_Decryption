
#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include "param.h"

static const uint64_t pk_table1_avx2[7] = {
    1ULL,
    86ULL,
    7396ULL,
    636056ULL,
    54700816ULL,
    4704270176ULL,
    404567235136ULL
};

static const uint32_t pk_table2_avx2[17] = {
    1U,
    3U,
    9U,
    27U,
    81U,
    243U,
    729U,
    2187U,
    6561U,
    19683U,
    59049U,
    177147U,
    531441U,
    1594323U,
    4782969U,
    14348907U,
    43046721U
};

static inline uint32_t ct_ge_mask_u32_avx2(uint32_t x, uint32_t y)
{
    return 0U - (1U ^ ((x - y) >> 31));
}

static inline uint64_t ct_ge_mask_u64_avx2(uint64_t x, uint64_t y)
{
    return 0ULL - (uint64_t)(1U ^ (uint32_t)((x - y) >> 63));
}

static inline uint32_t ct_sub_if_ge_u32_avx2(uint32_t x, uint32_t y)
{
    uint32_t m = ct_ge_mask_u32_avx2(x, y);
    return x - (y & m);
}

static inline uint64_t ct_sub_if_ge_u64_avx2(uint64_t x, uint64_t y)
{
    uint64_t m = ct_ge_mask_u64_avx2(x, y);
    return x - (y & m);
}

static inline uint16_t ct_mod258_u32_small_avx2(uint32_t x)
{
    x = ct_sub_if_ge_u32_avx2(x, 4128U);
    x = ct_sub_if_ge_u32_avx2(x, 2064U);
    x = ct_sub_if_ge_u32_avx2(x, 1032U);
    x = ct_sub_if_ge_u32_avx2(x, 516U);
    x = ct_sub_if_ge_u32_avx2(x, 258U);
    return (uint16_t)x;
}

static inline void ct_divmod86_u64_avx2(uint64_t x, uint64_t *q_out, uint16_t *r_out)
{
    uint64_t magic = 0x02fa0be82fa0be82ULL; /* floor(2^64 / 86) */

    uint64_t q = (uint64_t)(((__uint128_t)x * magic) >> 64);
    uint64_t r = x - q * 86ULL;

    uint64_t m = ct_ge_mask_u64_avx2(r, 86ULL);
    r -= 86ULL & m;
    q += m & 1ULL;

    *q_out = q;
    *r_out = (uint16_t)r;
}

static inline void ct_divmod3_u32_avx2(uint32_t x, uint32_t *q_out, uint16_t *r_out)
{
    uint32_t q = (uint32_t)(((uint64_t)x * 0xAAAAAAABULL) >> 33);
    uint32_t r = x - q * 3U;

    uint32_t m = ct_ge_mask_u32_avx2(r, 3U);
    r -= 3U & m;
    q += m & 1U;

    *q_out = q;
    *r_out = (uint16_t)r;
}

static inline void store64_le_avx2(uint8_t *out, uint64_t x)
{
    out[0] = (uint8_t)x;
    out[1] = (uint8_t)(x >> 8);
    out[2] = (uint8_t)(x >> 16);
    out[3] = (uint8_t)(x >> 24);
    out[4] = (uint8_t)(x >> 32);
    out[5] = (uint8_t)(x >> 40);
    out[6] = (uint8_t)(x >> 48);
    out[7] = (uint8_t)(x >> 56);
}

static inline uint64_t load64_le_avx2(uint8_t *in)
{
    return ((uint64_t)in[0])
         | ((uint64_t)in[1] << 8)
         | ((uint64_t)in[2] << 16)
         | ((uint64_t)in[3] << 24)
         | ((uint64_t)in[4] << 32)
         | ((uint64_t)in[5] << 40)
         | ((uint64_t)in[6] << 48)
         | ((uint64_t)in[7] << 56);
}

static inline void store32_le_avx2(uint8_t *out, uint32_t x)
{
    out[0] = (uint8_t)x;
    out[1] = (uint8_t)(x >> 8);
    out[2] = (uint8_t)(x >> 16);
    out[3] = (uint8_t)(x >> 24);
}

static inline uint32_t load32_le_avx2(uint8_t *in)
{
    return ((uint32_t)in[0])
         | ((uint32_t)in[1] << 8)
         | ((uint32_t)in[2] << 16)
         | ((uint32_t)in[3] << 24);
}

static inline __m256i avx2_sub_if_ge_u32(__m256i x, uint32_t y)
{
    __m256i yy = _mm256_set1_epi32((int)y);
    __m256i ym1 = _mm256_set1_epi32((int)(y - 1U));
    __m256i mask = _mm256_cmpgt_epi32(x, ym1); /* x >= y, signed is safe for these ranges */
    return _mm256_sub_epi32(x, _mm256_and_si256(yy, mask));
}

static inline __m256i avx2_sub_if_ge_u16(__m256i x, uint16_t y)
{
    __m256i yy = _mm256_set1_epi16((int16_t)y);
    __m256i ym1 = _mm256_set1_epi16((int16_t)(y - 1U));
    __m256i mask = _mm256_cmpgt_epi16(x, ym1); /* x >= y, signed is safe for x <= 257 */
    return _mm256_sub_epi16(x, _mm256_and_si256(yy, mask));
}

static inline __m256i avx2_mod258_8x32(__m256i x)
{
    x = _mm256_add_epi32(x, _mm256_set1_epi32(33024));
    x = avx2_sub_if_ge_u32(x, 33024U);
    x = avx2_sub_if_ge_u32(x, 16512U);
    x = avx2_sub_if_ge_u32(x, 8256U);
    x = avx2_sub_if_ge_u32(x, 4128U);
    x = avx2_sub_if_ge_u32(x, 2064U);
    x = avx2_sub_if_ge_u32(x, 1032U);
    x = avx2_sub_if_ge_u32(x, 516U);
    x = avx2_sub_if_ge_u32(x, 258U);
    return x;
}

static inline __m256i avx2_mod258_i16x16(__m256i a)
{
    __m128i lo128 = _mm256_castsi256_si128(a);
    __m128i hi128 = _mm256_extracti128_si256(a, 1);

    __m256i lo = _mm256_cvtepi16_epi32(lo128);
    __m256i hi = _mm256_cvtepi16_epi32(hi128);

    lo = avx2_mod258_8x32(lo);
    hi = avx2_mod258_8x32(hi);

    __m256i packed = _mm256_packus_epi32(lo, hi);
    return _mm256_permute4x64_epi64(packed, 0xD8);
}

static inline __m256i avx2_mod86_u16x16(__m256i x)
{
    x = avx2_sub_if_ge_u16(x, 172U);
    x = avx2_sub_if_ge_u16(x, 86U);
    return x;
}

static inline __m256i avx2_mod3_u16x16(__m256i x)
{
    x = avx2_sub_if_ge_u16(x, 192U);
    x = avx2_sub_if_ge_u16(x, 96U);
    x = avx2_sub_if_ge_u16(x, 48U);
    x = avx2_sub_if_ge_u16(x, 24U);
    x = avx2_sub_if_ge_u16(x, 12U);
    x = avx2_sub_if_ge_u16(x, 6U);
    x = avx2_sub_if_ge_u16(x, 3U);
    return x;
}

/*
 * Compute both residue streams with AVX2:
 *   r86[i] = c[i] mod 258 mod 86
 *   r3[i]  = c[i] mod 258 mod 3
 */
static inline void precompute_pk_residues_avx2(int16_t *c,
                                               uint16_t *r86,
                                               uint16_t *r3)
{
    size_t i;

    for (i = 0; i < DIM_N; i += 16) {
        __m256i v = _mm256_loadu_si256((__m256i *)(void *)(c + i));
        __m256i x = avx2_mod258_i16x16(v);
        __m256i a = avx2_mod86_u16x16(x);
        __m256i b = avx2_mod3_u16x16(x);

        _mm256_storeu_si256((__m256i *)(void *)(r86 + i), a);
        _mm256_storeu_si256((__m256i *)(void *)(r3 + i), b);
    }
}

/*
 * AVX2-accelerated pk encoder.
 * The expensive coefficient normalization and mod-86/mod-3 projections are vectorized.
 * The mixed-radix packing is kept scalar because it is a short fixed-layout bit packing step.
 */
void encode_pk(int16_t *c, uint8_t *code_c)
{
    size_t i, j, idx;

    uint64_t tmp[73];
    uint32_t tmp32[31];
    uint32_t t;

    __attribute__((aligned(32))) uint16_t r86[DIM_N];
    __attribute__((aligned(32))) uint16_t r3[DIM_N];

    precompute_pk_residues_avx2(c, r86, r3);

    idx = 0;
    for (i = 0; i < DIM_N - 1; i += 7) {
        uint64_t acc =
            (uint64_t)r86[i + 0] * pk_table1_avx2[0] +
            (uint64_t)r86[i + 1] * pk_table1_avx2[1] +
            (uint64_t)r86[i + 2] * pk_table1_avx2[2] +
            (uint64_t)r86[i + 3] * pk_table1_avx2[3] +
            (uint64_t)r86[i + 4] * pk_table1_avx2[4] +
            (uint64_t)r86[i + 5] * pk_table1_avx2[5] +
            (uint64_t)r86[i + 6] * pk_table1_avx2[6];

        tmp[idx++] = acc;
    }

    j = 0;
    for (i = 0; i < 72; i += 3) {
        uint64_t r0 = (tmp[i + 1] << 19) | (tmp[i] & 0x7FFFFULL);
        uint64_t r1 = (tmp[i + 2] << 19) | ((tmp[i] >> 19) & 0x7FFFFULL);

        store64_le_avx2(code_c + 8 * j, r0);
        store64_le_avx2(code_c + 8 * (j + 1), r1);
        j += 2;
    }

    idx = 384;
    for (i = 0; i < 72; i += 24) {
        code_c[idx + 0] = (uint8_t)(((tmp[i + 0]  >> 38) << 1) | ((tmp[i + 21] >> 38) & 1ULL));
        code_c[idx + 1] = (uint8_t)(((tmp[i + 3]  >> 38) << 1) | ((tmp[i + 21] >> 39) & 1ULL));
        code_c[idx + 2] = (uint8_t)(((tmp[i + 6]  >> 38) << 1) | ((tmp[i + 21] >> 40) & 1ULL));
        code_c[idx + 3] = (uint8_t)(((tmp[i + 9]  >> 38) << 1) | ((tmp[i + 21] >> 41) & 1ULL));
        code_c[idx + 4] = (uint8_t)(((tmp[i + 12] >> 38) << 1) | ((tmp[i + 21] >> 42) & 1ULL));
        code_c[idx + 5] = (uint8_t)(((tmp[i + 15] >> 38) << 1) | ((tmp[i + 21] >> 43) & 1ULL));
        code_c[idx + 6] = (uint8_t)(((tmp[i + 18] >> 38) << 1) | ((tmp[i + 21] >> 44) & 1ULL));
        idx += 7;
    }

    for (i = 0; i < 6; i++) {
        code_c[idx++] = (uint8_t)(tmp[72] >> (8 * i));
    }
    code_c[idx++] = (uint8_t)r86[DIM_N - 1];

    idx = 0;
    for (i = 0; i < DIM_N - 2; i += 17) {
        uint32_t acc =
            (uint32_t)r3[i + 0]  * pk_table2_avx2[0]  +
            (uint32_t)r3[i + 1]  * pk_table2_avx2[1]  +
            (uint32_t)r3[i + 2]  * pk_table2_avx2[2]  +
            (uint32_t)r3[i + 3]  * pk_table2_avx2[3]  +
            (uint32_t)r3[i + 4]  * pk_table2_avx2[4]  +
            (uint32_t)r3[i + 5]  * pk_table2_avx2[5]  +
            (uint32_t)r3[i + 6]  * pk_table2_avx2[6]  +
            (uint32_t)r3[i + 7]  * pk_table2_avx2[7]  +
            (uint32_t)r3[i + 8]  * pk_table2_avx2[8]  +
            (uint32_t)r3[i + 9]  * pk_table2_avx2[9]  +
            (uint32_t)r3[i + 10] * pk_table2_avx2[10] +
            (uint32_t)r3[i + 11] * pk_table2_avx2[11] +
            (uint32_t)r3[i + 12] * pk_table2_avx2[12] +
            (uint32_t)r3[i + 13] * pk_table2_avx2[13] +
            (uint32_t)r3[i + 14] * pk_table2_avx2[14] +
            (uint32_t)r3[i + 15] * pk_table2_avx2[15] +
            (uint32_t)r3[i + 16] * pk_table2_avx2[16];

        tmp32[idx++] = acc;
    }

    tmp32[30] = (uint32_t)r3[DIM_N - 2] + (uint32_t)r3[DIM_N - 1] * 3U;

    j = 0;
    t = 0;
    for (i = 0; i < 30; i += 6) {
        uint32_t r0 = tmp32[i + 0] | ((tmp32[i + 5] & 31U) << 27);
        uint32_t r1 = tmp32[i + 1] | (((tmp32[i + 5] >> 5) & 31U) << 27);
        uint32_t r2 = tmp32[i + 2] | (((tmp32[i + 5] >> 10) & 31U) << 27);
        uint32_t r3w = tmp32[i + 3] | (((tmp32[i + 5] >> 15) & 31U) << 27);
        uint32_t r4 = tmp32[i + 4] | (((tmp32[i + 5] >> 20) & 31U) << 27);

        store32_le_avx2(code_c + 412 + 4 * j, r0);
        store32_le_avx2(code_c + 412 + 4 * (j + 1), r1);
        store32_le_avx2(code_c + 412 + 4 * (j + 2), r2);
        store32_le_avx2(code_c + 412 + 4 * (j + 3), r3w);
        store32_le_avx2(code_c + 412 + 4 * (j + 4), r4);

        j += 5;
        t = (t << 2) | (tmp32[i + 5] >> 25);
    }

    code_c[512] = (uint8_t)t;
    code_c[513] = (uint8_t)((t >> 8) | (tmp32[30] << 2));
}

static inline __m256i avx2_mod258_8x32_small(__m256i x)
{
    x = avx2_sub_if_ge_u32(x, 4128U);
    x = avx2_sub_if_ge_u32(x, 2064U);
    x = avx2_sub_if_ge_u32(x, 1032U);
    x = avx2_sub_if_ge_u32(x, 516U);
    x = avx2_sub_if_ge_u32(x, 258U);
    return x;
}

static inline void crt_reconstruct_avx2(uint16_t *r86,
                                        uint16_t *r3,
                                        int16_t *c)
{
    size_t i;

    __m256i mul86 = _mm256_set1_epi32(87);
    __m256i mul3  = _mm256_set1_epi32(172);

    for (i = 0; i < DIM_N; i += 16) {
        __m256i a16 = _mm256_loadu_si256((__m256i *)(void *)(r86 + i));
        __m256i b16 = _mm256_loadu_si256((__m256i *)(void *)(r3 + i));

        __m256i a_lo = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a16));
        __m256i a_hi = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a16, 1));
        __m256i b_lo = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(b16));
        __m256i b_hi = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(b16, 1));

        __m256i x_lo = _mm256_add_epi32(_mm256_mullo_epi32(a_lo, mul86),
                                        _mm256_mullo_epi32(b_lo, mul3));
        __m256i x_hi = _mm256_add_epi32(_mm256_mullo_epi32(a_hi, mul86),
                                        _mm256_mullo_epi32(b_hi, mul3));

        x_lo = avx2_mod258_8x32_small(x_lo);
        x_hi = avx2_mod258_8x32_small(x_hi);

        __m256i packed = _mm256_packus_epi32(x_lo, x_hi);
        packed = _mm256_permute4x64_epi64(packed, 0xD8);

        _mm256_storeu_si256((__m256i *)(void *)(c + i), packed);
    }
}

/*
 * AVX2-accelerated pk decoder.
 * The base-86/base-3 digit extraction is scalar due to reciprocal division,
 * while CRT reconstruction is vectorized.
 */
void decode_pk(uint8_t *code_c, int16_t *c)
{
    size_t i, j, idx;

    uint64_t tmp[73];
    uint32_t tmp32[31];
    uint32_t t;

    __attribute__((aligned(32))) uint16_t r86[DIM_N];
    __attribute__((aligned(32))) uint16_t r3[DIM_N];

    for (i = 0; i < 73; i++) {
        tmp[i] = 0;
    }

    r86[DIM_N - 1] = code_c[411];

    idx = 384;
    for (i = 0; i < 72; i += 24) {
        tmp[i + 0]  = ((uint64_t)(code_c[idx + 0] >> 1)) << 38;
        tmp[i + 3]  = ((uint64_t)(code_c[idx + 1] >> 1)) << 38;
        tmp[i + 6]  = ((uint64_t)(code_c[idx + 2] >> 1)) << 38;
        tmp[i + 9]  = ((uint64_t)(code_c[idx + 3] >> 1)) << 38;
        tmp[i + 12] = ((uint64_t)(code_c[idx + 4] >> 1)) << 38;
        tmp[i + 15] = ((uint64_t)(code_c[idx + 5] >> 1)) << 38;
        tmp[i + 18] = ((uint64_t)(code_c[idx + 6] >> 1)) << 38;

        tmp[i + 21] =
            ((uint64_t)(code_c[idx + 0] & 1U) << 38) |
            ((uint64_t)(code_c[idx + 1] & 1U) << 39) |
            ((uint64_t)(code_c[idx + 2] & 1U) << 40) |
            ((uint64_t)(code_c[idx + 3] & 1U) << 41) |
            ((uint64_t)(code_c[idx + 4] & 1U) << 42) |
            ((uint64_t)(code_c[idx + 5] & 1U) << 43) |
            ((uint64_t)(code_c[idx + 6] & 1U) << 44);

        idx += 7;
    }

    tmp[72] = 0;
    for (i = 0; i < 6; i++) {
        tmp[72] |= ((uint64_t)code_c[405 + i]) << (8 * i);
    }

    j = 0;
    for (i = 0; i < 72; i += 3) {
        uint64_t r0 = load64_le_avx2(code_c + 8 * j);
        uint64_t r1 = load64_le_avx2(code_c + 8 * (j + 1));

        tmp[i + 0] += ((r1 & 0x7FFFFULL) << 19) | (r0 & 0x7FFFFULL);
        tmp[i + 1]  = r0 >> 19;
        tmp[i + 2]  = r1 >> 19;

        j += 2;
    }

    idx = 0;
    for (i = 0; i < DIM_N - 1; i += 7) {
        uint64_t x = tmp[idx++];

        uint64_t q;
        uint16_t r;

        ct_divmod86_u64_avx2(x, &q, &r); r86[i + 0] = r; x = q;
        ct_divmod86_u64_avx2(x, &q, &r); r86[i + 1] = r; x = q;
        ct_divmod86_u64_avx2(x, &q, &r); r86[i + 2] = r; x = q;
        ct_divmod86_u64_avx2(x, &q, &r); r86[i + 3] = r; x = q;
        ct_divmod86_u64_avx2(x, &q, &r); r86[i + 4] = r; x = q;
        ct_divmod86_u64_avx2(x, &q, &r); r86[i + 5] = r; x = q;
        ct_divmod86_u64_avx2(x, &q, &r); r86[i + 6] = r;
    }

    tmp32[30] = code_c[513] >> 2;
    t = (((uint32_t)code_c[513] & 3U) << 8) | (uint32_t)code_c[512];

    j = 0;
    for (i = 0; i < 30; i += 6) {
        uint32_t r0 = load32_le_avx2(code_c + 412 + 4 * j);
        uint32_t r1 = load32_le_avx2(code_c + 412 + 4 * (j + 1));
        uint32_t r2 = load32_le_avx2(code_c + 412 + 4 * (j + 2));
        uint32_t r3w = load32_le_avx2(code_c + 412 + 4 * (j + 3));
        uint32_t r4 = load32_le_avx2(code_c + 412 + 4 * (j + 4));

        tmp32[i + 0] = r0 & 0x07FFFFFFU;
        tmp32[i + 1] = r1 & 0x07FFFFFFU;
        tmp32[i + 2] = r2 & 0x07FFFFFFU;
        tmp32[i + 3] = r3w & 0x07FFFFFFU;
        tmp32[i + 4] = r4 & 0x07FFFFFFU;
        tmp32[i + 5] =
            ((r0  >> 27)      ) |
            ((r1  >> 27) << 5 ) |
            ((r2  >> 27) << 10) |
            ((r3w >> 27) << 15) |
            ((r4  >> 27) << 20);

        j += 5;
    }

    tmp32[5]  += ((t >> 8)       ) << 25;
    tmp32[11] += ((t >> 6) & 3U  ) << 25;
    tmp32[17] += ((t >> 4) & 3U  ) << 25;
    tmp32[23] += ((t >> 2) & 3U  ) << 25;
    tmp32[29] += ((t     ) & 3U  ) << 25;

    {
        uint32_t x = tmp32[30];
        uint32_t q;
        uint16_t r;

        ct_divmod3_u32_avx2(x, &q, &r); r3[DIM_N - 2] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[DIM_N - 1] = r;
    }

    idx = 0;
    for (i = 0; i < DIM_N - 2; i += 17) {
        uint32_t x = tmp32[idx++];

        uint32_t q;
        uint16_t r;

        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 0] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 1] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 2] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 3] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 4] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 5] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 6] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 7] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 8] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 9] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 10] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 11] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 12] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 13] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 14] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 15] = r; x = q;
        ct_divmod3_u32_avx2(x, &q, &r); r3[i + 16] = r;
    }

    crt_reconstruct_avx2(r86, r3, c);
}

static inline uint32_t pack4_u6_avx2_scalar(uint32_t c0, uint32_t c1,
                                            uint32_t c2, uint32_t c3)
{
    return ((c0 & 0x3FU)      ) |
           ((c3 & 0x03U) << 6 ) |
           ((c1 & 0x3FU) << 8 ) |
           ((c3 & 0x0CU) << 12) |
           ((c2 & 0x3FU) << 16) |
           ((c3 & 0x30U) << 18);
}

/*
 * Fast constant-layout c encoder.
 * This routine is intentionally unrolled and writes 24-bit groups.
 */
void encode_c(int16_t *c, uint8_t *code_c)
{
    size_t i, j = 0;

    for (i = 0; i < DIM_N; i += 16) {
        uint32_t w0 = pack4_u6_avx2_scalar((uint16_t)c[i + 0],  (uint16_t)c[i + 1],
                                           (uint16_t)c[i + 2],  (uint16_t)c[i + 3]);
        uint32_t w1 = pack4_u6_avx2_scalar((uint16_t)c[i + 4],  (uint16_t)c[i + 5],
                                           (uint16_t)c[i + 6],  (uint16_t)c[i + 7]);
        uint32_t w2 = pack4_u6_avx2_scalar((uint16_t)c[i + 8],  (uint16_t)c[i + 9],
                                           (uint16_t)c[i + 10], (uint16_t)c[i + 11]);
        uint32_t w3 = pack4_u6_avx2_scalar((uint16_t)c[i + 12], (uint16_t)c[i + 13],
                                           (uint16_t)c[i + 14], (uint16_t)c[i + 15]);

        code_c[j + 0] = (uint8_t)w0;
        code_c[j + 1] = (uint8_t)(w0 >> 8);
        code_c[j + 2] = (uint8_t)(w0 >> 16);

        code_c[j + 3] = (uint8_t)w1;
        code_c[j + 4] = (uint8_t)(w1 >> 8);
        code_c[j + 5] = (uint8_t)(w1 >> 16);

        code_c[j + 6] = (uint8_t)w2;
        code_c[j + 7] = (uint8_t)(w2 >> 8);
        code_c[j + 8] = (uint8_t)(w2 >> 16);

        code_c[j + 9]  = (uint8_t)w3;
        code_c[j + 10] = (uint8_t)(w3 >> 8);
        code_c[j + 11] = (uint8_t)(w3 >> 16);

        j += 12;
    }
}

static inline void unpack4_u6_avx2_scalar(uint8_t *in, int16_t *out)
{
    uint32_t b0 = in[0];
    uint32_t b1 = in[1];
    uint32_t b2 = in[2];

    out[0] = (int16_t)(b0 & 0x3FU);
    out[1] = (int16_t)(b1 & 0x3FU);
    out[2] = (int16_t)(b2 & 0x3FU);
    out[3] = (int16_t)(((b2 >> 6) << 4) | (((b1 >> 6) & 3U) << 2) | (b0 >> 6));
}

/*
 * Fast c decoder.
 */
void decode_c(uint8_t *code_c, int16_t *c)
{
    size_t i, j = 0;

    for (i = 0; i < CIPHLEN; i += 12) {
        unpack4_u6_avx2_scalar(code_c + i + 0, c + j + 0);
        unpack4_u6_avx2_scalar(code_c + i + 3, c + j + 4);
        unpack4_u6_avx2_scalar(code_c + i + 6, c + j + 8);
        unpack4_u6_avx2_scalar(code_c + i + 9, c + j + 12);
        j += 16;
    }
}
