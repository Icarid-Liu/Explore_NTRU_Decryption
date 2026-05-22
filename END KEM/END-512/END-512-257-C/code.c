#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "param.h"

static const uint64_t pk_table1_ct[7] = {
    1ULL,
    86ULL,
    7396ULL,
    636056ULL,
    54700816ULL,
    4704270176ULL,
    404567235136ULL
};

static const uint32_t pk_table2_ct[17] = {
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

/*
 * Return all-ones mask if x >= y, otherwise 0.
 * Precondition: x and y are small enough that x - y has a meaningful top bit.
 */
static inline uint32_t ct_ge_mask_u32(uint32_t x, uint32_t y)
{
    return 0U - (1U ^ ((x - y) >> 31));
}

static inline uint64_t ct_ge_mask_u64(uint64_t x, uint64_t y)
{
    return 0ULL - (uint64_t)(1U ^ (uint32_t)((x - y) >> 63));
}

static inline uint32_t ct_sub_if_ge_u32(uint32_t x, uint32_t y)
{
    uint32_t m = ct_ge_mask_u32(x, y);
    return x - (y & m);
}

static inline uint64_t ct_sub_if_ge_u64(uint64_t x, uint64_t y)
{
    uint64_t m = ct_ge_mask_u64(x, y);
    return x - (y & m);
}

/*
 * Canonical reduction modulo 258 for signed int16_t input.
 * Adds 128 * 258 = 33024 first, so negative coefficients are interpreted
 * canonically modulo 258.
 */
static inline uint16_t ct_mod258_i16(int16_t a)
{
    uint32_t x = (uint32_t)((int32_t)a + 33024);

    x = ct_sub_if_ge_u32(x, 33024U);
    x = ct_sub_if_ge_u32(x, 16512U);
    x = ct_sub_if_ge_u32(x, 8256U);
    x = ct_sub_if_ge_u32(x, 4128U);
    x = ct_sub_if_ge_u32(x, 2064U);
    x = ct_sub_if_ge_u32(x, 1032U);
    x = ct_sub_if_ge_u32(x, 516U);
    x = ct_sub_if_ge_u32(x, 258U);

    return (uint16_t)x;
}

static inline uint16_t ct_mod258_u32_small(uint32_t x)
{
    x = ct_sub_if_ge_u32(x, 4128U);
    x = ct_sub_if_ge_u32(x, 2064U);
    x = ct_sub_if_ge_u32(x, 1032U);
    x = ct_sub_if_ge_u32(x, 516U);
    x = ct_sub_if_ge_u32(x, 258U);

    return (uint16_t)x;
}

static inline uint16_t ct_mod86_u16(uint16_t a)
{
    uint32_t x = a;

    x = ct_sub_if_ge_u32(x, 172U);
    x = ct_sub_if_ge_u32(x, 86U);

    return (uint16_t)x;
}

static inline uint16_t ct_mod3_u16(uint16_t a)
{
    uint32_t x = a;

    x = ct_sub_if_ge_u32(x, 192U);
    x = ct_sub_if_ge_u32(x, 96U);
    x = ct_sub_if_ge_u32(x, 48U);
    x = ct_sub_if_ge_u32(x, 24U);
    x = ct_sub_if_ge_u32(x, 12U);
    x = ct_sub_if_ge_u32(x, 6U);
    x = ct_sub_if_ge_u32(x, 3U);

    return (uint16_t)x;
}

/*
 * x is always less than 86^7 in this decoder.
 * Division is implemented with reciprocal multiplication and one correction.
 */
static inline void ct_divmod86_u64(uint64_t x, uint64_t *q_out, uint16_t *r_out)
{
    const uint64_t magic = 0x02fa0be82fa0be82ULL; /* floor(2^64 / 86) */

    uint64_t q = (uint64_t)(((__uint128_t)x * magic) >> 64);
    uint64_t r = x - q * 86ULL;

    uint64_t m = ct_ge_mask_u64(r, 86ULL);
    r -= 86ULL & m;
    q += m & 1ULL;

    *q_out = q;
    *r_out = (uint16_t)r;
}

/*
 * Exact unsigned division by 3 using reciprocal multiplication.
 */
static inline void ct_divmod3_u32(uint32_t x, uint32_t *q_out, uint16_t *r_out)
{
    uint32_t q = (uint32_t)(((uint64_t)x * 0xAAAAAAABULL) >> 33);
    uint32_t r = x - q * 3U;

    uint32_t m = ct_ge_mask_u32(r, 3U);
    r -= 3U & m;
    q += m & 1U;

    *q_out = q;
    *r_out = (uint16_t)r;
}

static inline void store64_le(uint8_t *out, uint64_t x)
{
    out[0] = (uint8_t)(x);
    out[1] = (uint8_t)(x >> 8);
    out[2] = (uint8_t)(x >> 16);
    out[3] = (uint8_t)(x >> 24);
    out[4] = (uint8_t)(x >> 32);
    out[5] = (uint8_t)(x >> 40);
    out[6] = (uint8_t)(x >> 48);
    out[7] = (uint8_t)(x >> 56);
}

static inline uint64_t load64_le(const uint8_t *in)
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

static inline void store32_le(uint8_t *out, uint32_t x)
{
    out[0] = (uint8_t)(x);
    out[1] = (uint8_t)(x >> 8);
    out[2] = (uint8_t)(x >> 16);
    out[3] = (uint8_t)(x >> 24);
}

static inline uint32_t load32_le(const uint8_t *in)
{
    return ((uint32_t)in[0])
         | ((uint32_t)in[1] << 8)
         | ((uint32_t)in[2] << 16)
         | ((uint32_t)in[3] << 24);
}

void encode_pk(int16_t *c, uint8_t *code_c)
{
    size_t i, j, idx;
    uint64_t tmp[73];
    uint32_t tmp32[31];
    uint32_t t;

    /*
     * First component: residues modulo 86.
     * Pack 511 coefficients into 73 base-86 words, each containing 7 digits.
     */
    idx = 0;
    for (i = 0; i < DIM_N - 1; i += 7) {
        uint64_t acc = 0;

        for (j = 0; j < 7; j++) {
            uint16_t x = ct_mod258_i16(c[i + j]);
            uint16_t r = ct_mod86_u16(x);
            acc += (uint64_t)r * pk_table1_ct[j];
        }

        tmp[idx++] = acc;
    }

    /*
     * Store the lower 38 bits of tmp[0..71] through 48 little-endian words.
     */
    j = 0;
    for (i = 0; i < 72; i += 3) {
        uint64_t r0 = (tmp[i + 1] << 19) | (tmp[i] & 0x7FFFFULL);
        uint64_t r1 = (tmp[i + 2] << 19) | ((tmp[i] >> 19) & 0x7FFFFULL);

        store64_le(code_c + 8 * j, r0);
        store64_le(code_c + 8 * (j + 1), r1);
        j += 2;
    }

    /*
     * Store the remaining high bits of selected tmp words.
     */
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

    /*
     * Store tmp[72] in 6 bytes and the final modulo-86 residue in one byte.
     */
    for (i = 0; i < 6; i++) {
        code_c[idx++] = (uint8_t)(tmp[72] >> (8 * i));
    }

    {
        uint16_t x = ct_mod258_i16(c[DIM_N - 1]);
        code_c[idx++] = (uint8_t)ct_mod86_u16(x);
    }

    /*
     * Second component: residues modulo 3.
     * Pack 510 coefficients into 30 base-3 words, each containing 17 digits.
     * The last two coefficients are stored in tmp32[30].
     */
    idx = 0;
    for (i = 0; i < DIM_N - 2; i += 17) {
        uint32_t acc = 0;

        for (j = 0; j < 17; j++) {
            uint16_t x = ct_mod258_i16(c[i + j]);
            uint16_t r = ct_mod3_u16(x);
            acc += (uint32_t)r * pk_table2_ct[j];
        }

        tmp32[idx++] = acc;
    }

    tmp32[30] = 0;
    for (i = 0; i < 2; i++) {
        uint16_t x = ct_mod258_i16(c[DIM_N - 2 + i]);
        uint16_t r = ct_mod3_u16(x);
        tmp32[30] += (uint32_t)r * pk_table2_ct[i];
    }

    /*
     * Store tmp32[0..29] in 25 little-endian words.
     */
    j = 0;
    t = 0;
    for (i = 0; i < 30; i += 6) {
        uint32_t r0 = tmp32[i + 0] | ((tmp32[i + 5] & 31U) << 27);
        uint32_t r1 = tmp32[i + 1] | (((tmp32[i + 5] >> 5) & 31U) << 27);
        uint32_t r2 = tmp32[i + 2] | (((tmp32[i + 5] >> 10) & 31U) << 27);
        uint32_t r3 = tmp32[i + 3] | (((tmp32[i + 5] >> 15) & 31U) << 27);
        uint32_t r4 = tmp32[i + 4] | (((tmp32[i + 5] >> 20) & 31U) << 27);

        store32_le(code_c + 412 + 4 * j, r0);
        store32_le(code_c + 412 + 4 * (j + 1), r1);
        store32_le(code_c + 412 + 4 * (j + 2), r2);
        store32_le(code_c + 412 + 4 * (j + 3), r3);
        store32_le(code_c + 412 + 4 * (j + 4), r4);

        j += 5;
        t = (t << 2) | (tmp32[i + 5] >> 25);
    }

    code_c[512] = (uint8_t)t;
    code_c[513] = (uint8_t)((t >> 8) | (tmp32[30] << 2));
}

void decode_pk(uint8_t *code_c, int16_t *c)
{
    size_t i, j, idx;
    uint64_t tmp[73];
    uint32_t tmp32[31];
    uint32_t t;
    uint16_t r86[DIM_N];
    uint16_t r3[DIM_N];

    /*
     * Recover high bits for base-86 words.
     */
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

    /*
     * Recover tmp[72] from 6 little-endian bytes.
     */
    tmp[72] = 0;
    for (i = 0; i < 6; i++) {
        tmp[72] |= ((uint64_t)code_c[405 + i]) << (8 * i);
    }

    /*
     * Recover low bits from the first 384 bytes.
     */
    j = 0;
    for (i = 0; i < 72; i += 3) {
        uint64_t r0 = load64_le(code_c + 8 * j);
        uint64_t r1 = load64_le(code_c + 8 * (j + 1));

        tmp[i + 0] += ((r1 & 0x7FFFFULL) << 19) | (r0 & 0x7FFFFULL);
        tmp[i + 1]  = r0 >> 19;
        tmp[i + 2]  = r1 >> 19;

        j += 2;
    }

    /*
     * Extract base-86 digits.
     */
    idx = 0;
    for (i = 0; i < DIM_N - 1; i += 7) {
        uint64_t x = tmp[idx++];

        for (j = 0; j < 7; j++) {
            uint64_t q;
            uint16_t r;

            ct_divmod86_u64(x, &q, &r);
            r86[i + j] = r;
            x = q;
        }
    }

    /*
     * Recover base-3 packed words.
     */
    tmp32[30] = code_c[513] >> 2;
    t = (((uint32_t)code_c[513] & 3U) << 8) | (uint32_t)code_c[512];

    j = 0;
    for (i = 0; i < 30; i += 6) {
        uint32_t r0 = load32_le(code_c + 412 + 4 * j);
        uint32_t r1 = load32_le(code_c + 412 + 4 * (j + 1));
        uint32_t r2 = load32_le(code_c + 412 + 4 * (j + 2));
        uint32_t r3 = load32_le(code_c + 412 + 4 * (j + 3));
        uint32_t r4 = load32_le(code_c + 412 + 4 * (j + 4));

        tmp32[i + 0] = r0 & 0x07FFFFFFU;
        tmp32[i + 1] = r1 & 0x07FFFFFFU;
        tmp32[i + 2] = r2 & 0x07FFFFFFU;
        tmp32[i + 3] = r3 & 0x07FFFFFFU;
        tmp32[i + 4] = r4 & 0x07FFFFFFU;
        tmp32[i + 5] =
            ((r0 >> 27)      ) |
            ((r1 >> 27) << 5 ) |
            ((r2 >> 27) << 10) |
            ((r3 >> 27) << 15) |
            ((r4 >> 27) << 20);

        j += 5;
    }

    tmp32[5]  += ((t >> 8)       ) << 25;
    tmp32[11] += ((t >> 6) & 3U  ) << 25;
    tmp32[17] += ((t >> 4) & 3U  ) << 25;
    tmp32[23] += ((t >> 2) & 3U  ) << 25;
    tmp32[29] += ((t     ) & 3U  ) << 25;

    /*
     * Extract last two base-3 digits.
     */
    {
        uint32_t x = tmp32[30];

        for (i = 0; i < 2; i++) {
            uint32_t q;
            uint16_t r;

            ct_divmod3_u32(x, &q, &r);
            r3[DIM_N - 2 + i] = r;
            x = q;
        }
    }

    /*
     * Extract remaining base-3 digits.
     */
    idx = 0;
    for (i = 0; i < DIM_N - 2; i += 17) {
        uint32_t x = tmp32[idx++];

        for (j = 0; j < 17; j++) {
            uint32_t q;
            uint16_t r;

            ct_divmod3_u32(x, &q, &r);
            r3[i + j] = r;
            x = q;
        }
    }

    /*
     * CRT reconstruction modulo 258:
     *
     *   x = r86 * 3 * 29 + r3 * 86 * 2 mod 258
     *     = r86 * 87 + r3 * 172 mod 258
     */
    for (i = 0; i < DIM_N; i++) {
        uint32_t x = (uint32_t)r86[i] * 87U + (uint32_t)r3[i] * 172U;
        c[i] = (int16_t)ct_mod258_u32_small(x);
    }
}

void encode_c(int16_t *c, uint8_t *code_c)
{
    int i, j;

    j = 0;
    for(i = 0; i < DIM_N; i += 4)
    {
        code_c[j] = ((c[i + 3] & 3) << 6) | c[i];
        code_c[j + 1] = (((c[i + 3] >> 2) & 3) << 6) | c[i + 1];
        code_c[j + 2] = (((c[i + 3] >> 4) & 3) << 6) | c[i + 2];
        j += 3;
    }

}

void decode_c(uint8_t *code_c, int16_t *c)
{
    int i, j;

    j = 0;
    for(i = 0; i < CIPHLEN; i += 3)
    {
        c[j] = code_c[i] & 0x3f;
        c[j + 1] = code_c[i + 1] & 0x3f;
        c[j + 2] = code_c[i + 2] & 0x3f;
        c[j + 3] = (((code_c[i + 2] >> 6) & 3) << 4) | (((code_c[i + 1] >> 6) & 3) << 2) | ((code_c[i] >> 6) & 3);
        j += 4;
    }
}