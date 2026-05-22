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

static const uint64_t c_encode_table_ct[8] = {
    1ULL,
    90ULL,
    8100ULL,
    729000ULL,
    65610000ULL,
    5904900000ULL,
    531441000000ULL,
    47829690000000ULL
};

/* ------------------------------------------------------------------------- */
/* Constant-time-style scalar helpers                                         */
/* ------------------------------------------------------------------------- */

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
 * This maps negative values to their canonical representative modulo 258.
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
 * x < 86^7 in decode_pk.
 * Reciprocal division by 86 with branchless correction.
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
 * x < 3^17 in decode_pk, or x < 3^4 for the final word.
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

/*
 * x < 90^8 in decode_c.
 * Reciprocal division by 90 with branchless correction.
 */
static inline void ct_divmod90_u64(uint64_t x, uint64_t *q_out, uint16_t *r_out)
{
    const uint64_t magic = 0x02d82d82d82d82d8ULL; /* floor(2^64 / 90) */

    uint64_t q = (uint64_t)(((__uint128_t)x * magic) >> 64);
    uint64_t r = x - q * 90ULL;

    uint64_t m = ct_ge_mask_u64(r, 90ULL);
    r -= 90ULL & m;
    q += m & 1ULL;

    *q_out = q;
    *r_out = (uint16_t)r;
}

/* ------------------------------------------------------------------------- */
/* Explicit little-endian I/O                                                 */
/* ------------------------------------------------------------------------- */

static inline void store64_le(uint8_t *out, uint64_t x)
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
    out[0] = (uint8_t)x;
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

/* ------------------------------------------------------------------------- */
/* Public-key codec                                                           */
/* ------------------------------------------------------------------------- */

void encode_pk(int16_t *c, uint8_t *code_c)
{
    size_t i, j, idx;
    uint64_t tmp[147];
    uint32_t tmp32[61];
    uint32_t tail;

    /*
     * First component: residues modulo 86.
     * 1022 coefficients -> 146 base-86 words with 7 digits each.
     * The last 2 coefficients are packed into tmp[146].
     */
    idx = 0;
    for (i = 0; i < DIM_N - 2; i += 7) {
        uint64_t acc = 0;

        for (j = 0; j < 7; j++) {
            uint16_t x = ct_mod258_i16(c[i + j]);
            uint16_t r = ct_mod86_u16(x);
            acc += (uint64_t)r * pk_table1_ct[j];
        }

        tmp[idx++] = acc;
    }

    tmp[146] = 0;
    for (i = 0; i < 2; i++) {
        uint16_t x = ct_mod258_i16(c[DIM_N - 2 + i]);
        uint16_t r = ct_mod86_u16(x);
        tmp[146] += (uint64_t)r * pk_table1_ct[i];
    }

    /*
     * Store lower 38 bits of tmp[0..143] into 96 little-endian uint64 words.
     */
    j = 0;
    for (i = 0; i < 144; i += 3) {
        uint64_t r0 = (tmp[i + 1] << 19) | (tmp[i] & 0x7FFFFULL);
        uint64_t r1 = (tmp[i + 2] << 19) | ((tmp[i] >> 19) & 0x7FFFFULL);

        store64_le(code_c + 8 * j, r0);
        store64_le(code_c + 8 * (j + 1), r1);
        j += 2;
    }

    /*
     * Store high bits for tmp[0..143].
     */
    idx = 768;
    for (i = 0; i < 144; i += 24) {
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
     * Store tmp[144], tmp[145], and tmp[146] tail.
     */
    for (i = 0; i < 5; i++) {
        code_c[idx++] = (uint8_t)(tmp[144] >> (8 * i));
    }
    for (i = 0; i < 5; i++) {
        code_c[idx++] = (uint8_t)(tmp[145] >> (8 * i));
    }

    code_c[idx++] = (uint8_t)((((tmp[144] >> 40) & 31ULL) << 3) | (tmp[146] & 7ULL));
    code_c[idx++] = (uint8_t)((((tmp[145] >> 40) & 31ULL) << 3) | ((tmp[146] >> 3) & 7ULL));
    code_c[idx++] = (uint8_t)(tmp[146] >> 6);

    /*
     * Second component: residues modulo 3.
     * 1020 coefficients -> 60 base-3 words with 17 digits each.
     * The last 4 coefficients are packed into tmp32[60].
     */
    idx = 0;
    for (i = 0; i < DIM_N - 4; i += 17) {
        uint32_t acc = 0;

        for (j = 0; j < 17; j++) {
            uint16_t x = ct_mod258_i16(c[i + j]);
            uint16_t r = ct_mod3_u16(x);
            acc += (uint32_t)r * pk_table2_ct[j];
        }

        tmp32[idx++] = acc;
    }

    tmp32[60] = 0;
    for (i = 0; i < 4; i++) {
        uint16_t x = ct_mod258_i16(c[DIM_N - 4 + i]);
        uint16_t r = ct_mod3_u16(x);
        tmp32[60] += (uint32_t)r * pk_table2_ct[i];
    }

    /*
     * Store tmp32[0..59] in 50 little-endian uint32 words.
     * Store the high 2-bit fragments in the high part of tail.
     */
    j = 0;
    tail = tmp32[60];

    for (i = 0; i < 60; i += 6) {
        uint32_t r0 = tmp32[i + 0] | ((tmp32[i + 5] & 31U) << 27);
        uint32_t r1 = tmp32[i + 1] | (((tmp32[i + 5] >> 5) & 31U) << 27);
        uint32_t r2 = tmp32[i + 2] | (((tmp32[i + 5] >> 10) & 31U) << 27);
        uint32_t r3 = tmp32[i + 3] | (((tmp32[i + 5] >> 15) & 31U) << 27);
        uint32_t r4 = tmp32[i + 4] | (((tmp32[i + 5] >> 20) & 31U) << 27);

        store32_le(code_c + 823 + 4 * j, r0);
        store32_le(code_c + 823 + 4 * (j + 1), r1);
        store32_le(code_c + 823 + 4 * (j + 2), r2);
        store32_le(code_c + 823 + 4 * (j + 3), r3);
        store32_le(code_c + 823 + 4 * (j + 4), r4);

        j += 5;
        tail = (tail << 2) | (tmp32[i + 5] >> 25);
    }

    code_c[1023] = (uint8_t)tail;
    code_c[1024] = (uint8_t)(tail >> 8);
    code_c[1025] = (uint8_t)(tail >> 16);
    code_c[1026] = (uint8_t)(tail >> 24);
}

void decode_pk(uint8_t *code_c, int16_t *c)
{
    size_t i, j, idx;
    uint64_t tmp[147];
    uint32_t tmp32[61];
    uint16_t r86[DIM_N];
    uint16_t r3[DIM_N];

    for (i = 0; i < 147; i++) {
        tmp[i] = 0;
    }

    /*
     * Recover tmp[144], tmp[145], and tmp[146].
     *
     * Encoding layout:
     *   code[810..814] = low 40 bits of tmp[144], little-endian
     *   code[815..819] = low 40 bits of tmp[145], little-endian
     *   code[820].bits[3..7] = tmp[144].bits[40..44]
     *   code[821].bits[3..7] = tmp[145].bits[40..44]
     *   code[820].bits[0..2], code[821].bits[0..2], code[822] = tmp[146]
     */
    tmp[144] =
        ((uint64_t)code_c[810]) |
        ((uint64_t)code_c[811] << 8) |
        ((uint64_t)code_c[812] << 16) |
        ((uint64_t)code_c[813] << 24) |
        ((uint64_t)code_c[814] << 32) |
        (((uint64_t)(code_c[820] >> 3) & 31ULL) << 40);

    tmp[145] =
        ((uint64_t)code_c[815]) |
        ((uint64_t)code_c[816] << 8) |
        ((uint64_t)code_c[817] << 16) |
        ((uint64_t)code_c[818] << 24) |
        ((uint64_t)code_c[819] << 32) |
        (((uint64_t)(code_c[821] >> 3) & 31ULL) << 40);

    tmp[146] =
        ((uint64_t)code_c[820] & 7ULL) |
        (((uint64_t)code_c[821] & 7ULL) << 3) |
        ((uint64_t)code_c[822] << 6);

    /*
     * Recover high bits for tmp[0..143].
     */
    idx = 768;
    for (i = 0; i < 144; i += 24) {
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
     * Recover lower bits from the first 768 bytes.
     */
    j = 0;
    for (i = 0; i < 144; i += 3) {
        uint64_t r0 = load64_le(code_c + 8 * j);
        uint64_t r1 = load64_le(code_c + 8 * (j + 1));

        tmp[i + 0] += ((r1 & 0x7FFFFULL) << 19) | (r0 & 0x7FFFFULL);
        tmp[i + 1]  = r0 >> 19;
        tmp[i + 2]  = r1 >> 19;

        j += 2;
    }

    /*
     * Extract final two base-86 digits.
     */
    {
        uint64_t x = tmp[146];

        for (i = 0; i < 2; i++) {
            uint64_t q;
            uint16_t r;

            ct_divmod86_u64(x, &q, &r);
            r86[DIM_N - 2 + i] = r;
            x = q;
        }
    }

    /*
     * Extract base-86 digits from tmp[0..145].
     */
    idx = 0;
    for (i = 0; i < DIM_N - 2; i += 7) {
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
     * Recover packed base-3 words.
     */
    tmp32[60] =
        ((uint32_t)code_c[1023]) |
        ((uint32_t)code_c[1024] << 8) |
        ((uint32_t)code_c[1025] << 16) |
        ((uint32_t)code_c[1026] << 24);

    j = 0;
    for (i = 0; i < 60; i += 6) {
        uint32_t r0 = load32_le(code_c + 823 + 4 * j);
        uint32_t r1 = load32_le(code_c + 823 + 4 * (j + 1));
        uint32_t r2 = load32_le(code_c + 823 + 4 * (j + 2));
        uint32_t r3w = load32_le(code_c + 823 + 4 * (j + 3));
        uint32_t r4 = load32_le(code_c + 823 + 4 * (j + 4));

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

    tmp32[5]  += (((tmp32[60] >> 18) & 3U) << 25);
    tmp32[11] += (((tmp32[60] >> 16) & 3U) << 25);
    tmp32[17] += (((tmp32[60] >> 14) & 3U) << 25);
    tmp32[23] += (((tmp32[60] >> 12) & 3U) << 25);
    tmp32[29] += (((tmp32[60] >> 10) & 3U) << 25);
    tmp32[35] += (((tmp32[60] >> 8)  & 3U) << 25);
    tmp32[41] += (((tmp32[60] >> 6)  & 3U) << 25);
    tmp32[47] += (((tmp32[60] >> 4)  & 3U) << 25);
    tmp32[53] += (((tmp32[60] >> 2)  & 3U) << 25);
    tmp32[59] += (((tmp32[60]      ) & 3U) << 25);

    tmp32[60] >>= 20;

    /*
     * Extract final four base-3 digits.
     */
    {
        uint32_t x = tmp32[60];

        for (i = 0; i < 4; i++) {
            uint32_t q;
            uint16_t r;

            ct_divmod3_u32(x, &q, &r);
            r3[DIM_N - 4 + i] = r;
            x = q;
        }
    }

    /*
     * Extract base-3 digits from tmp32[0..59].
     */
    idx = 0;
    for (i = 0; i < DIM_N - 4; i += 17) {
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
     * x = r86 * 87 + r3 * 172 mod 258.
     */
    for (i = 0; i < DIM_N; i++) {
        uint32_t x = (uint32_t)r86[i] * 87U + (uint32_t)r3[i] * 172U;
        c[i] = (int16_t)ct_mod258_u32_small(x);
    }
}

/* ------------------------------------------------------------------------- */
/* c codec                                                                    */
/* ------------------------------------------------------------------------- */

void encode_c(int16_t *c, uint8_t *code_c)
{
    size_t i, index;
    uint64_t tmp[128];
    uint64_t res[104];

    /*
     * Valid input is expected to be canonical in [0, 89].
     * This matches the original codec for valid input.
     */
    index = 0;
    for (i = 0; i < DIM_N; i += 8) {
        uint64_t acc =
            (uint64_t)(uint16_t)c[i + 0] * c_encode_table_ct[0] +
            (uint64_t)(uint16_t)c[i + 1] * c_encode_table_ct[1] +
            (uint64_t)(uint16_t)c[i + 2] * c_encode_table_ct[2] +
            (uint64_t)(uint16_t)c[i + 3] * c_encode_table_ct[3] +
            (uint64_t)(uint16_t)c[i + 4] * c_encode_table_ct[4] +
            (uint64_t)(uint16_t)c[i + 5] * c_encode_table_ct[5] +
            (uint64_t)(uint16_t)c[i + 6] * c_encode_table_ct[6] +
            (uint64_t)(uint16_t)c[i + 7] * c_encode_table_ct[7];

        tmp[index++] = acc;
    }

    for (i = 0; i < 104; i++) {
        res[i] = 0;
    }

    index = 0;
    res[64] = 0;
    for (i = 0; i < 80; i += 5) {
        res[index + 0] = (tmp[i + 0] << 12) | (tmp[i + 4] & 0xFFFULL);
        res[index + 1] = (tmp[i + 1] << 12) | ((tmp[i + 4] >> 12) & 0xFFFULL);
        res[index + 2] = (tmp[i + 2] << 12) | ((tmp[i + 4] >> 24) & 0xFFFULL);
        res[index + 3] = (tmp[i + 3] << 12) | ((tmp[i + 4] >> 36) & 0xFFFULL);
        index += 4;
        res[64] = (res[64] << 4) | (tmp[i + 4] >> 48);
    }

    index = 65;
    res[101] = 0;
    for (i = 80; i < 125; i += 5) {
        res[index + 0] = (tmp[i + 0] << 12) | (tmp[i + 4] & 0xFFFULL);
        res[index + 1] = (tmp[i + 1] << 12) | ((tmp[i + 4] >> 12) & 0xFFFULL);
        res[index + 2] = (tmp[i + 2] << 12) | ((tmp[i + 4] >> 24) & 0xFFFULL);
        res[index + 3] = (tmp[i + 3] << 12) | ((tmp[i + 4] >> 36) & 0xFFFULL);
        index += 4;
        res[101] = (res[101] << 4) | (tmp[i + 4] >> 48);
    }

    res[101] = (res[101] << 28) | (tmp[125] & 0x0FFFFFFFULL);
    res[102] = (tmp[126] << 12) | ((tmp[125] >> 28) & 0xFFFULL);
    res[103] = (tmp[127] << 12) | ((tmp[125] >> 40) & 0xFFFULL);

    for (i = 0; i < 104; i++) {
        store64_le(code_c + 8 * i, res[i]);
    }
}

void decode_c(uint8_t *code_c, int16_t *c)
{
    size_t i, j, index;
    uint64_t res[104];
    uint64_t tmp[128];
    uint64_t t;

    for (i = 0; i < 104; i++) {
        res[i] = load64_le(code_c + 8 * i);
    }

    tmp[125] = ((res[103] & 0xFFFULL) << 40) |
               ((res[102] & 0xFFFULL) << 28) |
               (res[101] & 0x0FFFFFFFULL);
    tmp[126] = (res[102] >> 12) & 0xFFFFFFFFFFFFFULL;
    tmp[127] = (res[103] >> 12) & 0xFFFFFFFFFFFFFULL;

    index = 65;
    t = res[101] >> 28;
    j = 0;
    for (i = 80; i < 125; i += 5) {
        tmp[i + 0] = res[index + 0] >> 12;
        tmp[i + 1] = res[index + 1] >> 12;
        tmp[i + 2] = res[index + 2] >> 12;
        tmp[i + 3] = res[index + 3] >> 12;
        tmp[i + 4] =
            (((t >> (32 - (j * 4))) & 15ULL) << 48) |
            ((res[index + 3] & 0xFFFULL) << 36) |
            ((res[index + 2] & 0xFFFULL) << 24) |
            ((res[index + 1] & 0xFFFULL) << 12) |
            (res[index + 0] & 0xFFFULL);
        index += 4;
        j++;
    }

    index = 0;
    t = res[64];
    j = 0;
    for (i = 0; i < 80; i += 5) {
        tmp[i + 0] = res[index + 0] >> 12;
        tmp[i + 1] = res[index + 1] >> 12;
        tmp[i + 2] = res[index + 2] >> 12;
        tmp[i + 3] = res[index + 3] >> 12;
        tmp[i + 4] =
            (((t >> (60 - (j * 4))) & 15ULL) << 48) |
            ((res[index + 3] & 0xFFFULL) << 36) |
            ((res[index + 2] & 0xFFFULL) << 24) |
            ((res[index + 1] & 0xFFFULL) << 12) |
            (res[index + 0] & 0xFFFULL);
        index += 4;
        j++;
    }

    index = 0;
    for (i = 0; i < DIM_N; i += 8) {
        uint64_t x = tmp[index++];

        for (j = 0; j < 8; j++) {
            uint64_t q;
            uint16_t r;

            ct_divmod90_u64(x, &q, &r);
            c[i + j] = (int16_t)r;
            x = q;
        }
    }
}