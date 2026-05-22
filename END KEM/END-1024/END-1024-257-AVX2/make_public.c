/*******************************************************************************
 * 
 *    Original source:
 *    ntrugen by Thomas Pornin,
 *    https://github.com/pornin/ntrugen
 *
 *******************************************************************************/

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <immintrin.h>
#include "param.h"
#include "make_public.h"
#define Q2   ((uint32_t)((uint32_t)Q * (uint32_t)Q))

#define MUL_16to32(u0, u1, a, b)   do { \
		__m256i t_mul_16to32_a = (a); \
		__m256i t_mul_16to32_b = (b); \
		__m256i t_mul_16to32_lo = _mm256_mullo_epi16( \
			t_mul_16to32_a, t_mul_16to32_b); \
		__m256i t_mul_16to32_hi = _mm256_mulhi_epu16( \
			t_mul_16to32_a, t_mul_16to32_b); \
		(u0) = _mm256_unpacklo_epi16( \
			t_mul_16to32_lo, t_mul_16to32_hi); \
		(u1) = _mm256_unpackhi_epi16( \
			t_mul_16to32_lo, t_mul_16to32_hi); \
	} while (0)

#define REPACK(lo, hi, u0, u1)   do { \
		__m256i t_repack_m = _mm256_set1_epi32(0xFFFF); \
		__m256i t_repack_u0 = (u0); \
		__m256i t_repack_u1 = (u1); \
		(lo) = _mm256_packus_epi32( \
			_mm256_and_si256(t_repack_u0, t_repack_m), \
			_mm256_and_si256(t_repack_u1, t_repack_m)); \
		(hi) = _mm256_packus_epi32( \
			_mm256_srli_epi32(t_repack_u0, 16), \
			_mm256_srli_epi32(t_repack_u1, 16)); \
	} while (0)

__attribute__((aligned(32)))
static const uint16_t NX[] = {
	 365,  404,  440,  329,   87,  682,  758,   11,  727,   42,
	 297,  472,  424,  345,  627,  142,  104,  665,  473,  296,
	 305,  464,  315,  454,  224,  545,  723,   46,  302,  467,
	 501,  268,  290,  479,  476,  293,  185,  584,   65,  704,
	 388,  381,  552,  217,  221,  548,  140,  629,  504,  265,
	 281,  488,  295,  474,  166,  603,  494,  275,  132,  637,
	 103,  666,  535,  234,  156,  613,  325,  444,   73,  696,
	  88,  681,  336,  433,  700,   69,  453,  316,  367,  402,
	 706,   63,   61,  708,  636,  133,  556,  213,  515,  254,
	 368,  401,  660,  109,  606,  163,  756,   13,   37,  732,
	  58,  711,  249,  520,  741,   28,  198,  571,  539,  230,
	 418,  351,  582,  187,   59,  710,  716,   53,  210,  559,
	 662,  107,  482,  287,  714,   55,  334,  435
};

__attribute__((aligned(32)))
static const uint16_t GM[] = {
	  19,  360,  211,  760,  455,  243,  277,  513,  155,  387,
	 669,   48,  393,  242,  317,  340,  447,  739,  431,  193,
	 667,  172,   41,  534,  692,  160,  521,  765,  544,  108,
	 294,  228,  617,  196,  619,   72,  205,  363,   91,  510,
	 298,  749,   31,  385,  701,  371,  540,  356,  269,  240,
	 397,  763,   47,  162,  441,  342,  616,  258,  446,   32,
	 262,  674,  724,  483,  365,  440,   87,  758,  727,  297,
	 424,  627,  104,  473,  305,  315,  224,  723,  302,  501,
	 290,  476,  185,   65,  388,  552,  221,  140,  504,  281,
	 295,  166,  494,  132,  103,  535,  156,  325,   73,   88,
	 336,  700,  453,  367,  706,   61,  636,  556,  515,  368,
	 660,  606,  756,   37,   58,  249,  741,  198,  539,  418,
	 582,   59,  716,  210,  662,  482,  714,  334
};

__attribute__((aligned(32)))
static const union {
	uint16_t w16[64];
	__m256i w256[4];
} vGM8 = {
	{
		  155,   155,   155,   155,   155,   155,   155,   155,
		  387,   387,   387,   387,   387,   387,   387,   387,
		  669,   669,   669,   669,   669,   669,   669,   669,
		   48,    48,    48,    48,    48,    48,    48,    48,
		  393,   393,   393,   393,   393,   393,   393,   393,
		  242,   242,   242,   242,   242,   242,   242,   242,
		  317,   317,   317,   317,   317,   317,   317,   317,
		  340,   340,   340,   340,   340,   340,   340,   340
	}
};

__attribute__((aligned(32)))
static const union {
	uint16_t w16[64];
	__m256i w256[4];
} vGM16 = {
	{
		  447,   447,   447,   447,   431,   431,   431,   431,
		  739,   739,   739,   739,   193,   193,   193,   193,
		  667,   667,   667,   667,    41,    41,    41,    41,
		  172,   172,   172,   172,   534,   534,   534,   534,
		  692,   692,   692,   692,   521,   521,   521,   521,
		  160,   160,   160,   160,   765,   765,   765,   765,
		  544,   544,   544,   544,   294,   294,   294,   294,
		  108,   108,   108,   108,   228,   228,   228,   228
	}
};

__attribute__((aligned(32)))
static const union {
	uint16_t w16[64];
	__m256i w256[4];
} vGM32 = {
	{
		  617,   617,   205,   205,   196,   196,   363,   363,
		  619,   619,    91,    91,    72,    72,   510,   510,
		  298,   298,   701,   701,   749,   749,   371,   371,
		   31,    31,   540,   540,   385,   385,   356,   356,
		  269,   269,    47,    47,   240,   240,   162,   162,
		  397,   397,   441,   441,   763,   763,   342,   342,
		  616,   616,   262,   262,   258,   258,   674,   674,
		  446,   446,   724,   724,    32,    32,   483,   483
	}
};

__attribute__((aligned(32)))
static const union {
	uint16_t w16[64];
	__m256i w256[4];
} vGM64 = {
	{
		  365,   104,   440,   473,    87,   305,   758,   315,
		  727,   224,   297,   723,   424,   302,   627,   501,
		  290,   504,   476,   281,   185,   295,    65,   166,
		  388,   494,   552,   132,   221,   103,   140,   535,
		  156,   706,   325,    61,    73,   636,    88,   556,
		  336,   515,   700,   368,   453,   660,   367,   606,
		  756,   582,    37,    59,    58,   716,   249,   210,
		  741,   662,   198,   482,   539,   714,   418,   334
	}
};

__attribute__((aligned(32)))
static const uint16_t iGM[] = {
	  19,  409,    9,  558,  256,  492,  526,  314,  429,  452,
	 527,  376,  721,  100,  382,  614,  541,  475,  661,  225,
	   4,  248,  609,   77,  235,  728,  597,  102,  576,  338,
	  30,  322,  286,   45,   95,  507,  737,  323,  511,  153,
	 427,  328,  607,  722,    6,  372,  529,  500,  413,  229,
	 398,   68,  384,  738,   20,  471,  259,  678,  406,  564,
	 697,  150,  573,  152,  435,   55,  287,  107,  559,   53,
	 710,  187,  351,  230,  571,   28,  520,  711,  732,   13,
	 163,  109,  401,  254,  213,  133,  708,   63,  402,  316,
	  69,  433,  681,  696,  444,  613,  234,  666,  637,  275,
	 603,  474,  488,  265,  629,  548,  217,  381,  704,  584,
	 293,  479,  268,  467,   46,  545,  454,  464,  296,  665,
	 142,  345,  472,   42,   11,  682,  329,  404
};

__attribute__((aligned(32)))
static const union {
	uint16_t w16[64];
	__m256i w256[4];
} viGM8 = {
	{
		  429,   429,   429,   429,   429,   429,   429,   429,
		  452,   452,   452,   452,   452,   452,   452,   452,
		  527,   527,   527,   527,   527,   527,   527,   527,
		  376,   376,   376,   376,   376,   376,   376,   376,
		  721,   721,   721,   721,   721,   721,   721,   721,
		  100,   100,   100,   100,   100,   100,   100,   100,
		  382,   382,   382,   382,   382,   382,   382,   382,
		  614,   614,   614,   614,   614,   614,   614,   614
	}
};

__attribute__((aligned(32)))
static const union {
	uint16_t w16[64];
	__m256i w256[4];
} viGM16 = {
	{
		  541,   541,   541,   541,   661,   661,   661,   661,
		  475,   475,   475,   475,   225,   225,   225,   225,
		    4,     4,     4,     4,   609,   609,   609,   609,
		  248,   248,   248,   248,    77,    77,    77,    77,
		  235,   235,   235,   235,   597,   597,   597,   597,
		  728,   728,   728,   728,   102,   102,   102,   102,
		  576,   576,   576,   576,    30,    30,    30,    30,
		  338,   338,   338,   338,   322,   322,   322,   322
	}
};

__attribute__((aligned(32)))
static const union {
	uint16_t w16[64];
	__m256i w256[4];
} viGM32 = {
	{
		  286,   286,   737,   737,    45,    45,   323,   323,
		   95,    95,   511,   511,   507,   507,   153,   153,
		  427,   427,     6,     6,   328,   328,   372,   372,
		  607,   607,   529,   529,   722,   722,   500,   500,
		  413,   413,   384,   384,   229,   229,   738,   738,
		  398,   398,    20,    20,    68,    68,   471,   471,
		  259,   259,   697,   697,   678,   678,   150,   150,
		  406,   406,   573,   573,   564,   564,   152,   152
	}
};

__attribute__((aligned(32)))
static const union {
	uint16_t w16[64];
	__m256i w256[4];
} viGM64 = {
	{
		  435,   351,    55,   230,   287,   571,   107,    28,
		  559,   520,    53,   711,   710,   732,   187,    13,
		  163,   402,   109,   316,   401,    69,   254,   433,
		  213,   681,   133,   696,   708,   444,    63,   613,
		  234,   629,   666,   548,   637,   217,   275,   381,
		  603,   704,   474,   584,   488,   293,   265,   479,
		  268,   142,   467,   345,    46,   472,   545,    42,
		  454,    11,   464,   682,   296,   329,   665,   404
	}
};

static inline uint32_t mq_montyred(uint32_t x)
{
	x *= 452395775;
	x = (x >> 16) * Q;
	return (x >> 16) + 1;
}

static inline __m256i mq_montyred_x16(__m256i lo, __m256i hi)
{

	__m256i Qx16 = _mm256_set1_epi16(769);
	__m256i Q1Ilox16 = _mm256_set1_epi16(767);
	__m256i Q1Ihix16 = _mm256_set1_epi16(6903);
	__m256i x;

	/*
	 * x = (uint32_t)(x * Q1I) >> 16
	 * Each x (32 bits) is split into its low 16 bits (in lo) and
	 * high 16 bits (in hi). Q1I is a 32-bit constant. The product
	 * is computed modulo 2^32, and we are only interested in its
	 * high 16 bits.
	 */
	x = _mm256_add_epi16(
		_mm256_add_epi16(
			_mm256_mulhi_epu16(lo, Q1Ilox16),
			_mm256_mullo_epi16(lo, Q1Ihix16)),
		_mm256_mullo_epi16(hi, Q1Ilox16));

	/*
	 * x = (x * Q) >> 16
	 * x is 16 bits here; Q also fits on 16 bits.
	 */
	x = _mm256_mulhi_epu16(x, Qx16);

	/*
	 * Result is x + 1.
	 */
	return _mm256_add_epi16(x, _mm256_set1_epi16(1));
}

static inline uint32_t mq_set(int32_t x)
{
	return mq_montyred((uint32_t)((int32_t)x
		+ (int32_t)Q * (1 + ((int32_t)503109 / Q))) * 361);
}

static inline uint32_t mq_montymul(uint32_t x, uint32_t y)
{
	return mq_montyred(x * y);
}

static inline __m256i mq_montymul_x16(__m256i x, __m256i y)
{
	return mq_montyred_x16(
		_mm256_mullo_epi16(x, y),
		_mm256_mulhi_epu16(x, y));
}

static inline uint32_t mq_unorm(uint32_t x)
{
	x = mq_montyred(x);
	x &= (uint32_t)(x - Q) >> 16;
	return x;
}

static inline uint32_t mq_add(uint32_t x, uint32_t y)
{
	/* Compute -(x+y) in the -q..q-2 range. */
	x = Q - (x + y);

	/* Add q if the value is strictly negative. Note that since
	   x <= q and y <= q, a negative value will have its
	   top 16 bits all equal to 1. */
	x += Q & (x >> 16);

	/* Since we have -(x+y) in the 0..q-1 range, we can get
	   x+y = -(-(x+y)) in the 1..q range. */
	return Q - x;
}

static inline __m256i mq_add_x16(__m256i x, __m256i y)
{
	__m256i Qx16 = _mm256_set1_epi16(769);

	x = _mm256_sub_epi16(
		Qx16,
		_mm256_add_epi16(x, y));
	x = _mm256_add_epi16(
		x,
		_mm256_and_si256(
			Qx16,
			_mm256_srai_epi16(x, 15)));
	return _mm256_sub_epi16(Qx16, x);
}

static inline __m256i mul2x16(__m256i x)
{
	return _mm256_add_epi16(x, x);
}

static inline uint32_t mq_sub(uint32_t x, uint32_t y)
{
	/* Get y-x in the -q+1..q-1 range. */
	y -= x;

	/* Add q if the value is strictly negative. New range is 0..q-1 */
	y += Q & (y >> 16);

	/* Return -(y-x) = x-y. */
	return Q - y;
}

static inline __m256i mq_sub_x16(__m256i x, __m256i y)
{
	__m256i Qx16 = _mm256_set1_epi16(769);

	y = _mm256_sub_epi16(y, x);
	y = _mm256_add_epi16(
		y,
		_mm256_and_si256(
			Qx16,
			_mm256_srai_epi16(y, 15)));
	return _mm256_sub_epi16(Qx16, y);
}

static inline uint32_t mq_inv(uint32_t x)
{
	/*
	 * We use Fermat's little theorem: 1/x = x^(q-2) mod q.
	 * An efficient addition chain on the exponent is used; the
	 * chain depends on the modulus.
	 */

	uint32_t x2, x3, x5, x10, x13;

	x2 = mq_montymul(x, x);
	x3 = mq_montymul(x2, x);
	x5 = mq_montymul(x3, x2);
	x10 = mq_montymul(x5, x5);
	x13 = mq_montymul(x10, x3);
	x = mq_montymul(x13, x10);   /* x^23 */
	x = mq_montymul(x, x);       /* x^46 */
	x = mq_montymul(x, x);       /* x^92 */
	x = mq_montymul(x, x);       /* x^184 */
	x = mq_montymul(x, x);       /* x^368 */
	x = mq_montymul(x, x13);     /* x^381 */
	x = mq_montymul(x, x);       /* x^762 */
	x = mq_montymul(x, x5);      /* x^767 = 1/x */
	return x;
}

static inline __m256i mq_inv_x16(__m256i x)
{
	/*
	 * We use Fermat's little theorem: 1/x = x^(q-2) mod q.
	 * An efficient addition chain on the exponent is used; the
	 * chain depends on the modulus.
	 */
	__m256i x2, x3, x5, x10, x13;

	x2 = mq_montymul_x16(x, x);
	x3 = mq_montymul_x16(x2, x);
	x5 = mq_montymul_x16(x3, x2);
	x10 = mq_montymul_x16(x5, x5);
	x13 = mq_montymul_x16(x10, x3);
	x = mq_montymul_x16(x13, x10);   /* x^23 */
	x = mq_montymul_x16(x, x);       /* x^46 */
	x = mq_montymul_x16(x, x);       /* x^92 */
	x = mq_montymul_x16(x, x);       /* x^184 */
	x = mq_montymul_x16(x, x);       /* x^368 */
	x = mq_montymul_x16(x, x13);     /* x^381 */
	x = mq_montymul_x16(x, x);       /* x^762 */
	x = mq_montymul_x16(x, x5);      /* x^767 = 1/x */

	return x;
}

static inline __m256i negx16(__m256i x)
{
	return _mm256_sub_epi16(_mm256_set1_epi16(2 * Q), x);
}

static void unpack1024(uint16_t *d, const uint16_t *a)
{
	__m256i a0, a1, a2, a3, a4, a5, a6, a7;
	__m256i e0, e1, e2, e3, e4, e5, e6, e7;
	__m256i f0, f1, f2, f3, f4, f5, f6, f7;
	__m256i dx[64];
	size_t u;

	for (u = 0; u < 8; u ++) {
		a0 = _mm256_loadu_si256((void *)(a + (u << 7)));
		a1 = _mm256_loadu_si256((void *)(a + (u << 7) + 16));
		a2 = _mm256_loadu_si256((void *)(a + (u << 7) + 32));
		a3 = _mm256_loadu_si256((void *)(a + (u << 7) + 48));
		a4 = _mm256_loadu_si256((void *)(a + (u << 7) + 64));
		a5 = _mm256_loadu_si256((void *)(a + (u << 7) + 80));
		a6 = _mm256_loadu_si256((void *)(a + (u << 7) + 96));
		a7 = _mm256_loadu_si256((void *)(a + (u << 7) + 112));

		e0 = _mm256_permute2x128_si256(a0, a4, 0x20);
		e1 = _mm256_permute2x128_si256(a0, a4, 0x31);
		e2 = _mm256_permute2x128_si256(a1, a5, 0x20);
		e3 = _mm256_permute2x128_si256(a1, a5, 0x31);
		e4 = _mm256_permute2x128_si256(a2, a6, 0x20);
		e5 = _mm256_permute2x128_si256(a2, a6, 0x31);
		e6 = _mm256_permute2x128_si256(a3, a7, 0x20);
		e7 = _mm256_permute2x128_si256(a3, a7, 0x31);

		f0 = _mm256_unpacklo_epi16(e0, e1);
		f1 = _mm256_unpackhi_epi16(e0, e1);
		f2 = _mm256_unpacklo_epi16(e2, e3);
		f3 = _mm256_unpackhi_epi16(e2, e3);
		f4 = _mm256_unpacklo_epi16(e4, e5);
		f5 = _mm256_unpackhi_epi16(e4, e5);
		f6 = _mm256_unpacklo_epi16(e6, e7);
		f7 = _mm256_unpackhi_epi16(e6, e7);

		e0 = _mm256_unpacklo_epi32(f0, f2);
		e1 = _mm256_unpackhi_epi32(f0, f2);
		e2 = _mm256_unpacklo_epi32(f1, f3);
		e3 = _mm256_unpackhi_epi32(f1, f3);
		e4 = _mm256_unpacklo_epi32(f4, f6);
		e5 = _mm256_unpackhi_epi32(f4, f6);
		e6 = _mm256_unpacklo_epi32(f5, f7);
		e7 = _mm256_unpackhi_epi32(f5, f7);

		dx[u +  0] = _mm256_unpacklo_epi64(e0, e4);
		dx[u +  8] = _mm256_unpackhi_epi64(e0, e4);
		dx[u + 16] = _mm256_unpacklo_epi64(e1, e5);
		dx[u + 24] = _mm256_unpackhi_epi64(e1, e5);
		dx[u + 32] = _mm256_unpacklo_epi64(e2, e6);
		dx[u + 40] = _mm256_unpackhi_epi64(e2, e6);
		dx[u + 48] = _mm256_unpacklo_epi64(e3, e7);
		dx[u + 56] = _mm256_unpackhi_epi64(e3, e7);
	}

	for (u = 0; u < 64; u ++) {
		_mm256_storeu_si256((void *)(d + (u << 4)), dx[u]);
	}
}

static inline __m256i mq_montyLC2_x16(__m256i a0, __m256i b0, __m256i a1, __m256i b1)
{
	__m256i u00, u01, u10, u11, lo, hi;

	MUL_16to32(u00, u01, a0, b0);
	MUL_16to32(u10, u11, a1, b1);
	REPACK(lo, hi,
		_mm256_add_epi32(u00, u10),
		_mm256_add_epi32(u01, u11));
	return mq_montyred_x16(lo, hi);
}

static inline __m256i mq_montyLC3_x16(__m256i a0, __m256i b0, __m256i a1, __m256i b1, __m256i a2, __m256i b2)
{
	__m256i u00, u01, u10, u11, u20, u21, lo, hi;

	MUL_16to32(u00, u01, a0, b0);
	MUL_16to32(u10, u11, a1, b1);
	MUL_16to32(u20, u21, a2, b2);
	REPACK(lo, hi,
		_mm256_add_epi32(
			_mm256_add_epi32(u00, u10),
			u20),
		_mm256_add_epi32(
			_mm256_add_epi32(u01, u11),
			u21));
	return mq_montyred_x16(lo, hi);
}

static inline __m256i mq_montyLC4_x16(__m256i a0, __m256i b0, __m256i a1, __m256i b1, __m256i a2, __m256i b2, __m256i a3, __m256i b3)
{
	__m256i u00, u01, u10, u11, u20, u21, u30, u31, lo, hi;

	MUL_16to32(u00, u01, a0, b0);
	MUL_16to32(u10, u11, a1, b1);
	MUL_16to32(u20, u21, a2, b2);
	MUL_16to32(u30, u31, a3, b3);
	REPACK(lo, hi,
		_mm256_add_epi32(
			_mm256_add_epi32(u00, u10),
			_mm256_add_epi32(u20, u30)),
		_mm256_add_epi32(
			_mm256_add_epi32(u01, u11),
			_mm256_add_epi32(u21, u31)));
	return mq_montyred_x16(lo, hi);
}

static inline __m256i mq_montyLC5_x16(__m256i a0, __m256i b0, __m256i a1, __m256i b1, __m256i a2, __m256i b2, __m256i a3, __m256i b3, __m256i a4, __m256i b4)
{
	__m256i u00, u01, u10, u11, u20, u21, u30, u31, u40, u41, lo, hi;

	MUL_16to32(u00, u01, a0, b0);
	MUL_16to32(u10, u11, a1, b1);
	MUL_16to32(u20, u21, a2, b2);
	MUL_16to32(u30, u31, a3, b3);
	MUL_16to32(u40, u41, a4, b4);
	REPACK(lo, hi,
		_mm256_add_epi32(
			_mm256_add_epi32(
				_mm256_add_epi32(u00, u10),
				_mm256_add_epi32(u20, u30)),
			u40),
		_mm256_add_epi32(
			_mm256_add_epi32(
				_mm256_add_epi32(u01, u11),
				_mm256_add_epi32(u21, u31)),
			u41));
	return mq_montyred_x16(lo, hi);
}

static void repack1024(uint16_t *d, const uint16_t *a)
{
	__m256i a0, a1, a2, a3, a4, a5, a6, a7;
	__m256i e0, e1, e2, e3, e4, e5, e6, e7;
	__m256i f0, f1, f2, f3, f4, f5, f6, f7;
	__m256i dx[64];
	size_t u;

	for (u = 0; u < 8; u ++) {
		a0 = _mm256_loadu_si256((void *)(a + (u << 4)));
		a1 = _mm256_loadu_si256((void *)(a + (u << 4) + 128));
		a2 = _mm256_loadu_si256((void *)(a + (u << 4) + 256));
		a3 = _mm256_loadu_si256((void *)(a + (u << 4) + 384));
		a4 = _mm256_loadu_si256((void *)(a + (u << 4) + 512));
		a5 = _mm256_loadu_si256((void *)(a + (u << 4) + 640));
		a6 = _mm256_loadu_si256((void *)(a + (u << 4) + 768));
		a7 = _mm256_loadu_si256((void *)(a + (u << 4) + 896));

		f0 = _mm256_unpacklo_epi16(a0, a1);
		f1 = _mm256_unpackhi_epi16(a0, a1);
		f2 = _mm256_unpacklo_epi16(a2, a3);
		f3 = _mm256_unpackhi_epi16(a2, a3);
		f4 = _mm256_unpacklo_epi16(a4, a5);
		f5 = _mm256_unpackhi_epi16(a4, a5);
		f6 = _mm256_unpacklo_epi16(a6, a7);
		f7 = _mm256_unpackhi_epi16(a6, a7);

		e0 = _mm256_unpacklo_epi32(f0, f2);
		e1 = _mm256_unpackhi_epi32(f0, f2);
		e2 = _mm256_unpacklo_epi32(f1, f3);
		e3 = _mm256_unpackhi_epi32(f1, f3);
		e4 = _mm256_unpacklo_epi32(f4, f6);
		e5 = _mm256_unpackhi_epi32(f4, f6);
		e6 = _mm256_unpacklo_epi32(f5, f7);
		e7 = _mm256_unpackhi_epi32(f5, f7);

		f0 = _mm256_unpacklo_epi64(e0, e4);
		f1 = _mm256_unpackhi_epi64(e0, e4);
		f2 = _mm256_unpacklo_epi64(e1, e5);
		f3 = _mm256_unpackhi_epi64(e1, e5);
		f4 = _mm256_unpacklo_epi64(e2, e6);
		f5 = _mm256_unpackhi_epi64(e2, e6);
		f6 = _mm256_unpacklo_epi64(e3, e7);
		f7 = _mm256_unpackhi_epi64(e3, e7);

		a0 = _mm256_permute2x128_si256(f0, f1, 0x20);
		a4 = _mm256_permute2x128_si256(f0, f1, 0x31);
		a1 = _mm256_permute2x128_si256(f2, f3, 0x20);
		a5 = _mm256_permute2x128_si256(f2, f3, 0x31);
		a2 = _mm256_permute2x128_si256(f4, f5, 0x20);
		a6 = _mm256_permute2x128_si256(f4, f5, 0x31);
		a3 = _mm256_permute2x128_si256(f6, f7, 0x20);
		a7 = _mm256_permute2x128_si256(f6, f7, 0x31);

		dx[(u << 3) + 0] = a0;
		dx[(u << 3) + 1] = a1;
		dx[(u << 3) + 2] = a2;
		dx[(u << 3) + 3] = a3;
		dx[(u << 3) + 4] = a4;
		dx[(u << 3) + 5] = a5;
		dx[(u << 3) + 6] = a6;
		dx[(u << 3) + 7] = a7;
	}

	for (u = 0; u < 64; u ++) {
		_mm256_storeu_si256((void *)(d + (u << 4)), dx[u]);
	}
}

static void NTT128(uint16_t *d, const uint16_t *a)
{
	__m256i s, u, v, e, f, pp;
	__m256i d0, d1, d2, d3, d4, d5, d6, d7;

	/* m = 1, t = 128 */
	s = _mm256_set1_epi16((short)GM[1]);
	u = _mm256_loadu_si256((void *)a);
	v = mq_montymul_x16(s, _mm256_loadu_si256((void *)(a + 64)));
	d0 = mq_add_x16(u, v);
	d4 = mq_sub_x16(u, v);
	u = _mm256_loadu_si256((void *)(a + 16));
	v = mq_montymul_x16(s, _mm256_loadu_si256((void *)(a + 80)));
	d1 = mq_add_x16(u, v);
	d5 = mq_sub_x16(u, v);
	u = _mm256_loadu_si256((void *)(a + 32));
	v = mq_montymul_x16(s, _mm256_loadu_si256((void *)(a + 96)));
	d2 = mq_add_x16(u, v);
	d6 = mq_sub_x16(u, v);
	u = _mm256_loadu_si256((void *)(a + 48));
	v = mq_montymul_x16(s, _mm256_loadu_si256((void *)(a + 112)));
	d3 = mq_add_x16(u, v);
	d7 = mq_sub_x16(u, v);

	/* m = 2, t = 64 */
	s = _mm256_set1_epi16((short)GM[2]);
	u = d0;
	v = mq_montymul_x16(s, d2);
	d0 = mq_add_x16(u, v);
	d2 = mq_sub_x16(u, v);
	u = d1;
	v = mq_montymul_x16(s, d3);
	d1 = mq_add_x16(u, v);
	d3 = mq_sub_x16(u, v);
	s = _mm256_set1_epi16((short)GM[3]);
	u = d4;
	v = mq_montymul_x16(s, d6);
	d4 = mq_add_x16(u, v);
	d6 = mq_sub_x16(u, v);
	u = d5;
	v = mq_montymul_x16(s, d7);
	d5 = mq_add_x16(u, v);
	d7 = mq_sub_x16(u, v);

	/* m = 4, t = 32 */
	s = _mm256_set1_epi16((short)GM[4]);
	u = d0;
	v = mq_montymul_x16(s, d1);
	d0 = mq_add_x16(u, v);
	d1 = mq_sub_x16(u, v);
	s = _mm256_set1_epi16((short)GM[5]);
	u = d2;
	v = mq_montymul_x16(s, d3);
	d2 = mq_add_x16(u, v);
	d3 = mq_sub_x16(u, v);
	s = _mm256_set1_epi16((short)GM[6]);
	u = d4;
	v = mq_montymul_x16(s, d5);
	d4 = mq_add_x16(u, v);
	d5 = mq_sub_x16(u, v);
	s = _mm256_set1_epi16((short)GM[7]);
	u = d6;
	v = mq_montymul_x16(s, d7);
	d6 = mq_add_x16(u, v);
	d7 = mq_sub_x16(u, v);

	/* m = 8, t = 16 */
	/*
	 * For this step, the two lanes of each 256-bit register should
	 * be combined together.
	 */
	u = _mm256_permute2x128_si256(d0, d1, 0x20);
	v = mq_montymul_x16(
		_mm256_permute2x128_si256(d0, d1, 0x31), vGM8.w256[0]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	d0 = _mm256_permute2x128_si256(e, f, 0x20);
	d1 = _mm256_permute2x128_si256(e, f, 0x31);
	u = _mm256_permute2x128_si256(d2, d3, 0x20);
	v = mq_montymul_x16(
		_mm256_permute2x128_si256(d2, d3, 0x31), vGM8.w256[1]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	d2 = _mm256_permute2x128_si256(e, f, 0x20);
	d3 = _mm256_permute2x128_si256(e, f, 0x31);
	u = _mm256_permute2x128_si256(d4, d5, 0x20);
	v = mq_montymul_x16(
		_mm256_permute2x128_si256(d4, d5, 0x31), vGM8.w256[2]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	d4 = _mm256_permute2x128_si256(e, f, 0x20);
	d5 = _mm256_permute2x128_si256(e, f, 0x31);
	u = _mm256_permute2x128_si256(d6, d7, 0x20);
	v = mq_montymul_x16(
		_mm256_permute2x128_si256(d6, d7, 0x31), vGM8.w256[3]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	d6 = _mm256_permute2x128_si256(e, f, 0x20);
	d7 = _mm256_permute2x128_si256(e, f, 0x31);

	/* m = 16, t = 8 */
	u = _mm256_unpacklo_epi64(d0, d1);
	v = mq_montymul_x16(_mm256_unpackhi_epi64(d0, d1), vGM16.w256[0]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	d0 = _mm256_unpacklo_epi64(e, f);
	d1 = _mm256_unpackhi_epi64(e, f);
	u = _mm256_unpacklo_epi64(d2, d3);
	v = mq_montymul_x16(_mm256_unpackhi_epi64(d2, d3), vGM16.w256[1]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	d2 = _mm256_unpacklo_epi64(e, f);
	d3 = _mm256_unpackhi_epi64(e, f);
	u = _mm256_unpacklo_epi64(d4, d5);
	v = mq_montymul_x16(_mm256_unpackhi_epi64(d4, d5), vGM16.w256[2]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	d4 = _mm256_unpacklo_epi64(e, f);
	d5 = _mm256_unpackhi_epi64(e, f);
	u = _mm256_unpacklo_epi64(d6, d7);
	v = mq_montymul_x16(_mm256_unpackhi_epi64(d6, d7), vGM16.w256[3]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	d6 = _mm256_unpacklo_epi64(e, f);
	d7 = _mm256_unpackhi_epi64(e, f);

	/* m = 32, t = 4 */
	e = _mm256_shuffle_epi32(d0, 0xD8);
	f = _mm256_shuffle_epi32(d1, 0xD8);
	u = _mm256_unpacklo_epi32(e, f);
	v = mq_montymul_x16(_mm256_unpackhi_epi32(e, f), vGM32.w256[0]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	e = _mm256_shuffle_epi32(e, 0xD8);
	f = _mm256_shuffle_epi32(f, 0xD8);
	d0 = _mm256_unpacklo_epi32(e, f);
	d1 = _mm256_unpackhi_epi32(e, f);

	e = _mm256_shuffle_epi32(d2, 0xD8);
	f = _mm256_shuffle_epi32(d3, 0xD8);
	u = _mm256_unpacklo_epi32(e, f);
	v = mq_montymul_x16(_mm256_unpackhi_epi32(e, f), vGM32.w256[1]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	e = _mm256_shuffle_epi32(e, 0xD8);
	f = _mm256_shuffle_epi32(f, 0xD8);
	d2 = _mm256_unpacklo_epi32(e, f);
	d3 = _mm256_unpackhi_epi32(e, f);

	e = _mm256_shuffle_epi32(d4, 0xD8);
	f = _mm256_shuffle_epi32(d5, 0xD8);
	u = _mm256_unpacklo_epi32(e, f);
	v = mq_montymul_x16(_mm256_unpackhi_epi32(e, f), vGM32.w256[2]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	e = _mm256_shuffle_epi32(e, 0xD8);
	f = _mm256_shuffle_epi32(f, 0xD8);
	d4 = _mm256_unpacklo_epi32(e, f);
	d5 = _mm256_unpackhi_epi32(e, f);

	e = _mm256_shuffle_epi32(d6, 0xD8);
	f = _mm256_shuffle_epi32(d7, 0xD8);
	u = _mm256_unpacklo_epi32(e, f);
	v = mq_montymul_x16(_mm256_unpackhi_epi32(e, f), vGM32.w256[3]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	e = _mm256_shuffle_epi32(e, 0xD8);
	f = _mm256_shuffle_epi32(f, 0xD8);
	d6 = _mm256_unpacklo_epi32(e, f);
	d7 = _mm256_unpackhi_epi32(e, f);

	/* m = 64, t = 2 */
	pp = _mm256_setr_epi8(
		0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15,
		0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15);

	e = _mm256_shuffle_epi8(d0, pp);
	f = _mm256_shuffle_epi8(d1, pp);
	u = _mm256_unpacklo_epi16(e, f);
	v = mq_montymul_x16(_mm256_unpackhi_epi16(e, f), vGM64.w256[0]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	e = _mm256_shuffle_epi8(e, pp);
	f = _mm256_shuffle_epi8(f, pp);
	d0 = _mm256_unpacklo_epi16(e, f);
	d1 = _mm256_unpackhi_epi16(e, f);

	e = _mm256_shuffle_epi8(d2, pp);
	f = _mm256_shuffle_epi8(d3, pp);
	u = _mm256_unpacklo_epi16(e, f);
	v = mq_montymul_x16(_mm256_unpackhi_epi16(e, f), vGM64.w256[1]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	e = _mm256_shuffle_epi8(e, pp);
	f = _mm256_shuffle_epi8(f, pp);
	d2 = _mm256_unpacklo_epi16(e, f);
	d3 = _mm256_unpackhi_epi16(e, f);

	e = _mm256_shuffle_epi8(d4, pp);
	f = _mm256_shuffle_epi8(d5, pp);
	u = _mm256_unpacklo_epi16(e, f);
	v = mq_montymul_x16(_mm256_unpackhi_epi16(e, f), vGM64.w256[2]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	e = _mm256_shuffle_epi8(e, pp);
	f = _mm256_shuffle_epi8(f, pp);
	d4 = _mm256_unpacklo_epi16(e, f);
	d5 = _mm256_unpackhi_epi16(e, f);

	e = _mm256_shuffle_epi8(d6, pp);
	f = _mm256_shuffle_epi8(d7, pp);
	u = _mm256_unpacklo_epi16(e, f);
	v = mq_montymul_x16(_mm256_unpackhi_epi16(e, f), vGM64.w256[3]);
	e = mq_add_x16(u, v);
	f = mq_sub_x16(u, v);
	e = _mm256_shuffle_epi8(e, pp);
	f = _mm256_shuffle_epi8(f, pp);
	d6 = _mm256_unpacklo_epi16(e, f);
	d7 = _mm256_unpackhi_epi16(e, f);

	_mm256_storeu_si256((void *)d, d0);
	_mm256_storeu_si256((void *)(d + 16), d1);
	_mm256_storeu_si256((void *)(d + 32), d2);
	_mm256_storeu_si256((void *)(d + 48), d3);
	_mm256_storeu_si256((void *)(d + 64), d4);
	_mm256_storeu_si256((void *)(d + 80), d5);
	_mm256_storeu_si256((void *)(d + 96), d6);
	_mm256_storeu_si256((void *)(d + 112), d7);
}

static void iNTT128(uint16_t *d, const uint16_t *a)
{
	__m256i s, s0, s1, u, v, e, f, pp;
	__m256i d0, d1, d2, d3, d4, d5, d6, d7;

	d0 = _mm256_loadu_si256((void *)a);
	d1 = _mm256_loadu_si256((void *)(a + 16));
	d2 = _mm256_loadu_si256((void *)(a + 32));
	d3 = _mm256_loadu_si256((void *)(a + 48));
	d4 = _mm256_loadu_si256((void *)(a + 64));
	d5 = _mm256_loadu_si256((void *)(a + 80));
	d6 = _mm256_loadu_si256((void *)(a + 96));
	d7 = _mm256_loadu_si256((void *)(a + 112));

	/* m = 64, t = 2 */
	pp = _mm256_setr_epi8(
		0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15,
		0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15);

	e = _mm256_shuffle_epi8(d0, pp);
	f = _mm256_shuffle_epi8(d1, pp);
	u = _mm256_unpacklo_epi16(e, f);
	v = _mm256_unpackhi_epi16(e, f);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM64.w256[0]);
	e = _mm256_shuffle_epi8(e, pp);
	f = _mm256_shuffle_epi8(f, pp);
	d0 = _mm256_unpacklo_epi16(e, f);
	d1 = _mm256_unpackhi_epi16(e, f);

	e = _mm256_shuffle_epi8(d2, pp);
	f = _mm256_shuffle_epi8(d3, pp);
	u = _mm256_unpacklo_epi16(e, f);
	v = _mm256_unpackhi_epi16(e, f);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM64.w256[1]);
	e = _mm256_shuffle_epi8(e, pp);
	f = _mm256_shuffle_epi8(f, pp);
	d2 = _mm256_unpacklo_epi16(e, f);
	d3 = _mm256_unpackhi_epi16(e, f);

	e = _mm256_shuffle_epi8(d4, pp);
	f = _mm256_shuffle_epi8(d5, pp);
	u = _mm256_unpacklo_epi16(e, f);
	v = _mm256_unpackhi_epi16(e, f);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM64.w256[2]);
	e = _mm256_shuffle_epi8(e, pp);
	f = _mm256_shuffle_epi8(f, pp);
	d4 = _mm256_unpacklo_epi16(e, f);
	d5 = _mm256_unpackhi_epi16(e, f);

	e = _mm256_shuffle_epi8(d6, pp);
	f = _mm256_shuffle_epi8(d7, pp);
	u = _mm256_unpacklo_epi16(e, f);
	v = _mm256_unpackhi_epi16(e, f);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM64.w256[3]);
	e = _mm256_shuffle_epi8(e, pp);
	f = _mm256_shuffle_epi8(f, pp);
	d6 = _mm256_unpacklo_epi16(e, f);
	d7 = _mm256_unpackhi_epi16(e, f);

	/* m = 32, t = 4 */
	e = _mm256_shuffle_epi32(d0, 0xD8);
	f = _mm256_shuffle_epi32(d1, 0xD8);
	u = _mm256_unpacklo_epi32(e, f);
	v = _mm256_unpackhi_epi32(e, f);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM32.w256[0]);
	e = _mm256_shuffle_epi32(e, 0xD8);
	f = _mm256_shuffle_epi32(f, 0xD8);
	d0 = _mm256_unpacklo_epi32(e, f);
	d1 = _mm256_unpackhi_epi32(e, f);

	e = _mm256_shuffle_epi32(d2, 0xD8);
	f = _mm256_shuffle_epi32(d3, 0xD8);
	u = _mm256_unpacklo_epi32(e, f);
	v = _mm256_unpackhi_epi32(e, f);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM32.w256[1]);
	e = _mm256_shuffle_epi32(e, 0xD8);
	f = _mm256_shuffle_epi32(f, 0xD8);
	d2 = _mm256_unpacklo_epi32(e, f);
	d3 = _mm256_unpackhi_epi32(e, f);

	e = _mm256_shuffle_epi32(d4, 0xD8);
	f = _mm256_shuffle_epi32(d5, 0xD8);
	u = _mm256_unpacklo_epi32(e, f);
	v = _mm256_unpackhi_epi32(e, f);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM32.w256[2]);
	e = _mm256_shuffle_epi32(e, 0xD8);
	f = _mm256_shuffle_epi32(f, 0xD8);
	d4 = _mm256_unpacklo_epi32(e, f);
	d5 = _mm256_unpackhi_epi32(e, f);

	e = _mm256_shuffle_epi32(d6, 0xD8);
	f = _mm256_shuffle_epi32(d7, 0xD8);
	u = _mm256_unpacklo_epi32(e, f);
	v = _mm256_unpackhi_epi32(e, f);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM32.w256[3]);
	e = _mm256_shuffle_epi32(e, 0xD8);
	f = _mm256_shuffle_epi32(f, 0xD8);
	d6 = _mm256_unpacklo_epi32(e, f);
	d7 = _mm256_unpackhi_epi32(e, f);

	/* m = 16, t = 8 */
	u = _mm256_unpacklo_epi64(d0, d1);
	v = _mm256_unpackhi_epi64(d0, d1);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM16.w256[0]);
	d0 = _mm256_unpacklo_epi64(e, f);
	d1 = _mm256_unpackhi_epi64(e, f);
	u = _mm256_unpacklo_epi64(d2, d3);
	v = _mm256_unpackhi_epi64(d2, d3);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM16.w256[1]);
	d2 = _mm256_unpacklo_epi64(e, f);
	d3 = _mm256_unpackhi_epi64(e, f);
	u = _mm256_unpacklo_epi64(d4, d5);
	v = _mm256_unpackhi_epi64(d4, d5);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM16.w256[2]);
	d4 = _mm256_unpacklo_epi64(e, f);
	d5 = _mm256_unpackhi_epi64(e, f);
	u = _mm256_unpacklo_epi64(d6, d7);
	v = _mm256_unpackhi_epi64(d6, d7);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM16.w256[3]);
	d6 = _mm256_unpacklo_epi64(e, f);
	d7 = _mm256_unpackhi_epi64(e, f);

	/* m = 8, t = 16 */
	/*
	 * For this step, the two lanes of each 256-bit register should
	 * be combined together.
	 */
	u = _mm256_permute2x128_si256(d0, d1, 0x20);
	v = _mm256_permute2x128_si256(d0, d1, 0x31);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM8.w256[0]);
	d0 = _mm256_permute2x128_si256(e, f, 0x20);
	d1 = _mm256_permute2x128_si256(e, f, 0x31);
	u = _mm256_permute2x128_si256(d2, d3, 0x20);
	v = _mm256_permute2x128_si256(d2, d3, 0x31);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM8.w256[1]);
	d2 = _mm256_permute2x128_si256(e, f, 0x20);
	d3 = _mm256_permute2x128_si256(e, f, 0x31);
	u = _mm256_permute2x128_si256(d4, d5, 0x20);
	v = _mm256_permute2x128_si256(d4, d5, 0x31);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM8.w256[2]);
	d4 = _mm256_permute2x128_si256(e, f, 0x20);
	d5 = _mm256_permute2x128_si256(e, f, 0x31);
	u = _mm256_permute2x128_si256(d6, d7, 0x20);
	v = _mm256_permute2x128_si256(d6, d7, 0x31);
	e = mq_add_x16(u, v);
	f = mq_montymul_x16(mq_sub_x16(u, v), viGM8.w256[3]);
	d6 = _mm256_permute2x128_si256(e, f, 0x20);
	d7 = _mm256_permute2x128_si256(e, f, 0x31);

	/* m = 4, t = 32 */
	s = _mm256_set1_epi16((short)iGM[4]);
	u = d0;
	v = d1;
	d0 = mq_add_x16(u, v);
	d1 = mq_montymul_x16(mq_sub_x16(u, v), s);
	s = _mm256_set1_epi16((short)iGM[5]);
	u = d2;
	v = d3;
	d2 = mq_add_x16(u, v);
	d3 = mq_montymul_x16(mq_sub_x16(u, v), s);
	s = _mm256_set1_epi16((short)iGM[6]);
	u = d4;
	v = d5;
	d4 = mq_add_x16(u, v);
	d5 = mq_montymul_x16(mq_sub_x16(u, v), s);
	s = _mm256_set1_epi16((short)iGM[7]);
	u = d6;
	v = d7;
	d6 = mq_add_x16(u, v);
	d7 = mq_montymul_x16(mq_sub_x16(u, v), s);

	/* m = 2, t = 64 */
	s = _mm256_set1_epi16((short)iGM[2]);
	u = d0;
	v = d2;
	d0 = mq_add_x16(u, v);
	d2 = mq_montymul_x16(mq_sub_x16(u, v), s);
	u = d1;
	v = d3;
	d1 = mq_add_x16(u, v);
	d3 = mq_montymul_x16(mq_sub_x16(u, v), s);
	s = _mm256_set1_epi16((short)iGM[3]);
	u = d4;
	v = d6;
	d4 = mq_add_x16(u, v);
	d6 = mq_montymul_x16(mq_sub_x16(u, v), s);
	u = d5;
	v = d7;
	d5 = mq_add_x16(u, v);
	d7 = mq_montymul_x16(mq_sub_x16(u, v), s);

	/* m = 1, t = 128 */
	/*
	 * We integrate the final division by 128 into the multipliers.
	 */
	s0 = _mm256_set1_epi16((short)655);
	s1 = _mm256_set1_epi16((short)622);
	u = d0;
	v = d4;
	d0 = mq_montymul_x16(mq_add_x16(u, v), s0);
	d4 = mq_montymul_x16(mq_sub_x16(u, v), s1);
	u = d1;
	v = d5;
	d1 = mq_montymul_x16(mq_add_x16(u, v), s0);
	d5 = mq_montymul_x16(mq_sub_x16(u, v), s1);
	u = d2;
	v = d6;
	d2 = mq_montymul_x16(mq_add_x16(u, v), s0);
	d6 = mq_montymul_x16(mq_sub_x16(u, v), s1);
	u = d3;
	v = d7;
	d3 = mq_montymul_x16(mq_add_x16(u, v), s0);
	d7 = mq_montymul_x16(mq_sub_x16(u, v), s1);

	_mm256_storeu_si256((void *)d, d0);
	_mm256_storeu_si256((void *)(d + 16), d1);
	_mm256_storeu_si256((void *)(d + 32), d2);
	_mm256_storeu_si256((void *)(d + 48), d3);
	_mm256_storeu_si256((void *)(d + 64), d4);
	_mm256_storeu_si256((void *)(d + 80), d5);
	_mm256_storeu_si256((void *)(d + 96), d6);
	_mm256_storeu_si256((void *)(d + 112), d7);
}

static void NTT(uint16_t *d, const uint16_t *a, unsigned logn)
{

	if (logn < 7) {
		unsigned n, t, m;

		n = 1u << logn;
		if (d != a) {
			memmove(d, a, n * sizeof *a);
		}
		t = n;
		for (m = 1; m < n; m <<= 1) {
			unsigned ht, i, j1;

			ht = t >> 1;
			for (i = 0, j1 = 0; i < m; i ++, j1 += t) {
				unsigned j, j2;
				uint32_t s;

				s = GM[m + i];
				j2 = j1 + ht;
				for (j = j1; j < j2; j ++) {
					uint32_t u, v;

					u = d[j];
					v = mq_montymul(d[j + ht], s);
					d[j] = mq_add(u, v);
					d[j + ht] = mq_sub(u, v);
				}
			}
			t = ht;
		}
		return;
	}

	unpack1024(d, a);
	NTT128(d, d);
	NTT128(d + 128, d + 128);
	NTT128(d + 256, d + 256);
	NTT128(d + 384, d + 384);
	NTT128(d + 512, d + 512);
	NTT128(d + 640, d + 640);
	NTT128(d + 768, d + 768);
	NTT128(d + 896, d + 896);
	return;

}

static int mq_poly_inv_ntt(uint16_t *d, const uint16_t *a, unsigned logn)
{
	size_t u, n;
	uint32_t z;
	__m256i Qx16, zz;

	if (logn <= 3) {
		z = (uint32_t)-1;
		n = (size_t)1 << logn;
		for (u = 0; u < n; u ++) {
			z &= a[u] - Q;
			d[u] = mq_inv(a[u]);
		}
		return (int)(z >> 31);
	}

	Qx16 = _mm256_set1_epi16(769);
	if (logn <= 7) {
		n = (size_t)1 << logn;
		zz = _mm256_setzero_si256();
		for (u = 0; u < n; u += 16) {
			__m256i a0;

			a0 = _mm256_loadu_si256((void *)(a + u));
			zz = _mm256_or_si256(zz, _mm256_cmpeq_epi16(a0, Qx16));
			_mm256_storeu_si256((void *)(d + u), mq_inv_x16(a0));
		}
		zz = _mm256_or_si256(zz, _mm256_bsrli_epi128(zz, 8));
		zz = _mm256_or_si256(zz, _mm256_bsrli_epi128(zz, 4));
		zz = _mm256_or_si256(zz, _mm256_bsrli_epi128(zz, 2));
		z = (uint32_t)_mm256_extract_epi32(zz, 0)
			| (uint32_t)_mm256_extract_epi32(zz, 4);
		return 1 - (int)(z & 1);
	}

	/*
	 * For larger degrees, we split polynomial a[] into its odd and
	 * even coefficients:
	 *    a = a0(X^2) + X*a1(X^2)
	 * With a0 and a1 being half-degree polynomials (they operate
	 * modulo X^(n/2)+1).
	 *
	 * We then define an adjoint:
	 *    a' = a0(X^2) - X*a1(X^2)
	 * This yields:
	 *    a*a' = (a0^2)(X^2) - X^2*(a1^2)(X^2)
	 *         = (a0^2 - X*a1^2)(X^2)
	 * i.e. a*a' is a half-degree polynomial (composed with X^2).
	 *
	 * If we can invert a*a', then:
	 *    1/a = (1/(a*a')) * a'
	 * It can be shown that a*a' is invertible if and only if a is
	 * invertible.
	 *
	 * Thus, to invert a polynomial modulo X^n+1, we just have to
	 * invert another polynomial modulo X^(n/2)+1. We can apply this
	 * process recursively to get down to degree 128 (logn = 7),
	 * that we can handle with coefficient-wise inversion in NTT
	 * representation.
	 */
	zz = _mm256_set1_epi16(0);

	for (u = 0; u < 128; u += 16) {
		__m256i a0, a1, a2, a3, a4, a5, a6, a7;
		__m256i b0, b1, b2, b3, c0, c1, e, x, nx;
		__m256i f0, f1, f2, f3;

		a0 = _mm256_loadu_si256((void *)(a + u));
		a1 = _mm256_loadu_si256((void *)(a + u + 128));
		a2 = _mm256_loadu_si256((void *)(a + u + 256));
		a3 = _mm256_loadu_si256((void *)(a + u + 384));
		a4 = _mm256_loadu_si256((void *)(a + u + 512));
		a5 = _mm256_loadu_si256((void *)(a + u + 640));
		a6 = _mm256_loadu_si256((void *)(a + u + 768));
		a7 = _mm256_loadu_si256((void *)(a + u + 896));
		x = _mm256_loadu_si256((void *)(NX + u));
		nx = negx16(x);

		b0 = mq_montyLC2_x16(
			a0, a0,
			x, mq_montyLC4_x16(
				a4, a4,
				mul2x16(a2), a6,
				mul2x16(a1), negx16(a7),
				mul2x16(a3), negx16(a5)));
		b1 = mq_montyLC3_x16(
			mul2x16(a0), a2,
			a1, negx16(a1),
			x, mq_montyLC3_x16(
				a5, negx16(a5),
				mul2x16(a4), a6,
				mul2x16(a3), negx16(a7)));
		b2 = mq_montyLC4_x16(
			a2, a2,
			mul2x16(a0), a4,
			mul2x16(a1), negx16(a3),
			x, mq_montyLC2_x16(
				a6, a6,
				mul2x16(a5), negx16(a7)));
		b3 = mq_montyLC5_x16(
			a3, negx16(a3),
			mul2x16(a0), a6,
			mul2x16(a2), a4,
			mul2x16(a1), negx16(a5),
			nx, mq_montymul_x16(a7, a7));

		/*
		b0 = mq_montyred(a0 * a0 + x * mq_montyred(
			4 * Q2 + a4 * a4
			+ 2 * (a2 * a6 - a1 * a7 - a3 * a5)));
		b1 = mq_montyred(Q2 + 2 * a0 * a2 - a1 * a1
			+ x * mq_montyred(3 * Q2 - a5 * a5
			+ 2 * (a4 * a6 - a3 * a7)));
		b2 = mq_montyred(2 * Q2 + a2 * a2
			+ 2 * (a0 * a4 - a1 * a3)
			+ x * mq_montyred(2 * Q2
			+ a6 * a6 - 2 * a5 * a7));
		b3 = mq_montyred(4 * Q2 - a3 * a3
			+ 2 * (a0 * a6 + a2 * a4 - a1 * a5)
			- x * mq_montyred(a7 * a7));
		*/

		c0 = mq_montyLC2_x16(
			b0, b0,
			x, mq_montyLC2_x16(
				b2, b2,
				mul2x16(b1), negx16(b3)));
		c1 = mq_montyLC3_x16(
			mul2x16(b0), b2,
			b1, negx16(b1),
			nx, mq_montymul_x16(b3, b3));

		/*
		c0 = mq_montyred(b0 * b0 + x * mq_montyred(
			2 * Q2 + b2 * b2 - 2 * b1 * b3));
		c1 = mq_montyred(2 * Q2
			+ 2 * b0 * b2 - b1 * b1
			- x * mq_montyred(b3 * b3));
		*/

		e = mq_inv_x16(mq_montyLC2_x16(
			c0, c0, nx, mq_montymul_x16(c1, c1)));
		zz = _mm256_or_si256(zz, _mm256_cmpeq_epi16(e, Qx16));

		/*
		e = mq_inv(mq_montyred(
			Q2 + c0 * c0 - x * mq_montyred(c1 * c1)));
		z &= e - Q;
		*/

		c0 = mq_montymul_x16(c0, e);
		c1 = mq_montymul_x16(c1, negx16(e));

		/*
		c0 = mq_montyred(c0 * e);
		c1 = mq_montyred(c1 * (2 * Q - e));
		*/

		f0 = mq_montyLC2_x16(
			b0, c0,
			x, mq_montymul_x16(b2, c1));
		f1 = mq_montyLC2_x16(
			b1, negx16(c0),
			nx, mq_montymul_x16(b3, c1));
		f2 = mq_montyLC2_x16(
			b2, c0,
			b0, c1);
		f3 = mq_montyLC2_x16(
			b3, negx16(c0),
			b1, negx16(c1));

		/*
		f0 = mq_montyred(b0 * c0 + x * mq_montyred(b2 * c1));
		f1 = mq_montyred(3 * Q2
			- b1 * c0 - x * mq_montyred(b3 * c1));
		f2 = mq_montyred(b2 * c0 + b0 * c1);
		f3 = mq_montyred(3 * Q2 - b3 * c0 - b1 * c1);
		*/

		_mm256_storeu_si256((void *)(d + u),
			mq_montyLC2_x16(
				a0, f0,
				x, mq_montyLC3_x16(
					a2, f3,
					a4, f2,
					a6, f1)));
		_mm256_storeu_si256((void *)(d + u + 128),
			mq_montyLC2_x16(
				a1, negx16(f0),
				nx, mq_montyLC3_x16(
					a3, f3,
					a5, f2,
					a7, f1)));
		_mm256_storeu_si256((void *)(d + u + 256),
			mq_montyLC3_x16(
				a0, f1,
				a2, f0,
				x, mq_montyLC2_x16(
					a4, f3,
					a6, f2)));
		_mm256_storeu_si256((void *)(d + u + 384),
			mq_montyLC3_x16(
				a1, negx16(f1),
				a3, negx16(f0),
				nx, mq_montyLC2_x16(
					a5, f3,
					a7, f2)));

		/*
		d[u] = mq_montyred(a0 * f0 + x * mq_montyred(
			a2 * f3 + a4 * f2 + a6 * f1));
		d[u + 128] = mq_montyred(3 * Q2 - a1 * f0
			- x * mq_montyred(a3 * f3 + a5 * f2 + a7 * f1));
		d[u + 256] = mq_montyred(a0 * f1 + a2 * f0
			+ x * mq_montyred(a4 * f3 + a6 * f2));
		d[u + 384] = mq_montyred(4 * Q2 - a1 * f1 - a3 * f0
			- x * mq_montyred(a5 * f3 + a7 * f2));
		*/

		_mm256_storeu_si256((void *)(d + u + 512),
			mq_montyLC4_x16(
				a0, f2,
				a2, f1,
				a4, f0,
				x, mq_montymul_x16(a6, f3)));
		_mm256_storeu_si256((void *)(d + u + 640),
			mq_montyLC4_x16(
				a1, negx16(f2),
				a3, negx16(f1),
				a5, negx16(f0),
				nx, mq_montymul_x16(a7, f3)));
		_mm256_storeu_si256((void *)(d + u + 768),
			mq_montyLC4_x16(
				a0, f3,
				a2, f2,
				a4, f1,
				a6, f0));
		_mm256_storeu_si256((void *)(d + u + 896),
			mq_montyLC4_x16(
				a1, negx16(f3),
				a3, negx16(f2),
				a5, negx16(f1),
				a7, negx16(f0)));

		/*
		d[u + 512] = mq_montyred(a0 * f2 + a2 * f1
			+ a4 * f0 + x * mq_montyred(a6 * f3));
		d[u + 640] = mq_montyred(5 * Q2 - a1 * f2 - a3 * f1
			- a5 * f0 - x * mq_montyred(a7 * f3));
		d[u + 768] = mq_montyred(a0 * f3 + a2 * f2
			+ a4 * f1 + a6 * f0);
		d[u + 896] = mq_montyred(5 * Q2 - a1 * f3 - a3 * f2
			- a5 * f1 - a7 * f0);
		*/
	}
	zz = _mm256_or_si256(zz, _mm256_bsrli_epi128(zz, 8));
	zz = _mm256_or_si256(zz, _mm256_bsrli_epi128(zz, 4));
	zz = _mm256_or_si256(zz, _mm256_bsrli_epi128(zz, 2));
	z = (uint32_t)_mm256_extract_epi32(zz, 0)
		| (uint32_t)_mm256_extract_epi32(zz, 4);
	return 1 - (int)(z & 1);
}

static inline __m256i mq_montyLC6_x16(__m256i a0, __m256i b0, __m256i a1, __m256i b1, __m256i a2, __m256i b2, __m256i a3, __m256i b3, __m256i a4, __m256i b4, __m256i a5, __m256i b5)
{
	__m256i u00, u01, u10, u11, u20, u21, u30, u31;
	__m256i u40, u41, u50, u51, lo, hi;

	MUL_16to32(u00, u01, a0, b0);
	MUL_16to32(u10, u11, a1, b1);
	MUL_16to32(u20, u21, a2, b2);
	MUL_16to32(u30, u31, a3, b3);
	MUL_16to32(u40, u41, a4, b4);
	MUL_16to32(u50, u51, a5, b5);
	REPACK(lo, hi,
		_mm256_add_epi32(
			_mm256_add_epi32(
				_mm256_add_epi32(u00, u10),
				_mm256_add_epi32(u20, u30)),
			_mm256_add_epi32(u40, u50)),
		_mm256_add_epi32(
			_mm256_add_epi32(
				_mm256_add_epi32(u01, u11),
				_mm256_add_epi32(u21, u31)),
			_mm256_add_epi32(u41, u51)));
	return mq_montyred_x16(lo, hi);
}

static inline __m256i mq_montyLC7_x16(__m256i a0, __m256i b0, __m256i a1, __m256i b1, __m256i a2, __m256i b2, __m256i a3, __m256i b3, __m256i a4, __m256i b4, __m256i a5, __m256i b5, __m256i a6, __m256i b6)
{
	__m256i u00, u01, u10, u11, u20, u21, u30, u31;
	__m256i u40, u41, u50, u51, u60, u61, lo, hi;

	MUL_16to32(u00, u01, a0, b0);
	MUL_16to32(u10, u11, a1, b1);
	MUL_16to32(u20, u21, a2, b2);
	MUL_16to32(u30, u31, a3, b3);
	MUL_16to32(u40, u41, a4, b4);
	MUL_16to32(u50, u51, a5, b5);
	MUL_16to32(u60, u61, a6, b6);
	REPACK(lo, hi,
		_mm256_add_epi32(
			_mm256_add_epi32(
				_mm256_add_epi32(u00, u10),
				_mm256_add_epi32(u20, u30)),
			_mm256_add_epi32(
				_mm256_add_epi32(u40, u50),
				u60)),
		_mm256_add_epi32(
			_mm256_add_epi32(
				_mm256_add_epi32(u01, u11),
				_mm256_add_epi32(u21, u31)),
			_mm256_add_epi32(
				_mm256_add_epi32(u41, u51),
				u61)));
	return mq_montyred_x16(lo, hi);
}

static inline __m256i mq_montyLC8_x16(__m256i a0, __m256i b0, __m256i a1, __m256i b1, __m256i a2, __m256i b2, __m256i a3, __m256i b3, __m256i a4, __m256i b4, __m256i a5, __m256i b5, __m256i a6, __m256i b6, __m256i a7, __m256i b7)
{
	__m256i u00, u01, u10, u11, u20, u21, u30, u31;
	__m256i u40, u41, u50, u51, u60, u61, u70, u71, lo, hi;

	MUL_16to32(u00, u01, a0, b0);
	MUL_16to32(u10, u11, a1, b1);
	MUL_16to32(u20, u21, a2, b2);
	MUL_16to32(u30, u31, a3, b3);
	MUL_16to32(u40, u41, a4, b4);
	MUL_16to32(u50, u51, a5, b5);
	MUL_16to32(u60, u61, a6, b6);
	MUL_16to32(u70, u71, a7, b7);
	REPACK(lo, hi,
		_mm256_add_epi32(
			_mm256_add_epi32(
				_mm256_add_epi32(u00, u10),
				_mm256_add_epi32(u20, u30)),
			_mm256_add_epi32(
				_mm256_add_epi32(u40, u50),
				_mm256_add_epi32(u60, u70))),
		_mm256_add_epi32(
			_mm256_add_epi32(
				_mm256_add_epi32(u01, u11),
				_mm256_add_epi32(u21, u31)),
			_mm256_add_epi32(
				_mm256_add_epi32(u41, u51),
				_mm256_add_epi32(u61, u71))));
	return mq_montyred_x16(lo, hi);
}

static void mq_poly_mul_ntt(uint16_t *d, const uint16_t *a, const uint16_t *b, unsigned logn)
{

	size_t u;

	if (logn <= 3) {
		size_t n;

		n = (size_t)1 << logn;
		for (u = 0; u < n; u ++) {
			d[u] = mq_montymul(a[u], b[u]);
		}
		return;
	} else if (logn <= 7) {
		size_t n;

		n = (size_t)1 << logn;
		for (u = 0; u < n; u += 16) {
			_mm256_storeu_si256((void *)(d + u),
				mq_montymul_x16(
					_mm256_loadu_si256((void *)(a + u)),
					_mm256_loadu_si256((void *)(b + u))));
		}
		return;
	}

	for (u = 0; u < 128; u += 16) {
		__m256i a0, a1, a2, a3, a4, a5, a6, a7;
		__m256i b0, b1, b2, b3, b4, b5, b6, b7;
		__m256i x;

		a0 = _mm256_loadu_si256((void *)(a + u));
		a1 = _mm256_loadu_si256((void *)(a + u + 128));
		a2 = _mm256_loadu_si256((void *)(a + u + 256));
		a3 = _mm256_loadu_si256((void *)(a + u + 384));
		a4 = _mm256_loadu_si256((void *)(a + u + 512));
		a5 = _mm256_loadu_si256((void *)(a + u + 640));
		a6 = _mm256_loadu_si256((void *)(a + u + 768));
		a7 = _mm256_loadu_si256((void *)(a + u + 896));
		b0 = _mm256_loadu_si256((void *)(b + u));
		b1 = _mm256_loadu_si256((void *)(b + u + 128));
		b2 = _mm256_loadu_si256((void *)(b + u + 256));
		b3 = _mm256_loadu_si256((void *)(b + u + 384));
		b4 = _mm256_loadu_si256((void *)(b + u + 512));
		b5 = _mm256_loadu_si256((void *)(b + u + 640));
		b6 = _mm256_loadu_si256((void *)(b + u + 768));
		b7 = _mm256_loadu_si256((void *)(b + u + 896));
		x = _mm256_loadu_si256((void *)(NX + u));
		_mm256_storeu_si256((void *)(d + u),
			mq_montyLC2_x16(
				a0, b0,
				x, mq_montyLC7_x16(
					a1, b7, a2, b6, a3, b5,
					a4, b4, a5, b3, a6, b2,
					a7, b1)));
		_mm256_storeu_si256((void *)(d + u + 128),
			mq_montyLC3_x16(
				a0, b1, a1, b0,
				x, mq_montyLC6_x16(
					a2, b7, a3, b6, a4, b5,
					a5, b4, a6, b3, a7, b2)));
		_mm256_storeu_si256((void *)(d + u + 256),
			mq_montyLC4_x16(
				a0, b2, a1, b1, a2, b0,
				x, mq_montyLC5_x16(
					a3, b7, a4, b6, a5, b5,
					a6, b4, a7, b3)));
		_mm256_storeu_si256((void *)(d + u + 384),
			mq_montyLC5_x16(
				a0, b3, a1, b2, a2, b1, a3, b0,
				x, mq_montyLC4_x16(
					a4, b7, a5, b6, a6, b5,
					a7, b4)));
		_mm256_storeu_si256((void *)(d + u + 512),
			mq_montyLC6_x16(
				a0, b4, a1, b3, a2, b2, a3, b1, a4, b0,
				x, mq_montyLC3_x16(
					a5, b7, a6, b6, a7, b5)));
		_mm256_storeu_si256((void *)(d + u + 640),
			mq_montyLC7_x16(
				a0, b5, a1, b4, a2, b3, a3, b2,
				a4, b1, a5, b0,
				x, mq_montyLC2_x16(a6, b7, a7, b6)));
		_mm256_storeu_si256((void *)(d + u + 768),
			mq_montyLC8_x16(
				a0, b6, a1, b5, a2, b4, a3, b3,
				a4, b2, a5, b1, a6, b0,
				x, mq_montymul_x16(a7, b7)));
		_mm256_storeu_si256((void *)(d + u + 896),
			mq_montyLC8_x16(
				a0, b7, a1, b6, a2, b5, a3, b4,
				a4, b3, a5, b2, a6, b1, a7, b0));
	}
}

static void iNTT(uint16_t *d, const uint16_t *a, unsigned logn)
{
	if (logn < 7) {
		unsigned n, t, m;
		uint32_t ni;

		n = 1u << logn;
		if (d != a) {
			memmove(d, a, n * sizeof *a);
		}
		t = 1;
		m = n;
		while (m > 1) {
			unsigned hm, dt, i, j1;

			hm = m >> 1;
			dt = t << 1;
			for (i = 0, j1 = 0; i < hm; i ++, j1 += dt) {
				unsigned j, j2;
				uint32_t s;

				j2 = j1 + t;
				s = iGM[hm + i];
				for (j = j1; j < j2; j ++) {
					uint32_t u, v;

					u = d[j];
					v = d[j + t];
					d[j] = mq_add(u, v);
					d[j + t] = mq_montyred((Q + u - v) * s);
				}
			}
			t = dt;
			m = hm;
		}

		/*
		 * We need to divide by n, which we do by a multiplication by
		 * 1/n. Since we use Montgomery representation with R = 2^32,
		 * we use 2^32/n = 2^(32 - logn). We start with 2^54 mod q,
		 * then do a left shift by 10-logn (thus yielding a value equal
		 * to 2^(64-logn) modulo q), and apply a Montgomery reduction.
		 */
		ni = mq_montyred(306 << (10 - logn));
		for (m = 0; m < n; m ++) {
			d[m] = mq_montymul(d[m], ni);
		}
		return;
	}

	iNTT128(d, a);
	iNTT128(d + 128, a + 128);
	iNTT128(d + 256, a + 256);
	iNTT128(d + 384, a + 384);
	iNTT128(d + 512, a + 512);
	iNTT128(d + 640, a + 640);
	iNTT128(d + 768, a + 768);
	iNTT128(d + 896, a + 896);
	repack1024(d, d);
	return;
}

int end_make_public_769(uint16_t *h, const int8_t *f, const int8_t *g, unsigned logn, uint32_t *tmp)
{
	size_t u, n;
	uint16_t *ft, *gt;

	n = (size_t)1 << logn;
	ft = (uint16_t *)tmp;
	gt = ft + n;
	for (u = 0; u < n; u ++) {
		ft[u] = mq_set(f[u]);
		gt[u] = mq_set(g[u]);
	}
	NTT(ft, ft, logn);
	if (!mq_poly_inv_ntt(ft, ft, logn)) {
		return 0;
	}
	NTT(gt, gt, logn);
	mq_poly_mul_ntt(ft, ft, gt, logn);
	iNTT(ft, ft, logn);
	for (u = 0; u < n; u ++) {
		h[u] = mq_unorm(ft[u]);
	}
	return 1;
}