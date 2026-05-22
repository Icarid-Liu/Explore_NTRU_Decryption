#ifndef MODULE_P_H
#define MODULE_P_H

#include <stdint.h>
#include <immintrin.h>

typedef struct {
	uint32_t p;
	uint32_t p0i;
	uint32_t R2;
	uint32_t g;
	uint32_t ig;
	uint32_t s;
} small_prime;

static inline uint32_t tbmask(uint32_t x)
{
	return (uint32_t)(*(int32_t *)&x >> 31);
}

/*
 * Get v mod p in the 0 to p-1 range; input v must be in the -(p-1) to +(p-1)
 * range.
 */
static inline uint32_t mp_set(int32_t v, uint32_t p)
{
	uint32_t w = (uint32_t)v;
	return w + (p & tbmask(w));
}

static inline __m256i
mp_set_x8(__m256i yv, __m256i yp)
{
	return _mm256_add_epi32(yv, _mm256_and_si256(yp,
		_mm256_srai_epi32(yv, 31)));
}

/*
 * Get the signed normalized value of x mod p.
 */
static inline int32_t mp_norm(uint32_t x, uint32_t p)
{
	uint32_t w = x - (p & tbmask((p >> 1) - x));
	return *(int32_t *)&w;
}

static inline __m256i
mp_norm_x8(__m256i yv, __m256i yp, __m256i yhp)
{
	return _mm256_sub_epi32(yv, _mm256_and_si256(yp,
		_mm256_cmpgt_epi32(yv, yhp)));
}

/*
 * Compute R = 2^32 mod p.
 */
static inline uint32_t mp_R(uint32_t p)
{
	/*
	 * Since 2*p < 2^32 < 3*p, we just subtract 2*p from 2^32.
	 */
	return -(p << 1);
}

/*
 * Compute R/2 = 2^31 mod p.
 */
static inline uint32_t mp_hR(uint32_t p)
{
	/*
	 * Since p < 2^31 < (3/2)*p, we just subtract p from 2^31.
	 */
	return ((uint32_t)1 << 31) - p;
}

/*
 * Addition modulo p.
 */
static inline uint32_t mp_add(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d = a + b - p;
	return d + (p & tbmask(d));
}

static inline __m256i
mp_add_x8(__m256i ya, __m256i yb, __m256i yp)
{
	__m256i yd = _mm256_sub_epi32(_mm256_add_epi32(ya, yb), yp);
	return _mm256_add_epi32(yd, _mm256_and_si256(yp,
		_mm256_srai_epi32(yd, 31)));
}

/*
 * Subtraction modulo p.
 */
static inline uint32_t mp_sub(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d = a - b;
	return d + (p & tbmask(d));
}

static inline __m256i
mp_sub_x8(__m256i ya, __m256i yb, __m256i yp)
{
	__m256i yd = _mm256_sub_epi32(ya, yb);
	return _mm256_add_epi32(yd, _mm256_and_si256(yp,
		_mm256_srai_epi32(yd, 31)));
}

/*
 * Halving modulo p.
 */
static inline uint32_t mp_half(uint32_t a, uint32_t p)
{
	return (a + (p & -(a & 1))) >> 1;
}

static inline __m256i
mp_half_x8(__m256i ya, __m256i yp)
{
	return _mm256_srli_epi32(
		_mm256_add_epi32(ya, _mm256_and_si256(yp,
			_mm256_sub_epi32(_mm256_setzero_si256(),
			_mm256_and_si256(ya, _mm256_set1_epi32(1))))), 1);
}

/*
 * Montgomery multiplication modulo p.
 *
 * Reduction computes (a*b + w*p)/(2^32) for some w <= 2^(32-1);
 * then p is conditionally subtracted. This process works as long as:
 *    (a*b + p*(2^32-1))/(2^32) <= 2*p-1
 * which holds if:
 *    a*b <= p*2^32 - 2^32 + p
 * This works if both a and b are proper integers modulo p (in the 0 to p-1
 * range), but also if, for instance, a is an integer modulo p, and b is an
 * arbitrary 32-bit integer.
 */
static inline uint32_t mp_montymul(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i)
{
	uint64_t z = (uint64_t)a * (uint64_t)b;
	uint32_t w = (uint32_t)z * p0i;
	uint32_t d = (uint32_t)((z + (uint64_t)w * (uint64_t)p) >> 32) - p;
	return d + (p & tbmask(d));
}

static inline __m256i
mp_montymul_x4(__m256i ya, __m256i yb, __m256i yp, __m256i yp0i)
{
	__m256i yd = _mm256_mul_epu32(ya, yb);
	__m256i ye = _mm256_mul_epu32(yd, yp0i);
	ye = _mm256_mul_epu32(ye, yp);
	yd = _mm256_srli_epi64(_mm256_add_epi64(yd, ye), 32);
	yd = _mm256_sub_epi32(yd, yp);
	return _mm256_add_epi32(yd, _mm256_and_si256(yp,
		_mm256_srai_epi32(yd, 31)));
}

static inline __m256i
mp_montymul_x8(__m256i ya, __m256i yb, __m256i yp, __m256i yp0i)
{
	/* yd0 <- a0*b0 : a2*b2 (+high lane) */
	__m256i yd0 = _mm256_mul_epu32(ya, yb);
	/* yd1 <- a1*b1 : a3*b3 (+high lane) */
	__m256i yd1 = _mm256_mul_epu32(
		_mm256_srli_epi64(ya, 32),
		_mm256_srli_epi64(yb, 32));

	__m256i ye0 = _mm256_mul_epu32(yd0, yp0i);
	__m256i ye1 = _mm256_mul_epu32(yd1, yp0i);
	ye0 = _mm256_mul_epu32(ye0, yp);
	ye1 = _mm256_mul_epu32(ye1, yp);
	yd0 = _mm256_add_epi64(yd0, ye0);
	yd1 = _mm256_add_epi64(yd1, ye1);

	/* yf0 <- lo(d0) : lo(d1) : hi(d0) : hi(d1) (+high lane) */
	__m256i yf0 = _mm256_unpacklo_epi32(yd0, yd1);
	/* yf1 <- lo(d2) : lo(d3) : hi(d2) : hi(d3) (+high lane) */
	__m256i yf1 = _mm256_unpackhi_epi32(yd0, yd1);
	/* yg <- hi(d0) : hi(d1) : hi(d2) : hi(d3) (+high lane) */
	__m256i yg = _mm256_unpackhi_epi64(yf0, yf1);
	/*
	 * Alternate version (instead of the three unpack above) but it
	 * seems to be slightly slower.
	__m256i yg = _mm256_blend_epi32(_mm256_srli_epi64(yd0, 32), yd1, 0xAA);
	 */

	yg = _mm256_sub_epi32(yg, yp);
	return _mm256_add_epi32(yg, _mm256_and_si256(yp,
		_mm256_srai_epi32(yg, 31)));
}

/*
 * Compute 2^(31*e) mod p.
 */
static inline uint32_t mp_Rx31(unsigned e, uint32_t p, uint32_t p0i, uint32_t R2)
{
	/* x <- 2^63 mod p = Montgomery representation of 2^31 */
	uint32_t x = mp_half(R2, p);
	uint32_t d = 1;
	for (;;) {
		if ((e & 1) != 0) {
			d = mp_montymul(d, x, p, p0i);
		}
		e >>= 1;
		if (e == 0) {
			return d;
		}
		x = mp_montymul(x, x, p, p0i);
	}
}

uint32_t mp_div(uint32_t x, uint32_t y, uint32_t p);
__m256i mp_div_x8(__m256i ynum, __m256i yden, __m256i yp);
void mp_mkgmigm(unsigned logn, uint32_t *__restrict gm, uint32_t *__restrict igm, uint32_t g, uint32_t ig, uint32_t p, uint32_t p0i);
void mp_mkgm(unsigned logn, uint32_t *__restrict gm, uint32_t g, uint32_t p, uint32_t p0i);
void mp_mkgm7(uint32_t *__restrict gm, uint32_t g, uint32_t p, uint32_t p0i);
void mp_mkigm(unsigned logn, uint32_t *__restrict igm, uint32_t ig, uint32_t p, uint32_t p0i);
void mp_NTT(unsigned logn, uint32_t *__restrict a, const uint32_t *__restrict gm, uint32_t p, uint32_t p0i);
void mp_iNTT(unsigned logn, uint32_t *__restrict a, const uint32_t *__restrict igm, uint32_t p, uint32_t p0i);
extern const small_prime PRIMES[];

#endif