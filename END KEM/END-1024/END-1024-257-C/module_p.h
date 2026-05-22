/*******************************************************************************
 * 
 *    Original source:
 *    ntrugen by Thomas Pornin,
 *    https://github.com/pornin/ntrugen
 *
 *******************************************************************************/

#ifndef MODULE_P_H
#define MODULE_P_H

#include <stdint.h>

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

/*
 * Get the signed normalized value of x mod p.
 */
static inline int32_t mp_norm(uint32_t x, uint32_t p)
{
	uint32_t w = x - (p & tbmask((p >> 1) - x));
	return *(int32_t *)&w;
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

/*
 * Subtraction modulo p.
 */
static inline uint32_t mp_sub(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d = a - b;
	return d + (p & tbmask(d));
}

/*
 * Halving modulo p.
 */
static inline uint32_t mp_half(uint32_t a, uint32_t p)
{
	return (a + (p & -(a & 1))) >> 1;
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
void mp_mkgmigm(unsigned logn, uint32_t *__restrict gm, uint32_t *__restrict igm, uint32_t g, uint32_t ig, uint32_t p, uint32_t p0i);
void mp_mkgm(unsigned logn, uint32_t *__restrict gm, uint32_t g, uint32_t p, uint32_t p0i);
void mp_mkgm7(uint32_t *__restrict gm, uint32_t g, uint32_t p, uint32_t p0i);
void mp_mkigm(unsigned logn, uint32_t *__restrict igm, uint32_t ig, uint32_t p, uint32_t p0i);
void mp_NTT(unsigned logn, uint32_t *__restrict a, const uint32_t *__restrict gm, uint32_t p, uint32_t p0i);
void mp_iNTT(unsigned logn, uint32_t *__restrict a, const uint32_t *__restrict igm, uint32_t p, uint32_t p0i);
extern const small_prime PRIMES[];

#endif