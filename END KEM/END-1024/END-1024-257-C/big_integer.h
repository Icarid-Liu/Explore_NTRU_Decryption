/*******************************************************************************
 * 
 *    Original source:
 *    ntrugen by Thomas Pornin,
 *    https://github.com/pornin/ntrugen
 *
 *******************************************************************************/

#ifndef BIG_INTEGER_H
#define BIG_INTEGER_H

#include <stdio.h>
#include <stdint.h>
/*
 * Get the number of leading zeros in a 32-bit value.
 */
static inline unsigned lzcnt(uint32_t x)
{
	uint32_t m = tbmask((x >> 16) - 1);
	uint32_t s = m & 16;
	x = (x >> 16) ^ (m & (x ^ (x >> 16)));
	m = tbmask((x >>  8) - 1);
	s |= m &  8;
	x = (x >>  8) ^ (m & (x ^ (x >>  8)));
	m = tbmask((x >>  4) - 1);
	s |= m &  4;
	x = (x >>  4) ^ (m & (x ^ (x >>  4)));
	m = tbmask((x >>  2) - 1);
	s |= m &  2;
	x = (x >>  2) ^ (m & (x ^ (x >>  2)));

	/*
	 * At this point, x fits on 2 bits. Number of leading zeros is
	 * then:
	 *    x = 0   -> 2
	 *    x = 1   -> 1
	 *    x = 2   -> 0
	 *    x = 3   -> 0
	 */
	return (unsigned)(s + ((2 - x) & tbmask(x - 3)));
}

/*
 * Identical to lzcnt(), except that the caller makes sure that the
 * operand is non-zero. On (old-ish) x86 systems, this function could be
 * specialized with the bsr opcode (which does not support a zero input).
 */
#define lzcnt_nonzero   lzcnt

#define DIVREM31(q, r, x)  { \
		uint32_t divrem31_q, divrem31_x; \
		divrem31_x = (x); \
		divrem31_q = (uint32_t)(divrem31_x * (uint32_t)67651) >> 21; \
		(q) = divrem31_q; \
		(r) = divrem31_x - 31 * divrem31_q; \
	} while (0)

uint32_t zint_mul_small(uint32_t *m, size_t len, uint32_t x);
uint32_t zint_mod_small_unsigned(const uint32_t *d, size_t len, size_t stride, uint32_t p, uint32_t p0i, uint32_t R2);
uint32_t zint_mod_small_signed(const uint32_t *d, size_t len, size_t stride, uint32_t p, uint32_t p0i, uint32_t R2, uint32_t Rx);
void zint_add_mul_small(uint32_t *__restrict x, size_t len, size_t xstride, const uint32_t *__restrict y, uint32_t s);
void zint_norm_zero(uint32_t *__restrict x, size_t len, size_t xstride, const uint32_t *__restrict p);
void zint_rebuild_CRT(uint32_t *__restrict xx, size_t xlen, size_t n, size_t num_sets, int normalize_signed, uint32_t *__restrict tmp);
void zint_negate(uint32_t *a, size_t len, uint32_t ctl);
static uint32_t zint_co_reduce(uint32_t *__restrict a, uint32_t *__restrict b, size_t len, int64_t xa, int64_t xb, int64_t ya, int64_t yb);
static void zint_finish_mod(uint32_t *__restrict a, const uint32_t *__restrict m, size_t len, uint32_t neg);
static void zint_co_reduce_mod(uint32_t *__restrict a, uint32_t *__restrict b, const uint32_t *__restrict m, size_t len, uint32_t m0i, int64_t xa, int64_t xb, int64_t ya, int64_t yb);
int zint_bezout(uint32_t *__restrict u, uint32_t *__restrict v, const uint32_t *__restrict x, const uint32_t *__restrict y, size_t len, uint32_t *__restrict tmp);
void zint_add_scaled_mul_small(uint32_t *__restrict x, size_t xlen, const uint32_t *__restrict y, size_t ylen, size_t stride, int32_t k, uint32_t sch, uint32_t scl);
void zint_sub_scaled(uint32_t *__restrict x, size_t xlen, const uint32_t *__restrict y, size_t ylen, size_t stride, uint32_t sch, uint32_t scl);

#endif
