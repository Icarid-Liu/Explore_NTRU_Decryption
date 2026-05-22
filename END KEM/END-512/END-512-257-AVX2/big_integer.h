#ifndef BIG_INTEGER_H
#define BIG_INTEGER_H

#include <stdio.h>
#include <stdint.h>

// #define __DEFAULT_FN_ATTRS __attribute__((__always_inline__, __nodebug__, __target__("lzcnt")))

// static __inline__ unsigned int __DEFAULT_FN_ATTRS
// __lzcnt32(unsigned int __X)
// {
//   return __builtin_ia32_lzcnt_u32(__X);
// }

// static __inline__ unsigned int __DEFAULT_FN_ATTRS
// _lzcnt_u32(unsigned int __X)
// {
//   return __builtin_ia32_lzcnt_u32(__X);
// }

/*
 * Get the number of leading zeros in a 32-bit value.
 */
static unsigned
lzcnt(uint32_t x)
{
	/*
	 * All AVX2-capable CPUs have lzcnt.
	 */
	return _lzcnt_u32(x);
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
__m256i zint_mod_small_signed_x8(const uint32_t *d, size_t len, size_t stride, __m256i yp, __m256i yp0i, __m256i yR2, __m256i yRx);
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
