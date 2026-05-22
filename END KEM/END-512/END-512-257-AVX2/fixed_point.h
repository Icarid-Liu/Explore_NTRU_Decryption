#ifndef FIXED_POINT_H
#define FIXED_POINT_H

#include <stdint.h>

typedef struct {
	uint64_t v;
} fxr;
#define FXR(x)   { (x) }

static inline fxr fxr_of(int32_t j)
{
	fxr x;

	x.v = (uint64_t)j << 32;
	return x;
}

static inline fxr fxr_of_scaled32(uint64_t t)
{
	fxr x;

	x.v = t;
	return x;
}

static inline fxr fxr_add(fxr x, fxr y)
{
	x.v += y.v;
	return x;
}

static inline fxr fxr_sub(fxr x, fxr y)
{
	x.v -= y.v;
	return x;
}

static inline fxr fxr_double(fxr x)
{
	x.v <<= 1;
	return x;
}

static inline fxr fxr_neg(fxr x)
{
	x.v = -x.v;
	return x;
}

static inline fxr fxr_abs(fxr x)
{
	x.v -= (x.v << 1) & (uint64_t)(*(int64_t *)&x.v >> 63);
	return x;
}

static inline fxr fxr_mul(fxr x, fxr y)
{
#if defined __GNUC__ && defined __SIZEOF_INT128__
	__int128 z;

	z = (__int128)*(int64_t *)&x.v * (__int128)*(int64_t *)&y.v;
	x.v = (uint64_t)(z >> 32);
	return x;
#else
	int32_t xh, yh;
	uint32_t xl, yl;
	uint64_t z0, z1, z2, z3;

	xl = (uint32_t)x.v;
	yl = (uint32_t)y.v;
	xh = (int32_t)(*(int64_t *)&x.v >> 32);
	yh = (int32_t)(*(int64_t *)&y.v >> 32);
	z0 = ((uint64_t)xl * (uint64_t)yl) >> 32;
	z1 = (uint64_t)((int64_t)xl * (int64_t)yh);
	z2 = (uint64_t)((int64_t)yl * (int64_t)xh);
	z3 = (uint64_t)((int64_t)xh * (int64_t)yh) << 32;
	x.v = z0 + z1 + z2 + z3;
	return x;
#endif
}

static inline __m256i
fxr_mul_x4(__m256i ya, __m256i yb)
{
	__m256i ya_hi = _mm256_srli_epi64(ya, 32);
	__m256i yb_hi = _mm256_srli_epi64(yb, 32);
	__m256i y1 = _mm256_mul_epu32(ya, yb);
	__m256i y2 = _mm256_mul_epu32(ya, yb_hi);
	__m256i y3 = _mm256_mul_epu32(ya_hi, yb);
	__m256i y4 = _mm256_mul_epu32(ya_hi, yb_hi);
	y1 = _mm256_srli_epi64(y1, 32);
	y4 = _mm256_slli_epi64(y4, 32);
	__m256i y5 = _mm256_add_epi64(
		_mm256_add_epi64(y1, y2),
		_mm256_add_epi64(y3, y4));
	__m256i yna = _mm256_srai_epi32(ya, 31);
	__m256i ynb = _mm256_srai_epi32(yb, 31);
	return _mm256_sub_epi64(y5,
		_mm256_add_epi64(
			_mm256_and_si256(_mm256_slli_epi64(yb, 32), yna),
			_mm256_and_si256(_mm256_slli_epi64(ya, 32), ynb)));
}

static inline fxr fxr_sqr(fxr x)
{
#if defined __GNUC__ && defined __SIZEOF_INT128__
	int64_t t;
	__int128 z;

	t = *(int64_t *)&x.v;
	z = (__int128)t * (__int128)t;
	x.v = (uint64_t)(z >> 32);
	return x;
#else
	int32_t xh;
	uint32_t xl;
	uint64_t z0, z1, z3;

	xl = (uint32_t)x.v;
	xh = (int32_t)(*(int64_t *)&x.v >> 32);
	z0 = ((uint64_t)xl * (uint64_t)xl) >> 32;
	z1 = (uint64_t)((int64_t)xl * (int64_t)xh);
	z3 = (uint64_t)((int64_t)xh * (int64_t)xh) << 32;
	x.v = z0 + (z1 << 1) + z3;
	return x;
#endif
}

static inline __m256i
fxr_sqr_x4(__m256i ya)
{
	__m256i ya_hi = _mm256_srli_epi64(ya, 32);
	__m256i y1 = _mm256_mul_epu32(ya, ya);
	__m256i y2 = _mm256_mul_epu32(ya, ya_hi);
	__m256i y3 = _mm256_mul_epu32(ya_hi, ya_hi);
	y1 = _mm256_srli_epi64(y1, 32);
	y2 = _mm256_add_epi64(y2, y2);
	y3 = _mm256_slli_epi64(y3, 32);
	__m256i y4 = _mm256_add_epi64(_mm256_add_epi64(y1, y2), y3);
	return _mm256_sub_epi64(y4,
		_mm256_and_si256(_mm256_slli_epi64(ya, 33),
		_mm256_srai_epi32(ya, 31)));
}


static inline int32_t fxr_round(fxr x)
{
	x.v += 0x80000000ul;
	return (int32_t)(*(int64_t *)&x.v >> 32);
}

static inline fxr fxr_div2e(fxr x, unsigned n)
{
	x.v += (((uint64_t)1 << n) >> 1);
	x.v = (uint64_t)(*(int64_t *)&x.v >> n);
	return x;
}

static inline __m256i
fxr_half_x4(__m256i ya)
{
	const __m256i y1 = _mm256_set1_epi64x(1);
	const __m256i yh = _mm256_set1_epi64x((uint64_t)1 << 63);
	ya = _mm256_add_epi64(ya, y1);
	return _mm256_or_si256(
		_mm256_srli_epi64(ya, 1),
		_mm256_and_si256(ya, yh));
}

static inline fxr fxr_mul2e(fxr x, unsigned n)
{
	x.v <<= n;
	return x;
}

uint64_t inner_fxr_div(uint64_t x, uint64_t y);

static inline fxr fxr_inv(fxr x)
{
	x.v = inner_fxr_div((uint64_t)1 << 32, x.v);
	return x;
}

static inline fxr fxr_div(fxr x, fxr y)
{
	x.v = inner_fxr_div(x.v, y.v);
	return x;
}

__m256i fxr_div_x4(__m256i yn, __m256i yd);

static inline void
fxr_div_x4_1(fxr *n0, fxr *n1, fxr *n2, fxr *n3, fxr d)
{
	__m256i yn = _mm256_setr_epi64x(n0->v, n1->v, n2->v, n3->v);
	__m256i yd = _mm256_set1_epi64x(d.v);
	union {
		__m256i y;
		uint64_t q[4];
	} z;
	z.y = fxr_div_x4(yn, yd);
	n0->v = z.q[0];
	n1->v = z.q[1];
	n2->v = z.q[2];
	n3->v = z.q[3];
}

static inline int fxr_lt(fxr x, fxr y)
{
	return *(int64_t *)&x.v < *(int64_t *)&y.v;
}

static const fxr fxr_zero = { 0 };
static const fxr fxr_sqrt2 = { 6074001000ull };

/*
 * A complex value.
 */
typedef struct {
	fxr re, im;
} fxc;
#define FXC(re, im)   { FXR(re), FXR(im) }

static inline fxc fxc_add(fxc x, fxc y)
{
	x.re = fxr_add(x.re, y.re);
	x.im = fxr_add(x.im, y.im);
	return x;
}

static inline fxc fxc_sub(fxc x, fxc y)
{
	x.re = fxr_sub(x.re, y.re);
	x.im = fxr_sub(x.im, y.im);
	return x;
}

static inline fxc fxc_half(fxc x)
{
	x.re = fxr_div2e(x.re, 1);
	x.im = fxr_div2e(x.im, 1);
	return x;
}

static inline fxc fxc_mul(fxc x, fxc y)
{
	/*
	 * We are computing r = (a + i*b)*(c + i*d) with:
	 *   z0 = a*c
	 *   z1 = b*d
	 *   z2 = (a + b)*(c + d)
	 *   r = (z0 - z1) + i*(z2 - (z0 + z1))
	 * Since the intermediate values are truncated to our precision,
	 * the imaginary value of r _may_ be slightly different from
	 * a*d + b*c (if we had calculated it directly). For full
	 * reproducibility, all implementations should use the formulas
	 * above.
	 */
	fxr z0 = fxr_mul(x.re, y.re);
	fxr z1 = fxr_mul(x.im, y.im);
	fxr z2 = fxr_mul(fxr_add(x.re, x.im), fxr_add(y.re, y.im));
	fxc z;
	z.re = fxr_sub(z0, z1);
	z.im = fxr_sub(z2, fxr_add(z0, z1));
	return z;
}

static inline void
fxc_mul_x4(__m256i *yd_re, __m256i *yd_im,
	__m256i ya_re, __m256i ya_im, __m256i yb_re, __m256i yb_im)
{
	__m256i y0 = fxr_mul_x4(ya_re, yb_re);
	__m256i y1 = fxr_mul_x4(ya_im, yb_im);
	__m256i y2 = fxr_mul_x4(
		_mm256_add_epi64(ya_re, ya_im),
		_mm256_add_epi64(yb_re, yb_im));
	*yd_re = _mm256_sub_epi64(y0, y1);
	*yd_im = _mm256_sub_epi64(y2, _mm256_add_epi64(y0, y1));
}

static inline fxc fxc_conj(fxc x)
{
	x.im = fxr_neg(x.im);
	return x;
}


uint64_t inner_fxr_div(uint64_t x, uint64_t y);
void vect_FFT(unsigned logn, fxr *f);
void vect_iFFT(unsigned logn, fxr *f);
void vect_set(unsigned logn, fxr *d, const int8_t *f);
void vect_add(unsigned logn, fxr *__restrict a, const fxr *__restrict b);
void vect_mul_realconst(unsigned logn, fxr *a, fxr c);
void vect_mul2e(unsigned logn, fxr *a, unsigned e);
void vect_mul_fft(unsigned logn, fxr *__restrict a, const fxr *__restrict b);
void vect_adj_fft(unsigned logn, fxr *a);
void vect_mul_autoadj_fft(unsigned logn, fxr *__restrict a, const fxr *__restrict b);
void vect_div_autoadj_fft(unsigned logn, fxr *__restrict a, const fxr *__restrict b);
void vect_invnorm_fft(unsigned logn, fxr *__restrict d, const fxr *__restrict a, const fxr *__restrict b, unsigned e);
void vect_inv_mul2e_fft(unsigned logn, fxr *a, unsigned e);
#endif

