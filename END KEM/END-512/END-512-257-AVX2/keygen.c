/*******************************************************************************
 * 
 * Original source:
 *    ntrugen by Thomas Pornin,
 *    https://github.com/pornin/ntrugen
 * 
 * Modifications:
 * 	  Modified the values of the "SOLVE_END_257_512" array and the polynomial 
 *    sampling function to align with our parameter set.
 *
 *******************************************************************************/

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "param.h"
#include "module_p.h"
#include "big_integer.h"
#include "ntru.h"
#include "mq_ntt.h"
#include "poly.h"
#include "pke.h"
#include "poly.h"
#include "mq_ntt_param.h"
#include "randombytes.h"
#include "code.h"
#include "fips202.h"

const ntru_profile SOLVE_END_257_512 = {
	257,
	9, 9,
	{ 1, 1, 1, 2, 2, 4,  8,  15, 28, 56, 0},
	{ 1, 1, 2, 3, 6, 11, 21, 42, 82, 0   },
	{ 1, 1, 2, 2, 3, 3, 3, 4, 4, 0 },
	4,
	{ 0, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};



/*
 * Compute the vector w:
 *   w = round(q'*(gamma^2*F*adj(f) + G*adj(g))/(gamma^2*f*adj(f) + g*adj(g)))
 * with:
 *   q' = 12289
 *   gamma = 6 (for logn = 8 or 9) or sqrt(5) (for logn = 10)
 * Returned value: 1 on success, or 0 on error.
 *
 * Value of w is written at the start of tmp[] (in plain 32-bit format).
 *
 * Errors include an invalid parameter (logn must be 8, 9 or 10), or a
 * computation overflow. In the latter case, the function may return early;
 * thus, for constant-time discipline, any error should induce rejection
 * of the (f,g,F,G) quadruplet.
 *
 * RAM USAGE: 4.5*n words
 */
static int
compute_w(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	uint32_t *restrict tmp)
{
	uint32_t qp = 12289;

	int32_t gamma2;
	gamma2 = 6;


	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/*
	 * We create four buffers of size n words and one extra buffer of
	 * size n/2 words.
	 */
	uint32_t *t1 = tmp;
	uint32_t *t2 = t1 + n;
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	uint32_t *t5 = t4 + n;   /* half-size */

	/*
	 * We want:
	 *   t1 <- (gamma^2)*F*adj(f) + G*adj(g)  (RNS+NTT)
	 *   t2 <- (gamma^2)*f*adj(f) + g*adj(g)  (RNS+NTT)
	 * We work modulo a big 31-bit prime, which is large enough to
	 * avoid wrap-arounds, and thus allows us to get the plain integer
	 * values (after normalization).
	 */
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;
	uint32_t *gm = t1;
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);

	/*
	 * gmv <- R*gamma^2  (Montgomery representation of gamma^2)
	 */
	uint32_t gmv = mp_montymul(R2, mp_set(gamma2, p), p, p0i);

	/*
	 * t2 <- f  (RNS+NTT)
	 * t3 <- g  (RNS+NTT)
	 * t4 <- F  (RNS+NTT)
	 */
	poly_mp_set_small(logn, t2, f, p);
	poly_mp_set_small(logn, t3, g, p);
	poly_mp_set_small(logn, t4, F, p);
	mp_NTT(logn, t2, gm, p, p0i);
	mp_NTT(logn, t3, gm, p, p0i);
	mp_NTT(logn, t4, gm, p, p0i);

	/*
	 * t5 <- (gamma^2)*f*adj(f) + g*adj(g)  (RNS+NTT)  (auto-adjoint)
	 * t4 <- (gamma^2)*F*adj(f)             (RNS+NTT)
	 */
	for (size_t u = 0; u < hn; u ++) {
		uint32_t xf = t2[u];
		uint32_t xfa = t2[n - 1 - u];
		uint32_t xg = t3[u];
		uint32_t xga = t3[n - 1 - u];
		uint32_t xF = t4[u];
		uint32_t xFa = t4[n - 1 - u];
		uint32_t gmvf = mp_montymul(gmv, xf, p, p0i);
		uint32_t gmvfa = mp_montymul(gmv, xfa, p, p0i);
		t5[u] = mp_montymul(R2, mp_add(
			mp_montymul(xf, gmvfa, p, p0i),
			mp_montymul(xg, xga, p, p0i), p), p, p0i);
		t4[u] = mp_montymul(R2,
			mp_montymul(xF, gmvfa, p, p0i), p, p0i);
		t4[n - 1 - u] = mp_montymul(R2,
			mp_montymul(xFa, gmvf, p, p0i), p, p0i);
	}

	/*
	 * t1 <- (gamma^2)*F*adj(f) + G*adj(g)  (RNS+NTT)
	 */
	poly_mp_set_small(logn, t2, G, p);
	mp_NTT(logn, t2, gm, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		uint32_t xG = t2[u];
		uint32_t xga = t3[n - 1 - u];
		uint32_t xz = t4[u];
		uint32_t mga = mp_montymul(xga, R2, p, p0i);
		t1[u] = mp_add(mp_montymul(mga, xG, p, p0i), xz, p);
	}

	/*
	 * t2 <- (gamma^2)*f*adj(f) + g*adj(j)  (RNS+NTT) (full-size)
	 */
	for (size_t u = 0; u < hn; u ++) {
		t2[u] = t2[n - 1 - u] = t5[u];
	}

	/*
	 * Convert t1 and t2 to plain (32-bit).
	 */
	uint32_t *igm = t3;
	mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t1, igm, p, p0i);
	mp_iNTT(logn, t2, igm, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = (uint32_t)mp_norm(t1[u], p);
		t2[u] = (uint32_t)mp_norm(t2[u], p);
	}

	/*
	 * For the division, we go into the FFT domain. We check that the
	 * FFT won't overflow. The dividend is scaled down by 10 bits to
	 * compensate for the multiplication by q'.
	 *
	 * Buffer reorganization:
	 *    t1    (gamma^2)*F*adj(f) + G*adj(g) (plain, 32-bit) (n)
	 *    t2    (gamma^2)*f*adj(f) + g*adj(g) (plain, 32-bit) (n)
	 *    rt1   receives the dividend (FFT) (n fxr = 2*n)
	 *
	 * rt2 is an alias on t2 (as hn fxr values) and recieves the
	 * divisor (in FFT). The divisor is auto-adjoint, so it uses
	 * only half of the space in FFT representation; we compute in
	 * rt1 then move it to rt2.
	 */
	fxr *rt1 = (fxr *)(t2 + n);
	for (size_t u = 0; u < n; u ++) {
		rt1[u] = fxr_of(*(int32_t *)&t2[u]);
	}
	vect_FFT(logn, rt1);
	fxr *rt2 = (fxr *)t2;
	memmove(rt2, rt1, hn * sizeof *rt1);

	/*
	 * For the dividend, we multiply by q' but also scale down
	 * by 2^10; we check that the operation won't overflow.
	 */
	int32_t lim1 = (int32_t)(((uint64_t)1 << (41 - logn)) / qp);
	for (size_t u = 0; u < n; u ++) {
		int32_t x = *(int32_t *)&t1[u];
		if (x <= -lim1 || x >= +lim1) {
			return 0;
		}
		rt1[u] = fxr_of_scaled32(((uint64_t)x * qp) << 20);
	}
	vect_FFT(logn, rt1);

	/*
	 * Divisor is auto-adjoint. We inline the division loop here because
	 * we also want to check on overflows.
	 */
	for (size_t u = 0; u < hn; u ++) {
		fxr z1r = rt1[u];
		fxr z1i = rt1[u + hn];
		fxr z2 = rt2[u];
		if (!fxr_lt(fxr_div2e(fxr_abs(z1r), 30 - logn), z2)
			|| !fxr_lt(fxr_div2e(fxr_abs(z1i), 30 - logn), z2))
		{
			return 0;
		}
		rt1[u] = fxr_div(z1r, z2);
		rt1[u + hn] = fxr_div(z1i, z2);
	}
	vect_iFFT(logn, rt1);

	/*
	 * The unrounded w is in rt1 (scaled down by 2^10); we just have to
	 * round the coefficients and check that they are all in the
	 * allowed [-2^16..+2^16] range.
	 */
	fxr lim2 = fxr_of(1 << 3);
	for (size_t u = 0; u < n; u ++) {
		if (fxr_lt(lim2, fxr_abs(rt1[u]))) {
			return 0;
		}
		t1[u] = (uint32_t)fxr_round(fxr_mul2e(rt1[u], 12));
	}

	/*
	 * Check that the norm of (gamma*Fd, Gd) is low enough, with:
	 *   Fd = q'*F - f*w
	 *   Gd = q'*G - G*w
	 *
	 * Buffer layout:
	 *   t1    w (plain, 32-bit) (n)
	 *   t2    free (n)
	 *   t3    free (n)
	 *   t4    free (n)
	 * gm and igm are set to both point to t4 (for space-saving reasons,
	 * we use the same buffer for both values).
	 */
	gm = igm = t4;

	/*
	 * Convert w to NTT + Montgomery.
	 */
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mp_montymul(R2, mp_set(*(int32_t *)&t1[u], p), p, p0i);
	}
	mp_NTT(logn, t1, gm, p, p0i);

	/*
	 * t2 <- f     (NTT)
	 * t3 <- q'*F  (NTT)
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_set(f[u], p);
		t3[u] = mp_set((int32_t)qp * F[u], p);
	}
	mp_NTT(logn, t2, gm, p, p0i);
	mp_NTT(logn, t3, gm, p, p0i);

	/*
	 * t2 <- Fd, and compute its squared norm.
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_sub(t3[u], mp_montymul(t2[u], t1[u], p, p0i), p);
	}
	mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t2, igm, p, p0i);
	uint64_t Fdnorm = 0;
	for (size_t u = 0; u < n; u ++) {
		int32_t x = mp_norm(t2[u], p);
		Fdnorm += (uint64_t)((int64_t)x * (int64_t)x);
	}

	/*
	 * t2 <- g     (NTT)
	 * t3 <- q'*G  (NTT)
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_set(g[u], p);
		t3[u] = mp_set((int32_t)qp * G[u], p);
	}
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	mp_NTT(logn, t2, gm, p, p0i);
	mp_NTT(logn, t3, gm, p, p0i);

	/*
	 * t2 <- Gd, and compute its squared norm.
	 */
	for (size_t u = 0; u < n; u ++) {
		t2[u] = mp_sub(t3[u], mp_montymul(t2[u], t1[u], p, p0i), p);
	}
	mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t2, igm, p, p0i);
	uint64_t Gdnorm = 0;
	for (size_t u = 0; u < n; u ++) {
		int32_t x = mp_norm(t2[u], p);
		Gdnorm += (uint64_t)((int64_t)x * (int64_t)x);
	}

	/*
	 * Convert back w to plain 32-bit format.
	 */
	mp_iNTT(logn, t1, igm, p, p0i);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = (uint32_t)mp_norm(t1[u], p);
	}

	return 1;
}


void sample_binary(int8_t *x, keccak_state *state)
{
	shake256_squeeze(x, MESSBYTES, state);
}



void END_keygen(int8_t *h, int8_t *f, int8_t *g, int8_t *F, int8_t *G, int16_t *w, int8_t *k, void *tmp, size_t tmp_len)
{
	uint8_t __attribute__((aligned(32))) seed[SEEDLEN];
	randombytes(seed, SEEDLEN);
	keccak_state state;
    shake256_absorb_once(&state, seed, SEEDLEN);

	int16_t norm;
	int16_t __attribute__((aligned(32))) f_inv[DIM_N], tmp1[DIM_N], tmp2[DIM_N];
	sample_binary(k, &state);

	const ntru_profile *prof;
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint32_t *tt32 = (void *)utmp2;

	int8_t *tG = (int8_t *)tt32;
	int8_t *tF = tG + DIM_N;
	uint32_t *tw = (uint32_t *)(tF + DIM_N);

	for (;;) 
	{
		do 
		{
			ternary_sample_poly(f, &state);
			poly_int8_to_int16(f, tmp1);
            mq_poly_ntt_257_769(tmp1, nttdata_257);
		} while ((mq_poly_inv_ntt_257(f_inv, tmp1) != 0));

		do 
		{
			ternary_sample_poly(g, &state);
			poly_int8_to_int16(g, tmp1);
			mq_poly_ntt_257_769(tmp1, nttdata_257);
		} while ((mq_poly_inv_ntt_257(tmp2, tmp1) != 0));

		mq_poly_mul_ntt_257_769_mq(tmp2, f_inv, tmp1, muldata_257);
		encode_pk(tmp2, h);

		prof = &SOLVE_END_257_512;

		/*
		 * Solve the NTRU equation.
		 */
		int err = solve_NTRU(prof, LOGN, g, f, tt32);
		switch (err) 
		{
			case SOLVE_OK:
				break;
			default:
				continue;
		}

		if (!compute_w(LOGN, f, g, tF, tG, tw)) 
		{
			continue;
		}

		memcpy(F, tF, DIM_N);
		memcpy(G, tG, DIM_N);

		for(int i = 0;i < DIM_N; i++)
		{
			w[i] = ((int32_t)tw[i]) / (20478);
		}

		break;
	}
}
