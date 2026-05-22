/*******************************************************************************
 * 
 *    Original source:
 *    ntrugen by Thomas Pornin,
 *    https://github.com/pornin/ntrugen
 *
 *******************************************************************************/

#ifndef NTRU_H
#define NTRU_H

#include <stdio.h>
#include <stdint.h>

typedef struct {
	uint32_t q;
	unsigned min_logn, max_logn;
	uint16_t max_bl_small[11];
	uint16_t max_bl_large[10];
	uint16_t word_win[10];
	uint32_t reduce_bits;
	uint8_t coeff_FG_limit[11];
	uint16_t min_save_fg[11];
} ntru_profile;

/* Error code: no error (so far) */
#define SOLVE_OK           0

/* Error code: GCD(Res(f,X^n+1), Res(g,X^n+1)) != 1 */
#define SOLVE_ERR_GCD      -1

/* Error code: reduction error (NTRU equation no longer fulfilled) */
#define SOLVE_ERR_REDUCE   -2

/* Error code: output (F,G) coefficients are off-limits */
#define SOLVE_ERR_LIMIT    -3

/*
 * We use poly_sub_scaled() when log(n) < MIN_LOGN_FGNTT, and
 * poly_sub_scaled_ntt() when log(n) >= MIN_LOGN_FGNTT. The NTT variant
 * is faster at large degrees, but not at small degrees.
 */
#define MIN_LOGN_FGNTT   4

static void make_fg_zero(unsigned logn, const int8_t *__restrict f, const int8_t *__restrict g, uint32_t *__restrict tmp);
static void make_fg_step(const ntru_profile *prof, unsigned logn_top, unsigned depth, uint32_t *tmp);
static void make_fg_intermediate(const ntru_profile *prof, unsigned logn_top, const int8_t *__restrict f, const int8_t *__restrict g, unsigned depth, uint32_t *tmp);
static int make_fg_deepest(const ntru_profile *prof, unsigned logn_top, const int8_t *__restrict f, const int8_t *__restrict g, uint32_t *tmp, size_t sav_off);
static int solve_NTRU_deepest(const ntru_profile *prof, unsigned logn_top, const int8_t *__restrict f, const int8_t *__restrict g, uint32_t *tmp);
static int solve_NTRU_intermediate(const ntru_profile *__restrict prof, unsigned logn_top, const int8_t *__restrict f, const int8_t *__restrict g, unsigned depth, uint32_t *__restrict tmp);
static int solve_NTRU_depth0(const ntru_profile *__restrict prof, unsigned logn, const int8_t *__restrict f, const int8_t *__restrict g, uint32_t *__restrict tmp);
int solve_NTRU(const ntru_profile *__restrict prof, unsigned logn, const int8_t *__restrict f, const int8_t *__restrict g, uint32_t *tmp);
int recover_G(unsigned logn, int32_t q, uint32_t ulim, const int8_t *__restrict f, const int8_t *__restrict g, const int8_t *__restrict F, uint32_t *__restrict tmp);


#endif