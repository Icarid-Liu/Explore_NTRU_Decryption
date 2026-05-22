#ifndef POLY_H
#define POLY_H

#include <stdio.h>
#include <stdint.h>
#include "fixed_point.h"


void poly_int8_to_int16(int8_t *a, int16_t *b);
void poly_int16_to_int8(int16_t *a, int8_t *b);
void poly_round(int16_t *a);
void poly_deround(int16_t *a);
int16_t poly_norm2_int8(int8_t *a);
int16_t poly_norm2(int16_t *a);
void poly_compute_adjoint(int8_t *fa, int8_t *f);
void poly_compute_adjoint_int16(int16_t *fa, int16_t *f);
int poly_check_distribution_binary(int16_t * s);
void poly_mp_set_small(unsigned logn, uint32_t *__restrict d, const int8_t *__restrict f, uint32_t p);
void poly_mp_set(unsigned logn, uint32_t *f, uint32_t p);
void poly_mp_norm(unsigned logn, uint32_t *f, uint32_t p);
int poly_big_to_small(unsigned logn, int8_t *__restrict d, const uint32_t *__restrict s, int lim);
uint32_t poly_max_bitlength(unsigned logn, const uint32_t *f, size_t flen);
void poly_big_to_fixed(unsigned logn, fxr *__restrict d, const uint32_t *__restrict f, size_t len, uint32_t sc);
void poly_sub_scaled(unsigned logn, uint32_t *__restrict F, size_t Flen, const uint32_t *__restrict f, size_t flen, const int32_t *__restrict k, uint32_t sc);
void poly_sub_scaled_ntt(unsigned logn, uint32_t *__restrict F, size_t Flen, const uint32_t *__restrict f, size_t flen, const int32_t *__restrict k, uint32_t sc, uint32_t *__restrict tmp);
void poly_sub_kf_scaled_depth1(unsigned logn_top, uint32_t *restrict F, size_t FGlen, uint32_t *restrict k, uint32_t sc, const int8_t *restrict f, uint32_t *restrict tmp);
int poly_is_invertible(unsigned logn, const int8_t *__restrict f, uint32_t p, uint32_t p0i, uint32_t s, uint32_t r, uint32_t rm, unsigned rs, uint32_t *__restrict tmp);
int poly_is_invertible_ext(unsigned logn, const int8_t *__restrict f, uint32_t r1, uint32_t r2, uint32_t p, uint32_t p0i, uint32_t s, uint32_t r1m, unsigned r1s, uint32_t r2m, unsigned r2s, uint32_t *__restrict tmp);
uint32_t poly_sqnorm(unsigned logn, const int8_t *f);

#endif


