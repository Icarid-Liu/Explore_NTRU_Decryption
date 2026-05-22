#ifndef MQ_NTT_H
#define MQ_NTT_H

#include <stdint.h>
#include <immintrin.h>
#include "param.h"

void mq_poly_ntt_257_769(int16_t *a, int16_t *nttdata);
void mq_poly_intt_257_769(int16_t *a, int16_t *inttdata);
void mq_poly_intt_257_769_n(int16_t *a, int16_t *inttdata);
void mq_poly_mul_ntt_257_769(int16_t *r, int16_t *a, int16_t *b, int16_t *nttdata);
void mq_poly_mul_ntt_257_769_mq(int16_t *r, int16_t *a, int16_t *b, int16_t *nttdata);
int mq_poly_inv_ntt_257(int16_t *f_inv,int16_t *f);
void mq_poly_adjoint_ntt_769(int16_t *r, int16_t *a);
void mq_poly_ntt_12289(int16_t *a);
void mq_poly_intt_12289(int16_t *a);
void mq_poly_intt_12289_n(int16_t *a);
void mq_poly_mul_ntt_12289(int16_t *r, int16_t *a, int16_t *b);
void mq_poly_ntt_25601(int16_t *a);
void mq_poly_intt_25601(int16_t *a);
void mq_poly_mul_ntt_25601(int16_t *r, int16_t *a, int16_t *b);

#endif