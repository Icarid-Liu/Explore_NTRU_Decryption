#include <string.h>
#include <stdio.h>
#include "param.h"
#include "mq_ntt.h"
#include "pke.h"
#include "poly.h"
#include "code.h"
#include "randombytes.h"
#include "fips202.h"

static inline int verify(int16_t *a, int16_t *b)
{
  int i;
  int16_t r = 0;

  for(i = 0; i < DIM_N; i++)
    r |= a[i] ^ b[i];

  return (-(uint64_t)r) >> 63;
}

static inline void cmov(int8_t *r, int8_t *x, int len, int b)
{
  int i;

  b = -b;
  for(i = 0; i < len; i++)
    r[i] ^= b & (r[i] ^ x[i]);
}


void END_cca_encrypt(int8_t *h, int8_t *c, int8_t *K)
{
    int8_t m[MESSBYTES+CIPHLEN];
    int8_t hash_m[MESSBYTES];
    int16_t tmp_h[DIM_N], tmp_c[DIM_N];

    randombytes(m, MESSBYTES);
    sha3_256(hash_m, m, MESSBYTES);
    decode_pk(h, tmp_h);
    END_encrypt(tmp_h, m, hash_m, tmp_c);
    encode_c(tmp_c, m + MESSBYTES);
    memcpy(c, m + MESSBYTES, CIPHLEN);
    sha3_256(K, m, MESSBYTES + CIPHLEN);
}

void END_cca_decrypt(int8_t *h, int8_t *c, int8_t *K, int8_t *f, int8_t *g, int8_t *F, int8_t *G, int8_t *k, int16_t *w)
{

    int fail;
    int8_t  m[MESSBYTES + CIPHLEN], r[MESSBYTES];
    int16_t cc[DIM_N], tmp_c[DIM_N], c_re[DIM_N], tmp_h[DIM_N];

    decode_c(c, tmp_c);
    memcpy(cc, tmp_c, DIM_N * sizeof(int16_t));
    decode_pk(h, tmp_h);
    END_decrypt(tmp_h, tmp_c, f, g, F, G, w, m);
    sha3_256(r, m, MESSBYTES);
    END_encrypt(tmp_h, m, r, c_re);
    fail = verify(cc, c_re);
    cmov(m, k, MESSBYTES, fail);
    memcpy(m + MESSBYTES, c, CIPHLEN);
    sha3_256(K, m, MESSBYTES + CIPHLEN);
}
