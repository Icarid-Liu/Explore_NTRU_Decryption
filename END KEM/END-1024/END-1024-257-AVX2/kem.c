#include <string.h>
#include <stdio.h>
#include "param.h"
#include "mq_ntt.h"
#include "keygen.h"
#include "pke.h"
#include "poly.h"
#include "mq_ntt_param.h"
#include "code.h"
#include "randombytes.h"
#include "fips202.h"

static inline int verify(int16_t *a, int16_t *b)
{
  int i;
  uint64_t r;
  __m256i vdata0, vdata1, vdata2;

  vdata0 = _mm256_setzero_si256();
  for(i = 0; i < DIM_N; i += 16) 
  {
    vdata1 = _mm256_loadu_si256((__m256i *)(a + i));
    vdata2 = _mm256_loadu_si256((__m256i *)(b + i));
    vdata1 = _mm256_xor_si256(vdata1, vdata2);
    vdata0 = _mm256_or_si256(vdata0, vdata1);
  }
  r = 1 - _mm256_testz_si256(vdata0, vdata0);

  r = (-r) >> 63;
  return r;
}

static inline void cmov(int8_t *r, int8_t *x, int len, int b)
{
  int i;
  __m256i vdata0, vdata1, vdata2;

  b = -b;
  vdata0 = _mm256_set1_epi8(b);

  for(i = 0; i + 32 <= len; i += 32) 
  {
    vdata1 = _mm256_loadu_si256((__m256i *)(r + i));
    vdata2 = _mm256_loadu_si256((__m256i *)(x + i));
    vdata2 = _mm256_xor_si256(vdata2, vdata1);
    vdata2 = _mm256_and_si256(vdata2, vdata0);
    vdata1 = _mm256_xor_si256(vdata1, vdata2);
    _mm256_storeu_si256((__m256i *)(r + i), vdata1);
  }

  while(i < len) {
    r[i] ^= b & (x[i] ^ r[i]);
    i += 1;
  }
}

void END_cca_encrypt(int8_t *h, int8_t *c, int8_t *K)
{
    int i,j;
    int8_t  __attribute__((aligned(32))) m[MESSBYTES + CIPHLEN], hash_m[MESSBYTES];
    int16_t __attribute__((aligned(32))) tmp_c[DIM_N], tmp_h[DIM_N];

    randombytes(m, MESSBYTES);
    sha3_512(hash_m, m, MESSBYTES);
    decode_pk(h, tmp_h);
    END_encrypt(tmp_h, m, hash_m, tmp_c);
    encode_c(tmp_c, m + MESSBYTES);
    memcpy(c, m + MESSBYTES, CIPHLEN);
    sha3_256(K, m, MESSBYTES + CIPHLEN);
}

int END_cca_decrypt(int8_t *h, int8_t *c, int8_t *K, int8_t *f, int8_t *g, int8_t *F, int8_t *G, int8_t *k, int16_t *w)
{

    int i, j, fail;
    int8_t  __attribute__((aligned(32))) m[MESSBYTES + CIPHLEN], r[MESSBYTES];
    int16_t __attribute__((aligned(32))) cc[DIM_N], tmp_c[DIM_N], c_re[DIM_N], tmp_h[DIM_N];

    decode_c(c, tmp_c);
    memcpy(cc, tmp_c, DIM_N * sizeof(int16_t));
    decode_pk(h, tmp_h);
    END_decrypt(tmp_h, tmp_c, f, g, F, G, w, m);
    sha3_512(r, m, MESSBYTES);
    END_encrypt(tmp_h, m, r, c_re);
    fail = verify(cc, c_re);
    cmov(m, k, MESSBYTES, fail);
    memcpy(m + MESSBYTES, c, CIPHLEN);
    sha3_256(K, m, MESSBYTES + CIPHLEN);
}