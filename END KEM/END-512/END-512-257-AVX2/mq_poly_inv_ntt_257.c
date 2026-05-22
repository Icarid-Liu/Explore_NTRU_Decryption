#include <stdint.h>
#include <stdio.h>
#include "param.h"
#include "mq_ntt.h"
#include "mq_ntt_param.h"

__m128i cvtepi32_epi16_avx2(__m256i v) 
{
    __m256i masked = _mm256_and_si256(v, _mm256_set1_epi32(0xFFFF));

    const __m256i shuf_mask = _mm256_setr_epi8(
         0,  1,  4,  5,  8,  9, 12, 13,  0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
         0,  1,  4,  5,  8,  9, 12, 13,  0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80
    );
    __m256i shuf = _mm256_shuffle_epi8(masked, shuf_mask);

    const __m256i low_mask = _mm256_setr_epi8(
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
         0,0,0,0,0,0,0,0,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
         0,0,0,0,0,0,0,0
    );
    __m256i low_bytes = _mm256_and_si256(shuf, low_mask);

    __m256i packed = _mm256_permute2x128_si256(low_bytes, low_bytes, 0x20);
    __m128i result = _mm256_castsi256_si128(packed);

    return result;
}

static int montgomery_mul(__m256i *tmp_c,__m256i *tmp_a,__m256i *tmp_b,__m256i *tmp_Q,__m256i *tmp_QINV)
{
  __m256i tmp_c0,tmp_c1,tmp_x;

  tmp_c0 = _mm256_mullo_epi16(*tmp_a, *tmp_b);
  tmp_c1 = _mm256_mulhi_epi16(*tmp_a, *tmp_b);
  tmp_x = _mm256_mullo_epi16(tmp_c0, *tmp_QINV);
  tmp_x = _mm256_mulhi_epi16(tmp_x, *tmp_Q);
  *tmp_c = _mm256_sub_epi16(tmp_c1, tmp_x);
  return 0;
}

// if f is reversible return 0
int mq_poly_inv_ntt_257(int16_t *f_inv,int16_t *f)
{
  int i=0;
  __m256i tmp_Q=_mm256_set1_epi16(257);
  __m256i tmp_QINV=_mm256_set1_epi16(-255);
  __m256i tmp_0=_mm256_set1_epi16(0);
  __m256i tmp_1=_mm256_set1_epi16(1);
  __m256i tmp_2=_mm256_set1_epi16(2);
  __m256i tmp_4=_mm256_set1_epi16(4);
  __m256i tmp_neg1=_mm256_set1_epi16(-1);
  __m256i tmp_zeta;
  __m256i tmp_a0,tmp_a1,tmp_a2,tmp_a3,tmp_b0,tmp_b1,tmp_b2,tmp_b3,det,tmp_t,tmp_k;
  __m128i det0,det1;
  __m256i det2, det3;
  int bad = 0;
  for(i=0;i<DIM_N;i+=64)
  {
    tmp_zeta=_mm256_load_si256((__m256i *)(inv_zeta_257+i/4));
    tmp_a0=_mm256_load_si256((__m256i *)(f+i));
    tmp_a1=_mm256_load_si256((__m256i *)(f+i+16));
    tmp_a2=_mm256_load_si256((__m256i *)(f+i+32));
    tmp_a3=_mm256_load_si256((__m256i *)(f+i+48));

    montgomery_mul(&tmp_b0,&tmp_a2,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a1,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_2,&tmp_Q,&tmp_QINV);
    tmp_b0=_mm256_add_epi16(tmp_b0,tmp_t);
    montgomery_mul(&tmp_b0,&tmp_b0,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a1,&tmp_a1,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a2,&tmp_Q,&tmp_QINV);
    tmp_b0=_mm256_sub_epi16(tmp_b0,tmp_t);
    montgomery_mul(&tmp_b0,&tmp_b0,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a3,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_zeta,&tmp_Q,&tmp_QINV);
    tmp_b0=_mm256_sub_epi16(tmp_b0,tmp_t);
    montgomery_mul(&tmp_t,&tmp_a0,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a0,&tmp_Q,&tmp_QINV);
    tmp_b0=_mm256_sub_epi16(tmp_b0,tmp_t);

    montgomery_mul(&tmp_b1,&tmp_a1,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a0,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_2,&tmp_Q,&tmp_QINV);
    tmp_b1=_mm256_sub_epi16(tmp_b1,tmp_t);
    montgomery_mul(&tmp_b1,&tmp_b1,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a1,&tmp_a1,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a3,&tmp_Q,&tmp_QINV);
    tmp_b1=_mm256_sub_epi16(tmp_b1,tmp_t);
    montgomery_mul(&tmp_b1,&tmp_b1,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a3,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_zeta,&tmp_Q,&tmp_QINV);
    tmp_b1=_mm256_add_epi16(tmp_b1,tmp_t);
    montgomery_mul(&tmp_t,&tmp_a0,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a1,&tmp_Q,&tmp_QINV);
    tmp_b1=_mm256_add_epi16(tmp_b1,tmp_t);

    montgomery_mul(&tmp_b2,&tmp_a1,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_b2,&tmp_b2,&tmp_2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a2,&tmp_a2,&tmp_Q,&tmp_QINV);
    tmp_b2=_mm256_sub_epi16(tmp_b2,tmp_t);
    montgomery_mul(&tmp_b2,&tmp_b2,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a3,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a0,&tmp_Q,&tmp_QINV);
    tmp_b2=_mm256_sub_epi16(tmp_b2,tmp_t);
    montgomery_mul(&tmp_b2,&tmp_b2,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a0,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a2,&tmp_Q,&tmp_QINV);
    tmp_b2=_mm256_add_epi16(tmp_b2,tmp_t);
    montgomery_mul(&tmp_t,&tmp_a1,&tmp_a1,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a0,&tmp_Q,&tmp_QINV);
    tmp_b2=_mm256_sub_epi16(tmp_b2,tmp_t);

    montgomery_mul(&tmp_b3,&tmp_a2,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a1,&tmp_a3,&tmp_Q,&tmp_QINV);
    tmp_b3=_mm256_sub_epi16(tmp_b3,tmp_t);
    montgomery_mul(&tmp_b3,&tmp_b3,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_b3,&tmp_b3,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a1,&tmp_a1,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a1,&tmp_Q,&tmp_QINV);
    tmp_b3=_mm256_add_epi16(tmp_b3,tmp_t);
    montgomery_mul(&tmp_t,&tmp_a0,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a1,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_2,&tmp_Q,&tmp_QINV);
    tmp_b3=_mm256_sub_epi16(tmp_b3,tmp_t);
    montgomery_mul(&tmp_t,&tmp_a0,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a3,&tmp_Q,&tmp_QINV);
    tmp_b3=_mm256_add_epi16(tmp_b3,tmp_t);

    montgomery_mul(&det,&tmp_a2,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a1,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_4,&tmp_Q,&tmp_QINV);
    det=_mm256_sub_epi16(det,tmp_t);
    montgomery_mul(&det,&det,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&det,&det,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a0,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_k,&tmp_a1,&tmp_a1,&tmp_Q,&tmp_QINV);
    tmp_t=_mm256_add_epi16(tmp_t,tmp_k);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_2,&tmp_Q,&tmp_QINV);
    det=_mm256_add_epi16(tmp_t,det);
    montgomery_mul(&det,&det,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&det,&det,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_a3,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_zeta,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_zeta,&tmp_Q,&tmp_QINV);
    det=_mm256_sub_epi16(tmp_t,det);
    montgomery_mul(&tmp_t,&tmp_a0,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a0,&tmp_Q,&tmp_QINV);
    det=_mm256_sub_epi16(det,tmp_t);
    montgomery_mul(&tmp_t,&tmp_a1,&tmp_a3,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_k,&tmp_a2,&tmp_a2,&tmp_Q,&tmp_QINV);
    tmp_t=_mm256_add_epi16(tmp_t,tmp_k);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_a0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_k,&tmp_a0,&tmp_a2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_k,&tmp_k,&tmp_4,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_k,&tmp_k,&tmp_neg1,&tmp_Q,&tmp_QINV);
    montgomery_mul(&det2,&tmp_a1,&tmp_a1,&tmp_Q,&tmp_QINV);
    tmp_k=_mm256_add_epi16(tmp_k,det2);
    montgomery_mul(&tmp_k,&tmp_k,&tmp_a1,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_k,&tmp_k,&tmp_a1,&tmp_Q,&tmp_QINV);
    tmp_t=_mm256_add_epi16(tmp_t,tmp_k);
    montgomery_mul(&tmp_t,&tmp_t,&tmp_zeta,&tmp_Q,&tmp_QINV);
    det=_mm256_add_epi16(det,tmp_t);
    montgomery_mul(&det,&det,&tmp_1,&tmp_Q,&tmp_QINV);
    det=_mm256_add_epi16(det,_mm256_and_si256(_mm256_srai_epi16(det,15),tmp_Q));
    bad += (_mm256_movemask_epi8(_mm256_cmpeq_epi16(det,tmp_0)) == 0xFFFFFFFF);
  
    det0=_mm256_extracti128_si256(det,0);
    det2=_mm256_cvtepi16_epi32(det0);
    det2=_mm256_i32gather_epi32(q_inv_257,det2,sizeof(int32_t));

    det1=_mm256_extracti128_si256(det,1);
    det3=_mm256_cvtepi16_epi32(det1);
    det3=_mm256_i32gather_epi32(q_inv_257,det3,sizeof(int32_t));
    det2=_mm256_packs_epi32(det2,det3);
    det=_mm256_permute4x64_epi64(det2,0xd8);


    montgomery_mul(&tmp_b0,&det,&tmp_b0,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_b1,&det,&tmp_b1,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_b2,&det,&tmp_b2,&tmp_Q,&tmp_QINV);
    montgomery_mul(&tmp_b3,&det,&tmp_b3,&tmp_Q,&tmp_QINV);

    _mm256_storeu_si256((__m256i *)(f_inv+i),tmp_b0);
    _mm256_storeu_si256((__m256i *)(f_inv+i+16),tmp_b1);
    _mm256_storeu_si256((__m256i *)(f_inv+i+32),tmp_b2);
    _mm256_storeu_si256((__m256i *)(f_inv+i+48),tmp_b3);

  }

  return bad;
}
