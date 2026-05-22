#include <stdint.h>
#include <stdio.h>
#include "param.h"
#include "mq_ntt.h"
#include "mq_ntt_param.h"

static void montgomery_reduce_769(__m256i *r, __m256i *a, __m256i *b, __m256i *b_qinv)
{
    __m256i tmp_c0, tmp_c1;

    __m256i tmp_Q = _mm256_set1_epi16(769);

    tmp_c0 = _mm256_mullo_epi16(*a, *b_qinv);
    tmp_c1 = _mm256_mulhi_epi16(*a, *b);

    tmp_c0  = _mm256_mulhi_epi16(tmp_c0, tmp_Q);
    *r     = _mm256_sub_epi16(tmp_c1, tmp_c0);
}


void mq_poly_adjoint_ntt_769(int16_t *r, int16_t *a)
{
    int i;
    __m256i tmp0, tmp1, tmp2, tmp3;
    __m256i tmp_zeta, tmp_zeta_qinv;
    for(i=0;i<DIM_N;i+=64)
    {
        tmp_zeta_qinv = _mm256_loadu_si256((__m256i *)(adjoint_769 + i/2));
        tmp_zeta = _mm256_loadu_si256((__m256i *)(adjoint_769 + i/2 + 16));

        tmp0 = _mm256_loadu_si256((__m256i *)(a+i));
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 64), tmp3);

        tmp0 = _mm256_loadu_si256((__m256i *)(a+48+i));
        montgomery_reduce_769(&tmp0, &tmp0, &tmp_zeta, &tmp_zeta_qinv);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 48), tmp3);


        tmp0 = _mm256_loadu_si256((__m256i *)(a+32+i));
        montgomery_reduce_769(&tmp0, &tmp0, &tmp_zeta, &tmp_zeta_qinv);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 32), tmp3);
        
        tmp0 = _mm256_loadu_si256((__m256i *)(a+16+i));
        montgomery_reduce_769(&tmp0, &tmp0, &tmp_zeta, &tmp_zeta_qinv);
        tmp1 = _mm256_shufflelo_epi16(tmp0, 0x1B);
        tmp2 = _mm256_shufflehi_epi16(tmp1, 0x1B);
        tmp3 = _mm256_permute4x64_epi64 (tmp2, 0x1B);
        _mm256_storeu_si256((__m256i *)(r + DIM_N - i - 16), tmp3);

    }
    
}


static int ntt_core_12289(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);

    tmp_c0 = _mm256_mullo_epi16(*V, *S);
    tmp_c1 = _mm256_mulhi_epi16(*V, *S);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *V     = _mm256_sub_epi16(tmp_c1, tmp_x);

    *V     = _mm256_add_epi16(*V, _mm256_and_si256(_mm256_srai_epi16(*V, 15), tmp_Q));

    *U     = _mm256_add_epi16(*U, _mm256_and_si256(_mm256_srai_epi16(*U, 15), tmp_Q));

    *A1    = _mm256_add_epi16(*U, *V);
    *A1    = _mm256_sub_epi16(*A1, tmp_Q);

    *A2   = _mm256_sub_epi16(*U, *V);

    return 0;
}

static int ntt_core_lazymod_12289(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;
    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);

    tmp_c0 = _mm256_mullo_epi16(*V, *S);
    tmp_c1 = _mm256_mulhi_epi16(*V, *S);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *V     = _mm256_sub_epi16(tmp_c1, tmp_x);

    *A1    = _mm256_add_epi16(*U, *V);

    *A2   = _mm256_sub_epi16(*U, *V);

    return 0;
}

static int intt_core_12289(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);

    *A1 = _mm256_add_epi16(*U, *V);
    *A1 = _mm256_sub_epi16(*A1, tmp_Q);
    *A1 = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));

    *A2 = _mm256_sub_epi16(*U, *V);
 
    tmp_c0 = _mm256_mullo_epi16(*A2, *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2, *S);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2    = _mm256_sub_epi16(tmp_c1, tmp_x);

    *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}

static int intt_core_lazymod_12289(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);
    *A1 = _mm256_add_epi16(*U, *V);
    *A1 = _mm256_sub_epi16(*A1, tmp_Q);

    *A2 = _mm256_sub_epi16(*U, *V);

    tmp_c0 = _mm256_mullo_epi16(*A2, *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2, *S);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2    = _mm256_sub_epi16(tmp_c1, tmp_x);

    return 0;
}

static int intt_core_final_12289(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);
    __m256i tmp_b = _mm256_set1_epi16(7510);

    *A1 = _mm256_add_epi16(*U, *V);

    tmp_c0 = _mm256_mullo_epi16(*A1, tmp_b);
    tmp_c1 = _mm256_mulhi_epi16(*A1, tmp_b);

    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A1   = _mm256_sub_epi16(tmp_c1, tmp_x);

    // *A1 = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));

    *A2 = _mm256_sub_epi16(*U, *V);

    tmp_c0 = _mm256_mullo_epi16(*A2,  *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2,  *S);

    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2   = _mm256_sub_epi16(tmp_c1, tmp_x);

    // *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}

static int intt_core_final_12289_n(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);
    __m256i tmp_b = _mm256_set1_epi16(7510);

    *A1 = _mm256_add_epi16(*U, *V);

    tmp_c0 = _mm256_mullo_epi16(*A1, tmp_b);
    tmp_c1 = _mm256_mulhi_epi16(*A1, tmp_b);

    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A1   = _mm256_sub_epi16(tmp_c1, tmp_x);

    *A1 = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));

    *A2 = _mm256_sub_epi16(*U, *V);

    tmp_c0 = _mm256_mullo_epi16(*A2,  *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2,  *S);

    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2   = _mm256_sub_epi16(tmp_c1, tmp_x);

    *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}


void mq_poly_ntt_12289(int16_t *a)
{
    int i, j, s;

    __m256i tmp_U, tmp_V, tmp_A1, tmp_A2, tmp_S;

    __m256i tmp_shuffle16_pre = _mm256_set_epi8(0x0f, 0x0e, 0x0b, 0x0a, 0x07, 0x06, 0x03, 0x02, 0x0d, 0x0c, 0x09, 0x08, 0x05, 0x04, 0x01, 0x00,
                                        0x0f, 0x0e, 0x0b, 0x0a, 0x07, 0x06, 0x03, 0x02, 0x0d, 0x0c, 0x09, 0x08, 0x05, 0x04, 0x01, 0x00);


    tmp_S = _mm256_set1_epi16(f_12289[1]);
    for(j=0;j<256;j+=16)
    {

        tmp_U = _mm256_load_si256((__m256i *)(a +j));
        tmp_V = _mm256_load_si256((__m256i *)(a  + 256 +j));

        ntt_core_lazymod_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        _mm256_store_si256((__m256i *)(a +j), tmp_A1);
        _mm256_store_si256((__m256i *)(a  + 256 +j), tmp_A2);
    }
    
    s=0;
    for (i = 0; i < DIM_N; i += 256)
    {
        tmp_S = _mm256_set1_epi16(f_12289[2 +(s++)]);
        for(j=0;j<128;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 128+j));

            ntt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 128+j), tmp_A2);
        }
    }

    s=0;
    for (i = 0; i < DIM_N; i += 128)
    {
        //load root
        tmp_S = _mm256_set1_epi16(f_12289[4 +(s++)]);
        for(j=0;j<64;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 64+j));

            ntt_core_lazymod_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 64 +j), tmp_A2);
        }
    }

    s=0;
    for (i = 0; i < DIM_N; i += 64)
    {
        tmp_S = _mm256_set1_epi16(f_12289[8 +(s++)]);
        for(j=0;j<32;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 32+j));

            ntt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 32+j), tmp_A2);
        }
    }


    s=0;
    for (i = 0; i < DIM_N; i += 32)
    {
        tmp_S = _mm256_set1_epi16(f_12289[16 +(s++)]);

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));

        ntt_core_lazymod_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }
    
    
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(f_12289[32 + s + 1], f_12289[32 + s + 1], f_12289[32 + s + 1], f_12289[32 + s + 1], f_12289[32 + s + 1], f_12289[32 + s + 1], f_12289[32 + s + 1], f_12289[32 + s + 1],
                                 f_12289[32 + s], f_12289[32 + s], f_12289[32 + s], f_12289[32 + s], f_12289[32 + s], f_12289[32 + s], f_12289[32 + s], f_12289[32 + s]);  //256长 16个item，每个长16bit
        s += 2;
        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        ntt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        tmp_S = _mm256_set_epi16(f_12289[64 + s + 3], f_12289[64 + s + 3], f_12289[64 + s + 3], f_12289[64 + s + 3], f_12289[64 + s + 1], f_12289[64 + s + 1], f_12289[64 + s + 1], f_12289[64 + s + 1],
                                 f_12289[64 + s + 2], f_12289[64 + s + 2], f_12289[64 + s + 2], f_12289[64 + s + 2], f_12289[64 + s], f_12289[64 + s], f_12289[64 + s], f_12289[64 + s]);
        s += 4;

        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));

        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        ntt_core_lazymod_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }


    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        tmp_S = _mm256_set_epi16(f_12289[128 + s + 7], f_12289[128 + s + 7], f_12289[128 + s + 3], f_12289[128 + s + 3], f_12289[128 + s + 5], f_12289[128 + s + 5], f_12289[128 + s + 1], f_12289[128 + s + 1],
                                 f_12289[128 + s + 6], f_12289[128 + s + 6], f_12289[128 + s + 2], f_12289[128 + s + 2], f_12289[128 + s + 4], f_12289[128 + s + 4], f_12289[128 + s], f_12289[128 + s]);
        s += 8;

        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));

        tmp_A1 = _mm256_shuffle_epi32(tmp_A1, 0xd8);
        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);

        tmp_A2 = _mm256_shuffle_epi32(tmp_A2, 0xd8);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        ntt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        tmp_S = _mm256_set_epi16(f_12289[256 + s + 15], f_12289[256 + s + 7], f_12289[256 + s + 11], f_12289[256 + s + 3], f_12289[256 + s + 13], f_12289[256 + s + 5], f_12289[256 + s + 9], f_12289[256 + s + 1],
                                 f_12289[256 + s + 14], f_12289[256 + s + 6], f_12289[256 + s + 10], f_12289[256 + s + 2], f_12289[256 + s + 12], f_12289[256 + s + 4], f_12289[256 + s + 8], f_12289[256 + s]);
        s += 16;

        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));

        tmp_A1 = _mm256_shuffle_epi8(tmp_A1, tmp_shuffle16_pre);
        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);

        tmp_A2 = _mm256_shuffle_epi8(tmp_A2, tmp_shuffle16_pre);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        ntt_core_lazymod_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        //store
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

}

void mq_poly_mul_ntt_12289(int16_t *c, int16_t *a, int16_t *b)
{
    int i;

    __m256i tmp_a, tmp_b, tmp_c1, tmp_c0, tmp_x;
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);
    __m256i tmp_Q = _mm256_set1_epi16(12289);

    for (i = 0; i < DIM_N; i += 16)
    {

        tmp_a = _mm256_load_si256((__m256i *)(a + i));
        tmp_b = _mm256_load_si256((__m256i *)(b + i));

        tmp_c0 = _mm256_mullo_epi16(tmp_a, tmp_b);
        tmp_c1 = _mm256_mulhi_epi16(tmp_a, tmp_b);
        tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
        tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
        tmp_c1 = _mm256_sub_epi16(tmp_c1, tmp_x);
        tmp_c1 = _mm256_add_epi16(tmp_c1, _mm256_and_si256(_mm256_srai_epi16(tmp_c1, 15), tmp_Q));

        _mm256_store_si256((__m256i *)(c + i), tmp_c1);
    }

}


void mq_poly_intt_12289(int16_t *a)
{
    int i, j, s;

    __m256i tmp_U, tmp_V, tmp_A1, tmp_A2, tmp_S;
    __m256i tmp_shuffle16_pos = _mm256_set_epi8(0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00,
                                        0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00);

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_12289[256 + s + 15], fn_12289[256 + s + 7], fn_12289[256 + s + 11], fn_12289[256 + s + 3], fn_12289[256 + s + 13], fn_12289[256 + s + 5], fn_12289[256 + s + 9], fn_12289[256 + s + 1],
                                 fn_12289[256 + s + 14], fn_12289[256 + s + 6], fn_12289[256 + s + 10], fn_12289[256 + s + 2], fn_12289[256 + s + 12], fn_12289[256 + s + 4], fn_12289[256 + s + 8], fn_12289[256 + s]);
        s += 16;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));

        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi8(tmp_U, tmp_shuffle16_pos);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi8(tmp_V, tmp_shuffle16_pos);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[128 + s + 7], fn_12289[128 + s + 7], fn_12289[128 + s + 3], fn_12289[128 + s + 3], fn_12289[128 + s + 5], fn_12289[128 + s + 5], fn_12289[128 + s + 1], fn_12289[128 + s + 1],
                                 fn_12289[128 + s + 6], fn_12289[128 + s + 6], fn_12289[128 + s + 2], fn_12289[128 + s + 2], fn_12289[128 + s + 4], fn_12289[128 + s + 4], fn_12289[128 + s], fn_12289[128 + s]);
        s += 8;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
 
        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi32(tmp_U, 0xd8);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi32(tmp_V, 0xd8);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }


    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_12289[64 + s + 3], fn_12289[64 + s + 3], fn_12289[64 + s + 3], fn_12289[64 + s + 3], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1],
                                 fn_12289[64 + s + 2], fn_12289[64 + s + 2], fn_12289[64 + s + 2], fn_12289[64 + s + 2], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s]);
        s += 4;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));

        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);


        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }


    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1],
                                 fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s]);
        s += 2;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set1_epi16(fn_12289[16 + s++]);

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
       
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    s=0;
    for (i = 0; i < DIM_N; i += 64)
    {

        tmp_S = _mm256_set1_epi16(fn_12289[8 +(s++)]);

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 32));

        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 32), tmp_A2);


        tmp_U = _mm256_load_si256((__m256i *)(a + i+16));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 48));

        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i+16), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 48), tmp_A2);
    }

    s=0;
    for (i = 0; i < DIM_N; i += 128)
    {

        tmp_S = _mm256_set1_epi16(fn_12289[4 +(s++)]);
        for(j=0;j<64;j+=16)
        {

            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 64+j));

            intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 64+j), tmp_A2);
        }
    }

    s=0;
    for (i = 0; i < DIM_N; i += 256)
    {

        tmp_S = _mm256_set1_epi16(fn_12289[2 +(s++)]);
        for(j=0;j<128;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 128+j));

            intt_core_lazymod_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 128+j), tmp_A2);
        }
    }


    tmp_S = _mm256_set1_epi16(1966);
    for(j=0;j<256;j+=16)
    {
        tmp_U = _mm256_load_si256((__m256i *)(a +j));
        tmp_V = _mm256_load_si256((__m256i *)(a  + 256 +j));

        intt_core_final_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a +j), tmp_A1);
        _mm256_store_si256((__m256i *)(a  + 256 +j), tmp_A2);
    }

}

void mq_poly_intt_12289_n(int16_t *a)
{
    int i, j, s;

    __m256i tmp_U, tmp_V, tmp_A1, tmp_A2, tmp_S;
    __m256i tmp_shuffle16_pos = _mm256_set_epi8(0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00,
                                        0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00);

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_12289[256 + s + 15], fn_12289[256 + s + 7], fn_12289[256 + s + 11], fn_12289[256 + s + 3], fn_12289[256 + s + 13], fn_12289[256 + s + 5], fn_12289[256 + s + 9], fn_12289[256 + s + 1],
                                 fn_12289[256 + s + 14], fn_12289[256 + s + 6], fn_12289[256 + s + 10], fn_12289[256 + s + 2], fn_12289[256 + s + 12], fn_12289[256 + s + 4], fn_12289[256 + s + 8], fn_12289[256 + s]);
        s += 16;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));

        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi8(tmp_U, tmp_shuffle16_pos);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi8(tmp_V, tmp_shuffle16_pos);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_12289[128 + s + 7], fn_12289[128 + s + 7], fn_12289[128 + s + 3], fn_12289[128 + s + 3], fn_12289[128 + s + 5], fn_12289[128 + s + 5], fn_12289[128 + s + 1], fn_12289[128 + s + 1],
                                 fn_12289[128 + s + 6], fn_12289[128 + s + 6], fn_12289[128 + s + 2], fn_12289[128 + s + 2], fn_12289[128 + s + 4], fn_12289[128 + s + 4], fn_12289[128 + s], fn_12289[128 + s]);
        s += 8;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
 
        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi32(tmp_U, 0xd8);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi32(tmp_V, 0xd8);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }


    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_12289[64 + s + 3], fn_12289[64 + s + 3], fn_12289[64 + s + 3], fn_12289[64 + s + 3], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1], fn_12289[64 + s + 1],
                                 fn_12289[64 + s + 2], fn_12289[64 + s + 2], fn_12289[64 + s + 2], fn_12289[64 + s + 2], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s], fn_12289[64 + s]);
        s += 4;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));

        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);


        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }


    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1], fn_12289[32 + s + 1],
                                 fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s], fn_12289[32 + s]);
        s += 2;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set1_epi16(fn_12289[16 + s++]);

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
       
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    s=0;
    for (i = 0; i < DIM_N; i += 64)
    {

        tmp_S = _mm256_set1_epi16(fn_12289[8 +(s++)]);

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 32));

        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 32), tmp_A2);


        tmp_U = _mm256_load_si256((__m256i *)(a + i+16));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 48));

        intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i+16), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 48), tmp_A2);
    }

    s=0;
    for (i = 0; i < DIM_N; i += 128)
    {

        tmp_S = _mm256_set1_epi16(fn_12289[4 +(s++)]);
        for(j=0;j<64;j+=16)
        {

            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 64+j));

            intt_core_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 64+j), tmp_A2);
        }
    }

    s=0;
    for (i = 0; i < DIM_N; i += 256)
    {

        tmp_S = _mm256_set1_epi16(fn_12289[2 +(s++)]);
        for(j=0;j<128;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 128+j));

            intt_core_lazymod_12289(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 128+j), tmp_A2);
        }
    }


    tmp_S = _mm256_set1_epi16(1966);
    for(j=0;j<256;j+=16)
    {
        tmp_U = _mm256_load_si256((__m256i *)(a +j));
        tmp_V = _mm256_load_si256((__m256i *)(a  + 256 +j));

        intt_core_final_12289_n(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a +j), tmp_A1);
        _mm256_store_si256((__m256i *)(a  + 256 +j), tmp_A2);
    }

}

static int ntt_core_first_level_25601(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(25601);
    __m256i tmp_QINV = _mm256_set1_epi16(-25599);

    tmp_c0 = _mm256_mullo_epi16(*V, *S);
    tmp_c1 = _mm256_mulhi_epi16(*V, *S);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *V     = _mm256_sub_epi16(tmp_c1, tmp_x);

    *V     = _mm256_add_epi16(*V, _mm256_and_si256(_mm256_srai_epi16(*V, 15), tmp_Q));

    *A1    = _mm256_add_epi16(*U, *V);
    *A1    = _mm256_sub_epi16(*A1, tmp_Q);

    *A2   = _mm256_sub_epi16(*U, *V);

    return 0;
}


static int ntt_core_25601(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(25601);
    __m256i tmp_QINV = _mm256_set1_epi16(-25599);

    tmp_c0 = _mm256_mullo_epi16(*V, *S);
    tmp_c1 = _mm256_mulhi_epi16(*V, *S);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *V     = _mm256_sub_epi16(tmp_c1, tmp_x);

    *V     = _mm256_add_epi16(*V, _mm256_and_si256(_mm256_srai_epi16(*V, 15), tmp_Q));

    *U     = _mm256_add_epi16(*U, _mm256_and_si256(_mm256_srai_epi16(*U, 15), tmp_Q));


    *A1    = _mm256_add_epi16(*U, *V);
    *A1    = _mm256_sub_epi16(*A1, tmp_Q);

    *A1    = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));

    *A2   = _mm256_sub_epi16(*U, *V);

    *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}

static int intt_core_25601(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(25601);
    __m256i tmp_QINV = _mm256_set1_epi16(-25599);

    *A1 = _mm256_add_epi16(*U, *V);
    *A1 = _mm256_sub_epi16(*A1, tmp_Q);
    *A1 = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));

    *A2 = _mm256_sub_epi16(*U, *V);
 
    tmp_c0 = _mm256_mullo_epi16(*A2, *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2, *S);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2    = _mm256_sub_epi16(tmp_c1, tmp_x);

    *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}

static int intt_core_final_25601(__m256i *U, __m256i *V, __m256i *S,__m256i *A1,__m256i *A2)
{
    __m256i tmp_c0, tmp_c1, tmp_x;

    __m256i tmp_Q = _mm256_set1_epi16(25601);
    __m256i tmp_QINV = _mm256_set1_epi16(-25599);
    __m256i tmp_b = _mm256_set1_epi16(17081);

    *A1 = _mm256_add_epi16(*U, *V);
    *A1 = _mm256_sub_epi16(*A1, tmp_Q);

    tmp_c0 = _mm256_mullo_epi16(*A1, tmp_b);
    tmp_c1 = _mm256_mulhi_epi16(*A1, tmp_b);

    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A1   = _mm256_sub_epi16(tmp_c1, tmp_x);

    // *A1 = _mm256_add_epi16(*A1, _mm256_and_si256(_mm256_srai_epi16(*A1, 15), tmp_Q));

    *A2 = _mm256_sub_epi16(*U, *V);

    tmp_c0 = _mm256_mullo_epi16(*A2,  *S);
    tmp_c1 = _mm256_mulhi_epi16(*A2,  *S);

    tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    *A2   = _mm256_sub_epi16(tmp_c1, tmp_x);

    // *A2    = _mm256_add_epi16(*A2, _mm256_and_si256(_mm256_srai_epi16(*A2, 15), tmp_Q));

    return 0;
}

void mq_poly_ntt_25601(int16_t *a)
{
    int i, j, s;

    __m256i tmp_U, tmp_V, tmp_A1, tmp_A2, tmp_S;

    __m256i tmp_shuffle16_pre = _mm256_set_epi8(0x0f, 0x0e, 0x0b, 0x0a, 0x07, 0x06, 0x03, 0x02, 0x0d, 0x0c, 0x09, 0x08, 0x05, 0x04, 0x01, 0x00,
                                        0x0f, 0x0e, 0x0b, 0x0a, 0x07, 0x06, 0x03, 0x02, 0x0d, 0x0c, 0x09, 0x08, 0x05, 0x04, 0x01, 0x00);


    tmp_S = _mm256_set1_epi16(f_25601[1]);
    for(j=0;j<256;j+=16)
    {

        tmp_U = _mm256_load_si256((__m256i *)(a +j));
        tmp_V = _mm256_load_si256((__m256i *)(a  + 256 +j));

        ntt_core_first_level_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
        _mm256_store_si256((__m256i *)(a +j), tmp_A1);
        _mm256_store_si256((__m256i *)(a  + 256 +j), tmp_A2);
    }
    
    s=0;
    for (i = 0; i < DIM_N; i += 256)
    {
        tmp_S = _mm256_set1_epi16(f_25601[2 +(s++)]);
        for(j=0;j<128;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 128+j));

            ntt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 128+j), tmp_A2);
        }
    }

    s=0;
    for (i = 0; i < DIM_N; i += 128)
    {
        tmp_S = _mm256_set1_epi16(f_25601[4 +(s++)]);
        for(j=0;j<64;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 64+j));

            ntt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 64 +j), tmp_A2);
        }
    }

    s=0;
    for (i = 0; i < DIM_N; i += 64)
    {
        tmp_S = _mm256_set1_epi16(f_25601[8 +(s++)]);
        for(j=0;j<32;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 32+j));

            ntt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 32+j), tmp_A2);
        }
    }


    s=0;
    for (i = 0; i < DIM_N; i += 32)
    {
        tmp_S = _mm256_set1_epi16(f_25601[16 +(s++)]);

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));

        ntt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }
    
    
    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(f_25601[32 + s + 1], f_25601[32 + s + 1], f_25601[32 + s + 1], f_25601[32 + s + 1], f_25601[32 + s + 1], f_25601[32 + s + 1], f_25601[32 + s + 1], f_25601[32 + s + 1],
                                 f_25601[32 + s], f_25601[32 + s], f_25601[32 + s], f_25601[32 + s], f_25601[32 + s], f_25601[32 + s], f_25601[32 + s], f_25601[32 + s]);
        s += 2;
        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);  //A1低128bit + A2低128bit
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);  //A1高128bit + A2高128bit

        ntt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        tmp_S = _mm256_set_epi16(f_25601[64 + s + 3], f_25601[64 + s + 3], f_25601[64 + s + 3], f_25601[64 + s + 3], f_25601[64 + s + 1], f_25601[64 + s + 1], f_25601[64 + s + 1], f_25601[64 + s + 1],
                                 f_25601[64 + s + 2], f_25601[64 + s + 2], f_25601[64 + s + 2], f_25601[64 + s + 2], f_25601[64 + s], f_25601[64 + s], f_25601[64 + s], f_25601[64 + s]);
        s += 4;

        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));

        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);
        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        ntt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }


    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        tmp_S = _mm256_set_epi16(f_25601[128 + s + 7], f_25601[128 + s + 7], f_25601[128 + s + 3], f_25601[128 + s + 3], f_25601[128 + s + 5], f_25601[128 + s + 5], f_25601[128 + s + 1], f_25601[128 + s + 1],
                                 f_25601[128 + s + 6], f_25601[128 + s + 6], f_25601[128 + s + 2], f_25601[128 + s + 2], f_25601[128 + s + 4], f_25601[128 + s + 4], f_25601[128 + s], f_25601[128 + s]);
        s += 8;

        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));

        tmp_A1 = _mm256_shuffle_epi32(tmp_A1, 0xd8);
        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);

        tmp_A2 = _mm256_shuffle_epi32(tmp_A2, 0xd8);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        ntt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        tmp_S = _mm256_set_epi16(f_25601[256 + s + 15], f_25601[256 + s + 7], f_25601[256 + s + 11], f_25601[256 + s + 3], f_25601[256 + s + 13], f_25601[256 + s + 5], f_25601[256 + s + 9], f_25601[256 + s + 1],
                                 f_25601[256 + s + 14], f_25601[256 + s + 6], f_25601[256 + s + 10], f_25601[256 + s + 2], f_25601[256 + s + 12], f_25601[256 + s + 4], f_25601[256 + s + 8], f_25601[256 + s]);
        s += 16;

        tmp_A1 = _mm256_load_si256((__m256i *)(a + i));
        tmp_A2 = _mm256_load_si256((__m256i *)(a + i + 16));

        tmp_A1 = _mm256_shuffle_epi8(tmp_A1, tmp_shuffle16_pre);
        tmp_A1 = _mm256_permute4x64_epi64(tmp_A1, 0xd8);

        tmp_A2 = _mm256_shuffle_epi8(tmp_A2, tmp_shuffle16_pre);
        tmp_A2 = _mm256_permute4x64_epi64(tmp_A2, 0xd8);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        ntt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

}

void mq_poly_mul_ntt_25601(int16_t *c, int16_t *a, int16_t *b)
{
    int i;

    __m256i tmp_a, tmp_b, tmp_c1, tmp_c0, tmp_x;
    __m256i tmp_Q = _mm256_set1_epi16(25601);
    __m256i tmp_QINV = _mm256_set1_epi16(-25599);

    for (i = 0; i < DIM_N; i += 16)
    {

        tmp_a = _mm256_load_si256((__m256i *)(a + i));
        tmp_b = _mm256_load_si256((__m256i *)(b + i));

        tmp_c0 = _mm256_mullo_epi16(tmp_a, tmp_b);
        tmp_c1 = _mm256_mulhi_epi16(tmp_a, tmp_b);
        tmp_x = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
        tmp_x = _mm256_mulhi_epi16(tmp_x, tmp_Q);
        tmp_c1 = _mm256_sub_epi16(tmp_c1, tmp_x);
        tmp_c1 = _mm256_add_epi16(tmp_c1, _mm256_and_si256(_mm256_srai_epi16(tmp_c1, 15), tmp_Q));

        _mm256_store_si256((__m256i *)(c + i), tmp_c1);
    }

}

void mq_poly_intt_25601(int16_t *a)
{
    int i, j, s;

    __m256i tmp_U, tmp_V, tmp_A1, tmp_A2, tmp_S;
    __m256i tmp_shuffle16_pos = _mm256_set_epi8(0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00,
                                        0x0f, 0x0e, 0x07, 0x06, 0x0d, 0x0c, 0x05, 0x04, 0x0b, 0x0a, 0x03, 0x02, 0x09, 0x08, 0x01, 0x00);

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_25601[256 + s + 15], fn_25601[256 + s + 7], fn_25601[256 + s + 11], fn_25601[256 + s + 3], fn_25601[256 + s + 13], fn_25601[256 + s + 5], fn_25601[256 + s + 9], fn_25601[256 + s + 1],
                                 fn_25601[256 + s + 14], fn_25601[256 + s + 6], fn_25601[256 + s + 10], fn_25601[256 + s + 2], fn_25601[256 + s + 12], fn_25601[256 + s + 4], fn_25601[256 + s + 8], fn_25601[256 + s]);
        s += 16;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));

        intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi8(tmp_U, tmp_shuffle16_pos);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi8(tmp_V, tmp_shuffle16_pos);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {
        //S=f[m+i];
        tmp_S = _mm256_set_epi16(fn_25601[128 + s + 7], fn_25601[128 + s + 7], fn_25601[128 + s + 3], fn_25601[128 + s + 3], fn_25601[128 + s + 5], fn_25601[128 + s + 5], fn_25601[128 + s + 1], fn_25601[128 + s + 1],
                                 fn_25601[128 + s + 6], fn_25601[128 + s + 6], fn_25601[128 + s + 2], fn_25601[128 + s + 2], fn_25601[128 + s + 4], fn_25601[128 + s + 4], fn_25601[128 + s], fn_25601[128 + s]);
        s += 8;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
 
        intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_U = _mm256_shuffle_epi32(tmp_U, 0xd8);

        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);
        tmp_V = _mm256_shuffle_epi32(tmp_V, 0xd8);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }


    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_25601[64 + s + 3], fn_25601[64 + s + 3], fn_25601[64 + s + 3], fn_25601[64 + s + 3], fn_25601[64 + s + 1], fn_25601[64 + s + 1], fn_25601[64 + s + 1], fn_25601[64 + s + 1],
                                 fn_25601[64 + s + 2], fn_25601[64 + s + 2], fn_25601[64 + s + 2], fn_25601[64 + s + 2], fn_25601[64 + s], fn_25601[64 + s], fn_25601[64 + s], fn_25601[64 + s]);
        s += 4;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));

        intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);


        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);
        tmp_U = _mm256_permute4x64_epi64(tmp_U, 0xd8);
        tmp_V = _mm256_permute4x64_epi64(tmp_V, 0xd8);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }


    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set_epi16(fn_25601[32 + s + 1], fn_25601[32 + s + 1], fn_25601[32 + s + 1], fn_25601[32 + s + 1], fn_25601[32 + s + 1], fn_25601[32 + s + 1], fn_25601[32 + s + 1], fn_25601[32 + s + 1],
                                 fn_25601[32 + s], fn_25601[32 + s], fn_25601[32 + s], fn_25601[32 + s], fn_25601[32 + s], fn_25601[32 + s], fn_25601[32 + s], fn_25601[32 + s]);
        s += 2;

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        tmp_U = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x20);
        tmp_V = _mm256_permute2x128_si256(tmp_A1, tmp_A2, 0x31);

        _mm256_store_si256((__m256i *)(a + i), tmp_U);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_V);
    }

    s = 0;
    for (i = 0; i < DIM_N; i += 32)
    {

        tmp_S = _mm256_set1_epi16(fn_25601[16 + s++]);

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 16));
       
        intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);
       
        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 16), tmp_A2);
    }

    s=0;
    for (i = 0; i < DIM_N; i += 64)
    {

        tmp_S = _mm256_set1_epi16(fn_25601[8 +(s++)]);

        tmp_U = _mm256_load_si256((__m256i *)(a + i));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 32));

        intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 32), tmp_A2);


        tmp_U = _mm256_load_si256((__m256i *)(a + i+16));
        tmp_V = _mm256_load_si256((__m256i *)(a + i + 48));

        intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a + i+16), tmp_A1);
        _mm256_store_si256((__m256i *)(a + i + 48), tmp_A2);
    }

    s=0;
    for (i = 0; i < DIM_N; i += 128)
    {

        tmp_S = _mm256_set1_epi16(fn_25601[4 +(s++)]);
        for(j=0;j<64;j+=16)
        {

            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 64+j));

            intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 64+j), tmp_A2);
        }
    }

    s=0;
    for (i = 0; i < DIM_N; i += 256)
    {

        tmp_S = _mm256_set1_epi16(fn_25601[2 +(s++)]);
        for(j=0;j<128;j+=16)
        {
            tmp_U = _mm256_load_si256((__m256i *)(a + i+j));
            tmp_V = _mm256_load_si256((__m256i *)(a + i + 128+j));

            intt_core_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

            _mm256_store_si256((__m256i *)(a + i+j), tmp_A1);
            _mm256_store_si256((__m256i *)(a + i + 128+j), tmp_A2);
        }
    }


    tmp_S = _mm256_set1_epi16(19254);
    for(j=0;j<256;j+=16)
    {
        tmp_U = _mm256_load_si256((__m256i *)(a +j));
        tmp_V = _mm256_load_si256((__m256i *)(a  + 256 +j));

        intt_core_final_25601(&tmp_U,&tmp_V,&tmp_S,&tmp_A1,&tmp_A2);

        _mm256_store_si256((__m256i *)(a +j), tmp_A1);
        _mm256_store_si256((__m256i *)(a  + 256 +j), tmp_A2);
    }

}
