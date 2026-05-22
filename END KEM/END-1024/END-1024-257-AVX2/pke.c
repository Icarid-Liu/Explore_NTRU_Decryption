#include <string.h>
#include <stdio.h>
#include <immintrin.h>
#include "param.h"
#include "mq_ntt.h"
#include "poly.h"
#include "mq_ntt_param.h"
#include "fips202.h"

static int16_t barrett_reduce_514(int32_t a) 
{
  int32_t t;
  t = ((int32_t)2040 * a);
  t = t >> 20;
  t *= 514;
  return a - t;
}

static int32_t barrett_reduce_6316546(int64_t a) 
{
    int64_t t;
    t = ((int64_t)11140383 * a);
    t = t >> 46; 
    t *= 6316546;
    return a - t;
}

static int16_t barrett_reduce_769(int32_t a) 
{
    int64_t t;
    t = ((int64_t)349070 * a);
    t = t >> 28;  
    t *= 769;
    return (int16_t)(a - t);
}

static int16_t div_6316546(int32_t a) 
{
    int64_t t;
    t = ((int64_t)11140383 * a);
    t = t >> 46; 
    t = (int16_t)t;
    t += (t>>15) & 1;
    return (int16_t)t;
}

static __m256i montgomery_reduce_257(__m256i a, __m256i b)
{
    __m256i tmp_c0, tmp_c1, tmp_x, r;

    __m256i tmp_Q = _mm256_set1_epi16(257);
    __m256i tmp_QINV = _mm256_set1_epi16(-255);

    tmp_c0 = _mm256_mullo_epi16(a, b);
    tmp_c1 = _mm256_mulhi_epi16(a, b);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    r     = _mm256_sub_epi16(tmp_c1, tmp_x);

    return r;
}

static __m256i montgomery_reduce_769(__m256i a, __m256i b)
{
    __m256i tmp_c0, tmp_c1, tmp_x, r;

    __m256i tmp_Q = _mm256_set1_epi16(769);
    __m256i tmp_QINV = _mm256_set1_epi16(-767);

    tmp_c0 = _mm256_mullo_epi16(a, b);
    tmp_c1 = _mm256_mulhi_epi16(a, b);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    r = _mm256_sub_epi16(tmp_c1, tmp_x);

    return r;
}

static __m256i montgomery_reduce_12289(__m256i a, __m256i b)
{
    __m256i tmp_c0, tmp_c1, tmp_x, r;

    __m256i tmp_Q = _mm256_set1_epi16(12289);
    __m256i tmp_QINV = _mm256_set1_epi16(-12287);

    tmp_c0 = _mm256_mullo_epi16(a, b);
    tmp_c1 = _mm256_mulhi_epi16(a, b);

    tmp_x  = _mm256_mullo_epi16(tmp_c0, tmp_QINV);
    tmp_x  = _mm256_mulhi_epi16(tmp_x, tmp_Q);
    r = _mm256_sub_epi16(tmp_c1, tmp_x);

    return r;
}

int END_encrypt(int16_t *h, int8_t *m, int8_t *r, int16_t *c)
{
    int8_t __attribute__((aligned(32))) m_prime[MESSBYTES];
    int16_t __attribute__((aligned(32))) s[DIM_N];
    int16_t i, j;
    __m256i vmax, vdata0, vdata1, vone32 = _mm256_set1_epi8(1); 
    __m128i vdata2, vdata3;

    sha3_512(m_prime, r, MESSBYTES);

    for(j = 0; j < 64; j += 32)
    {
        vdata0 = _mm256_load_si256((__m256i *)(r + j));
        for(i = 0; i < 8; i++)
        {
            vdata1 = _mm256_and_si256(vdata0, vone32);
            vdata0 = _mm256_srli_epi16(vdata0, 1);
            vdata2 = _mm256_extracti128_si256(vdata1, 0);
            vdata3 = _mm256_extracti128_si256(vdata1, 1);
            vdata1 = _mm256_cvtepi8_epi16(vdata2);
            _mm256_store_si256((__m256i *)(s + 512 + j * 8 + i * 32), vdata1);
            vdata1 = _mm256_cvtepi8_epi16(vdata3);
            _mm256_store_si256((__m256i *)(s + 512 + j * 8 + i * 32 + 16), vdata1);
        }
    }

    for(j = 0; j < 64; j += 32)
    {
        vdata0 = _mm256_load_si256((__m256i *)(m + j));
        vdata1 = _mm256_load_si256((__m256i *)(m_prime + j));
        vdata0 = _mm256_xor_si256(vdata0, vdata1);
        for(i = 0; i < 8; i++)
        {
            vdata1 = _mm256_and_si256(vdata0, vone32);
            vdata0 = _mm256_srli_epi16(vdata0, 1);
            vdata2 = _mm256_extracti128_si256(vdata1, 0);
            vdata3 = _mm256_extracti128_si256(vdata1, 1);
            vdata1 = _mm256_cvtepi8_epi16(vdata2);
            _mm256_store_si256((__m256i *)(s + j * 8 + i * 32), vdata1);
            vdata1 = _mm256_cvtepi8_epi16(vdata3);
            _mm256_store_si256((__m256i *)(s + j * 8 + i * 32 + 16), vdata1);
        }
    }
    //check e
    mq_poly_ntt_257_769(s, nttdata_257);
    mq_poly_mul_ntt_257(c, h, s, muldata_257);
    mq_poly_intt_257_769_n(c, inttdata_257);

    poly_round(c);

    return 0;
}

void END_decrypt(int16_t *h, int16_t *c, int8_t *f, int8_t *g, int8_t *F, int8_t *G, int16_t *w, int8_t *m)
{
    int16_t i, tmp, u1, idx0, idx1, idx2;
    uint32_t norm1, norm2, norm3, norm4;
    uint32_t mask1, mask2, mask3;
    int32_t t, num1, num2;
    int32_t __attribute__((aligned(32))) C2[DIM_N], Tmp0[DIM_N], Tmp1[DIM_N];
    int16_t __attribute__((aligned(32))) tmp1[DIM_N],tmp2[DIM_N], tmp3[DIM_N], tmp4[DIM_N], tmp5[DIM_N], tmp6[DIM_N], tmp7[DIM_N];
    int16_t __attribute__((aligned(32))) c1[DIM_N], c2[DIM_N], sp[DIM_N], ep[DIM_N], Fp[DIM_N], Gp[DIM_N];
    int8_t __attribute__((aligned(32))) r[MESSBYTES], r_hash[MESSBYTES];

    __m256i vtmp1, vtmp2, vtmp3, vtmp4, vtmp5, vtmp6, vtmp7, vone32 = _mm256_set1_epi8(1);
    __m128i vtmp0, vtmp8;

    poly_deround(c);

    vtmp1 = _mm256_set1_epi16(1);
    for(i = 0; i < DIM_N; i += 16)
    { 
        _mm256_store_si256((__m256i *)(tmp1 + i), vtmp1);
    }

    poly_int8_to_int16(f, tmp2); //tmp2 = f
    poly_int8_to_int16(g, tmp3); //tmp3 = g
    poly_int8_to_int16(F, tmp4); //tmp4 = F
    poly_int8_to_int16(G, tmp5); //tmp5 = G
    norm2 = poly_norm2(tmp4) + poly_norm2(tmp5);

    mq_poly_ntt_257_769(tmp1, nttdata_257); //tmp1 = vector(1)_769
    memcpy(tmp6, c, DIM_N * sizeof(int16_t)); // tmp6 = c
    mq_poly_ntt_257_769(tmp6, nttdata_257); //tmp6 = c_769

    mq_poly_ntt_257_769(tmp2, nttdata_257); //tmp2 = f_769
    mq_poly_ntt_257_769(tmp3, nttdata_257); //tmp3 = g_769
    mq_poly_mul_ntt_257(c1, tmp6, tmp2, muldata_257); //c1 = cf_769

    mq_poly_mul_ntt_257(tmp2, tmp1, tmp3, muldata_257); //tmp2 = g.1_769

    vtmp1 = _mm256_set1_epi16(2);
    vtmp4 = _mm256_set1_epi16(1);
    for(i = 0; i < DIM_N; i += 16)
    { 
	    vtmp2 = _mm256_load_si256((__m256i *)(c1 + i));
        vtmp2 = _mm256_mullo_epi16(vtmp2, vtmp1);
        vtmp3 = _mm256_load_si256((__m256i *)(tmp2 + i));
        vtmp2 = _mm256_sub_epi16(vtmp2, vtmp3);
        vtmp3 = montgomery_reduce_257(vtmp2, vtmp4);
        _mm256_store_si256((__m256i *)(c1 + i), vtmp3);
    }

    mask1 = 0;
    vtmp1 = _mm256_set1_epi8(1);
    vtmp2 = _mm256_setzero_si256();
    for(i = 0; i < DIM_N; i += 32)
    { 
	    vtmp3 = _mm256_load_si256((__m256i *)(g + i));
        vtmp3 = _mm256_and_si256(vtmp3, vtmp1);
        vtmp2 = _mm256_xor_si256(vtmp2, vtmp3);
    }
    _mm256_store_si256((__m256i *)(r), vtmp2);
    for(i = 0; i < 32; i++)
    {
        mask1 ^= r[i];
    }
    mq_poly_intt_257_769_n(c1, inttdata_257); //c1 = cf-g.1_mod257

    for(i = 0; i < DIM_N; i++)
    {
        idx0 = barrett_reduce_514(c1[i] * 258 - 257);
        idx1 = barrett_reduce_514(c1[i] * 258);
        idx2 = (idx0 & (-mask1)) | (idx1 & (~(-mask1)));
        u1 = idx2 > 256;
        c1[i] = idx2 - (u1 * 514);
    }

    mq_poly_ntt_257_769(tmp4, nttdata_257); //tmp4 = F_769
    mq_poly_ntt_257_769(tmp5, nttdata_257); //tmp5 = G_769
    mq_poly_mul_ntt_257(c2, tmp6, tmp4, muldata_257); //c2 = cF_769
    mq_poly_mul_ntt_257(tmp4, tmp1, tmp5, muldata_257); //tmp4 = G.1_769

    vtmp1 = _mm256_set1_epi16(2);
    vtmp4 = _mm256_set1_epi16(1);
    for(i = 0; i < DIM_N; i += 16)
    { 
	    vtmp2 = _mm256_load_si256((__m256i *)(c2 + i));
        vtmp2 = _mm256_mullo_epi16(vtmp2, vtmp1);
        vtmp3 = _mm256_load_si256((__m256i *)(tmp4 + i));
        vtmp2 = _mm256_sub_epi16(vtmp2, vtmp3);
        vtmp3 = montgomery_reduce_257(vtmp2, vtmp4);
        _mm256_store_si256((__m256i *)(c2 + i), vtmp3);
    }

    mask1 = 0;
    vtmp1 = _mm256_set1_epi8(1);
    vtmp2 = _mm256_setzero_si256();
    for(i = 0; i < DIM_N; i += 32)
    { 
	    vtmp3 = _mm256_load_si256((__m256i *)(G + i));
        vtmp3 = _mm256_and_si256(vtmp3, vtmp1);
        vtmp2 = _mm256_xor_si256(vtmp2, vtmp3);
    }
    _mm256_store_si256((__m256i *)(r), vtmp2);
    for(i = 0; i < 32; i++)
    {
        mask1 ^= r[i];
    }
    mq_poly_intt_257_769_n(c2, inttdata_257); //c2 = cF-G.1_mod257

    for(i = 0; i < DIM_N; i++)
    {
        idx0 = barrett_reduce_514(c2[i] * 258 - 257);
        idx1 = barrett_reduce_514(c2[i] * 258);
        idx2 = (idx0 & (-mask1)) | (idx1 & (~(-mask1)));
        u1 = idx2 > 256;
        c2[i] = idx2 - (u1 * 514);
    }

    vtmp1 = _mm256_set1_epi16(1);
    vtmp2 = _mm256_set1_epi16(4091);
    for(i = 0; i < DIM_N; i += 16)
    {
        vtmp3 = _mm256_load_si256((__m256i *)(w + i));
        vtmp4 = montgomery_reduce_257(vtmp3, vtmp1);
        vtmp5 = montgomery_reduce_12289(vtmp3, vtmp2);
        _mm256_store_si256((__m256i *)(tmp2 + i), vtmp4);
        _mm256_store_si256((__m256i *)(tmp3 + i), vtmp5);
    }

    memcpy(tmp1, c1, DIM_N * sizeof(int16_t)); // tmp1 = c1
    mq_poly_ntt_257_769(tmp1, nttdata_257); //tmp1 = c_257
    mq_poly_ntt_257_769(tmp2, nttdata_257); //tmp2 = w_257
    mq_poly_mul_ntt_257(tmp4, tmp1, tmp2, muldata_257); //tmp4 = c1w_257
    mq_poly_intt_257_769_n(tmp4, inttdata_257); // tmp4 = c1w % 257

    memcpy(tmp1, c1, DIM_N * sizeof(int16_t)); // tmp1 = c1
    mq_poly_ntt_12289(tmp1); // tmp1 = c1_12289
    mq_poly_ntt_12289(tmp3); // tmp3 = w_12289
    mq_poly_mul_ntt_12289(tmp2, tmp1, tmp3); // tmp2 = c1w_12289
    mq_poly_intt_12289_n(tmp2); // tmp2 = c1w % 12289

    idx0 = 0, idx1 = 0;
    for(i = 0; i < DIM_N; i++)
    {
        idx0 ^= (c1[i] & 1);
        idx1 ^= (w[i] & 1);
    }
    mask1 = idx0 ^ idx1;

    for(i = 0; i < DIM_N; i++)
    {
        int32_t t2 = barrett_reduce_6316546((12289 * c2[i]) - (((int64_t)(tmp4[i]) * 1007698) + ((int64_t)(tmp2[i]) * 2150576)) + (63165460000LL));
        int32_t t1 = barrett_reduce_6316546(t2 + 3158273);
        int32_t val = (t1 & (-mask1)) | (t2 & (~(-mask1)));
        c2[i] = barrett_reduce_769(val);
        int8_t over = val > 3158273;
        c2[i] = c2[i] + (over * 20);
        c2[i] = c2[i] % 769;
        val = val - (over * 6316546);
        C2[i] = val;
    }

    poly_int8_to_int16(f, tmp1); //tmp1 = f
    mq_poly_ntt_257_769(tmp1, nttdata_769); //tmp1 = f_769
    mq_poly_ntt_257_769(c2, nttdata_769); //c2 = c2_769
    mq_poly_mul_ntt_769(tmp2, tmp1, c2, muldata_769); //tmp2 = c2f_769

    vtmp1 = _mm256_set1_epi16(171);
    for(i = 0; i < DIM_N; i += 16)
    {
        vtmp2 = _mm256_load_si256((__m256i *)(w + i));
        vtmp3 = montgomery_reduce_769(vtmp2, vtmp1);
        _mm256_store_si256((__m256i *)(tmp3 + i), vtmp3);
    }//tmp3 = w_769;

    mq_poly_ntt_257_769(tmp3, nttdata_769); //tmp3 = w_769
    mq_poly_mul_ntt_769(tmp4, tmp1, tmp3, muldata_769); //tmp4 = wf_769
    mq_poly_intt_257_769(tmp4, inttdata_769); //tmp4 = wf

    poly_int8_to_int16(F, tmp3); //tmp3 = F

    vtmp1 = _mm256_set1_epi16(511);
    vtmp2 = _mm256_set1_epi16(171);
    for(i = 0; i < DIM_N; i += 16)
    {
        vtmp3 = _mm256_load_si256((__m256i *)(tmp3 + i)); //F
        vtmp4 = montgomery_reduce_769(vtmp3, vtmp1); //QF
        vtmp5 = _mm256_load_si256((__m256i *)(tmp4 + i)); //wf
        vtmp3 = _mm256_sub_epi16(vtmp4, vtmp5); //QF-wf
        vtmp4 = montgomery_reduce_769(vtmp3, vtmp2);
        _mm256_store_si256((__m256i *)(tmp3 + i), vtmp4);
    }//tmp3 = QF-wf

    mq_poly_ntt_257_769(c1, nttdata_769); //c1 = c1_769
    mq_poly_ntt_257_769(tmp3, nttdata_769); //tmp3 = QF-wf_769
    mq_poly_mul_ntt_769(tmp5, c1, tmp3, muldata_769); //tmp5 = c1(QF-wf)_769

    vtmp1 = _mm256_set1_epi16(675);
    for(i = 0; i < DIM_N; i += 16)
    {
        vtmp2 = _mm256_load_si256((__m256i *)(tmp5 + i));
        vtmp3 = _mm256_load_si256((__m256i *)(tmp2 + i));
        vtmp2 = _mm256_sub_epi16(vtmp2, vtmp3);
        vtmp3 = montgomery_reduce_769(vtmp2, vtmp1);
        _mm256_store_si256((__m256i *)(sp + i), vtmp3);
    }//sp = sp_257

    mq_poly_intt_257_769(sp, inttdata_769); //sp = sp

    vtmp1 = _mm256_set1_epi16(1);
    for(i = 0; i < DIM_N; i += 16)
    { 
	    vtmp2 = _mm256_load_si256((__m256i *)(sp + i));
        vtmp2 = _mm256_add_epi16(vtmp2, vtmp1);
        vtmp2 = _mm256_srli_epi16(vtmp2, 1);
        _mm256_store_si256((__m256i *)(sp + i), vtmp2);
    }

    memcpy(tmp3, sp, DIM_N * sizeof(int16_t)); // tmp3 = sp
    mq_poly_ntt_257_769(tmp3, nttdata_257); //tmp3 = sp_257
    mq_poly_mul_ntt_257(tmp2, h, tmp3, muldata_257); // tmp2 = hsp_257
    mq_poly_intt_257_769(tmp2, inttdata_257); //tmp2 = hsp

    vtmp1 = _mm256_set1_epi16(1);
    for(i = 0; i < DIM_N; i += 16)
    {
        vtmp2 = _mm256_load_si256((__m256i *)(c + i));
        vtmp3 = _mm256_load_si256((__m256i *)(tmp2 + i));
        vtmp2 = _mm256_sub_epi16(vtmp2, vtmp3);
        vtmp3 = montgomery_reduce_257(vtmp2, vtmp1);
        _mm256_store_si256((__m256i *)(ep + i), vtmp3);
    } //ep = c - hsp

    poly_int8_to_int16(F, Fp); //Fp = F

    for(i = 0; i < DIM_N; i++)
    {
        tmp3[i] = 128 - abs(c1[i]); //tmp3 = q/2 - abs(c1)
    }

    idx0 = 0, idx1 = 0;
    norm1 = 0xFFFF; norm2 = 0xFFFF;
    for(i = 0; i < DIM_N; i++)
    {
        mask1 = abs(tmp3[i]) < abs(norm1);
        mask2 = (abs(tmp3[i]) < abs(norm2)) & !mask1;

        norm2 = (norm1 & (-mask1)) | (norm2 & (~(-mask1)));
        idx1 = (idx0 & (-mask1)) | (idx1 & (~(-mask1)));

        norm1 = (tmp3[i] & (-mask1)) | (norm1 & (~(-mask1)));
        idx0 = (i & (-mask1)) | (idx0 & (~(-mask1)));

        norm2 = (tmp3[i] & (-mask2)) | (norm2 & (~(-mask2)));
        idx1 = (i & (-mask2)) | (idx1 & (~(-mask2)));
    }

    //Tmp0 = 2qw*x^idx0
    for(i = 0; i < DIM_N; i++)
    {
        mask1 = (i < idx0);
        Tmp0[i] = ((-w[(DIM_N - abs(idx0 - i)) & ((1<<LOGN)-1)] * 514) & (-mask1)) | ((w[abs(i - idx0)] * 514) & (~(-mask1)));
    }

    //Tmp1 = 2qw*x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        mask1 = (i < idx1);
        Tmp1[i] = ((-w[(DIM_N - abs(idx1 - i)) & ((1<<LOGN)-1)] * 514) & (-mask1)) | ((w[abs(i - idx1)] * 514) & (~(-mask1)));
    }
    
    //tmp4 = F*x^idx0
    for(i = 0; i < DIM_N; i++)
    {
        mask1 = (i < idx0);
        tmp4[i] = (-Fp[(DIM_N - abs(idx0 - i)) & ((1<<LOGN)-1)] & (-mask1)) | (Fp[abs(i - idx0)] & (~(-mask1)));
    }
    
    //tmp5 = F*x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        mask1 = (i < idx1);
        tmp5[i] = (-Fp[(DIM_N - abs(idx1 - i)) & ((1<<LOGN)-1)] & (-mask1)) | (Fp[abs(i - idx1)] & (~(-mask1)));
    }
    
    //0
    norm1 = poly_norm2(sp);


    //x^idx0
    for(i = 0; i < DIM_N; i++)
    {
        norm3 = div_6316546((C2[i] - Tmp0[i] - 3158273));
        norm4 = div_6316546((C2[i] - Tmp0[i] + 3158273));
        mask1 = (C2[i] - Tmp0[i]) < 0;
        tmp6[i] = (norm3 & (-mask1)) | (norm4 & (~(-mask1)));
    } //tmp6 = C2-2qw^idx0/2qQ
    mq_poly_ntt_257_769(tmp6, nttdata_257); //tmp2 = C2-2qw^idx0/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1, muldata_257); //tmp4 = C2-qwe1/qQ * f_257
    mq_poly_intt_257_769(tmp7, inttdata_257); //tmp7 = C2-2qw^idx0/2qQ * f
    for(i = 0; i < DIM_N; i++)
    {
        tmp7[i] = sp[i] - tmp7[i] - tmp4[i];
    }
    norm2 = poly_norm2(tmp7);
    mask1 = norm1 > norm2;
    norm1 = (norm2 & (-mask1)) | (norm1 & (~(-mask1)));
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp7[i] & (-mask1)) | (sp[i] & (~(-mask1)));
        C2[i] = (C2[i] - Tmp0[i] & (-mask1)) | (C2[i] & (~(-mask1)));
    }

    //-x^idx0
    for(i = 0; i < DIM_N; i++)
    {
        norm3 = div_6316546((C2[i] - (-Tmp0[i]) - 3158273));
        norm4 = div_6316546((C2[i] - (-Tmp0[i]) + 3158273));
        mask1 = (C2[i] - (-Tmp0[i])) < 0;
        tmp6[i] = (norm3 & (-mask1)) | (norm4 & (~(-mask1)));
    } //tmp6 = C2-2qw^idx0/2qQ
    mq_poly_ntt_257_769(tmp6, nttdata_257); //tmp2 = C2-2qw^idx0/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1, muldata_257); //tmp4 = C2-qwe1/qQ * f_257
    mq_poly_intt_257_769(tmp7, inttdata_257); //tmp7 = C2-2qw^idx0/2qQ * f
    for(i = 0; i < DIM_N; i++)
    {
        tmp7[i] = sp[i] - tmp7[i] - (-tmp4[i]);
    }
    norm2 = poly_norm2(tmp7);
    mask1 = norm1 > norm2;
    norm1 = (norm2 & (-mask1)) | (norm1 & (~(-mask1)));
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp7[i] & (-mask1)) | (sp[i] & (~(-mask1)));
        C2[i] = (C2[i] - (-Tmp0[i]) & (-mask1)) | (C2[i] & (~(-mask1)));
    }

    //x^idx0 + x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        norm3 = div_6316546((C2[i] - (Tmp0[i] + Tmp1[i]) - 3158273));
        norm4 = div_6316546((C2[i] - (Tmp0[i] + Tmp1[i]) + 3158273));
        mask1 = (C2[i] - (Tmp0[i] + Tmp1[i])) < 0;
        tmp6[i] = (norm3 & (-mask1)) | (norm4 & (~(-mask1)));
    } //tmp6 = C2-2q(w^idx0+w^idx1)/2qQ
    mq_poly_ntt_257_769(tmp6, nttdata_257); //tmp2 = C2-2qw^idx0/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1, muldata_257); //tmp4 = C2-qwe1/qQ * f_257
    mq_poly_intt_257_769(tmp7, inttdata_257); //tmp7 = C2-2qw^idx0/2qQ * f
    for(i = 0; i < DIM_N; i++)
    {
        tmp7[i] = sp[i] - tmp7[i] - (tmp4[i] + tmp5[i]);
    }
    norm2 = poly_norm2(tmp7);
    mask1 = norm1 > norm2;
    norm1 = (norm2 & (-mask1)) | (norm1 & (~(-mask1)));
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp7[i] & (-mask1)) | (sp[i] & (~(-mask1)));
        C2[i] = (C2[i] - (Tmp0[i] + Tmp1[i]) & (-mask1)) | (C2[i] & (~(-mask1)));
    }

    //x^idx0 - x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        norm3 = div_6316546((C2[i] - (Tmp0[i] - Tmp1[i]) - 3158273));
        norm4 = div_6316546((C2[i] - (Tmp0[i] - Tmp1[i]) + 3158273));
        mask1 = (C2[i] - (Tmp0[i] - Tmp1[i])) < 0;
        tmp6[i] = (norm3 & (-mask1)) | (norm4 & (~(-mask1)));
    } //tmp6 = C2-2q(w^idx0+w^idx1)/2qQ
    mq_poly_ntt_257_769(tmp6, nttdata_257); //tmp2 = C2-2qw^idx0/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1, muldata_257); //tmp4 = C2-qwe1/qQ * f_257
    mq_poly_intt_257_769(tmp7, inttdata_257); //tmp7 = C2-2qw^idx0/2qQ * f
    for(i = 0; i < DIM_N; i++)
    {
        tmp7[i] = sp[i] - tmp7[i] - (tmp4[i] - tmp5[i]);
    }
    norm2 = poly_norm2(tmp7);
    mask1 = norm1 > norm2;
    norm1 = (norm2 & (-mask1)) | (norm1 & (~(-mask1)));
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp7[i] & (-mask1)) | (sp[i] & (~(-mask1)));
        C2[i] = (C2[i] - (Tmp0[i] - Tmp1[i]) & (-mask1)) | (C2[i] & (~(-mask1)));
    }

    //-x^idx0 - x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        norm3 = div_6316546((C2[i] - (-Tmp0[i] - Tmp1[i]) - 3158273));
        norm4 = div_6316546((C2[i] - (-Tmp0[i] - Tmp1[i]) + 3158273));
        mask1 = (C2[i] - (-Tmp0[i] - Tmp1[i])) < 0;
        tmp6[i] = (norm3 & (-mask1)) | (norm4 & (~(-mask1)));
    } //tmp6 = C2-2q(w^idx0+w^idx1)/2qQ
    mq_poly_ntt_257_769(tmp6, nttdata_257); //tmp2 = C2-2qw^idx0/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1, muldata_257); //tmp4 = C2-qwe1/qQ * f_257
    mq_poly_intt_257_769(tmp7, inttdata_257); //tmp7 = C2-2qw^idx0/2qQ * f
    for(i = 0; i < DIM_N; i++)
    {
        tmp7[i] = sp[i] - tmp7[i] - (-tmp4[i] - tmp5[i]);
    }
    norm2 = poly_norm2(tmp7);
    mask1 = norm1 > norm2;
    norm1 = (norm2 & (-mask1)) | (norm1 & (~(-mask1)));
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp7[i] & (-mask1)) | (sp[i] & (~(-mask1)));
        C2[i] = (C2[i] - (-Tmp0[i] - Tmp1[i]) & (-mask1)) | (C2[i] & (~(-mask1)));
    }

    //-x^idx0 + x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        norm3 = div_6316546((C2[i] - (-Tmp0[i] + Tmp1[i]) - 3158273));
        norm4 = div_6316546((C2[i] - (-Tmp0[i] + Tmp1[i]) + 3158273));
        mask1 = (C2[i] - (-Tmp0[i] + Tmp1[i])) < 0;
        tmp6[i] = (norm3 & (-mask1)) | (norm4 & (~(-mask1)));
    } //tmp6 = C2-2q(w^idx0+w^idx1)/2qQ
    mq_poly_ntt_257_769(tmp6, nttdata_257); //tmp2 = C2-2qw^idx0/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1, muldata_257); //tmp4 = C2-qwe1/qQ * f_257
    mq_poly_intt_257_769(tmp7, inttdata_257); //tmp7 = C2-2qw^idx0/2qQ * f
    for(i = 0; i < DIM_N; i++)
    {
        tmp7[i] = sp[i] - tmp7[i] - (-tmp4[i] + tmp5[i]);
    }
    norm2 = poly_norm2(tmp7);
    mask1 = norm1 > norm2;
    norm1 = (norm2 & (-mask1)) | (norm1 & (~(-mask1)));
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp7[i] & (-mask1)) | (sp[i] & (~(-mask1)));
        C2[i] = (C2[i] - (-Tmp0[i] + Tmp1[i]) & (-mask1)) | (C2[i] & (~(-mask1)));
    }

/*---------------------------------------------------------------------------*/

    poly_int8_to_int16(f, tmp1); //tmp1 = f

    for(i = 0; i < DIM_N; i++)
    {
        Tmp0[i] = 3158273 - abs(C2[i]); //Tmp0 = q2Q/2 - abs(C2)
    }

    idx0 = 0, idx1 = 0;
    norm1 = 0xFFFFFFFF; norm2 = 0xFFFFFFFF;
    for(i = 0; i < DIM_N; i++)
    {
        mask1 = (uint32_t)Tmp0[i] < norm1;
        mask2 = ((uint32_t)Tmp0[i] < norm2) & !mask1;

        norm2 = (norm1 & (-mask1)) | (norm2 & (~(-mask1)));
        idx1 = (idx0 & (-mask1)) | (idx1 & (~(-mask1)));

        norm1 = (Tmp0[i] & (-mask1)) | (norm1 & (~(-mask1)));
        idx0 = (i & (-mask1)) | (idx0 & (~(-mask1)));

        norm2 = (Tmp0[i] & (-mask2)) | (norm2 & (~(-mask2)));
        idx1 = (i & (-mask2)) | (idx1 & (~(-mask2)));
    }
    
    //tmp2 = f*x^idx0
    for(i = 0; i < DIM_N; i++)
    {
        mask1 = (i < idx0);
        tmp2[i] = (-tmp1[(DIM_N - abs(idx0 - i)) & ((1<<LOGN)-1)] & (-mask1)) | (tmp1[abs(i - idx0)] & (~(-mask1)));
    }
    
    //tmp3 = f*x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        mask1 = (i < idx1);
        tmp3[i] = (-tmp1[(DIM_N - abs(idx1 - i)) & ((1<<LOGN)-1)] & (-mask1)) | (tmp1[abs(i - idx1)] & (~(-mask1)));
    }

    //0

    //x^idx0
    for(i = 0; i < DIM_N; i++)
    {
        tmp4[i] = -sp[i] + tmp2[i];
    }
    mask1 = poly_check_distribution_binary(tmp4);
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp4[i] & (-mask1)) | (sp[i] & (~(-mask1)));
    }

    //-x^idx0
    for(i = 0; i < DIM_N; i++)
    {
        tmp4[i] = -sp[i] - tmp2[i];
    }
    mask1 = poly_check_distribution_binary(tmp4);
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp4[i] & (-mask1)) | (sp[i] & (~(-mask1)));
    }

    //x^idx0 + x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        tmp4[i] = -sp[i] + (tmp2[i] + tmp3[i]);
    }
    mask1 = poly_check_distribution_binary(tmp4);
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp4[i] & (-mask1)) | (sp[i] & (~(-mask1)));
    }

    //x^idx0 - x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        tmp4[i] = -sp[i] + (tmp2[i] - tmp3[i]);
    }
    mask1 = poly_check_distribution_binary(tmp4);
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp4[i] & (-mask1)) | (sp[i] & (~(-mask1)));
    }

    //-x^idx0 - x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        tmp4[i] = -sp[i] + (-tmp2[i] - tmp3[i]);
    }
    mask1 = poly_check_distribution_binary(tmp4);
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp4[i] & (-mask1)) | (sp[i] & (~(-mask1)));
    }

    //-x^idx0 + x^idx1
    for(i = 0; i < DIM_N; i++)
    {
        tmp4[i] = -sp[i] + (-tmp2[i] + tmp3[i]);
    }
    mask1 = poly_check_distribution_binary(tmp4);
    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (tmp4[i] & (-mask1)) | (sp[i] & (~(-mask1)));
    }

    // for(i = 7; i >= 0; i--)
    // {
    //     for(tmp = 0; tmp < MESSBYTES; tmp++)
    //     {
    //         r[tmp] = (r[tmp] << 1);
    //         r[tmp] += sp[512 + i * MESSBYTES + tmp];
    //     }
    // }

    for(tmp = 32; tmp >= 0; tmp -= 32)
    {
        vtmp3 = _mm256_setzero_si256();
        for(i = 7; i >= 0; i--)
        {
            vtmp1 = _mm256_load_si256((__m256i *)(sp + 512 + tmp * 8 + i * 32));
            vtmp1 = _mm256_packs_epi16(vtmp1, vtmp1);
            vtmp1 =_mm256_permute4x64_epi64(vtmp1, 0xd8);
            vtmp0 = _mm256_castsi256_si128(vtmp1);
            vtmp2 = _mm256_load_si256((__m256i *)(sp + 512 + tmp * 8 + i * 32 + 16));
            vtmp2 = _mm256_packs_epi16(vtmp2, vtmp2);
            vtmp2 =_mm256_permute4x64_epi64(vtmp2, 0xd8);
            vtmp8 = _mm256_castsi256_si128(vtmp2);
            vtmp1 = _mm256_insertf128_si256(vtmp1, vtmp0, 0);
            vtmp1 = _mm256_insertf128_si256(vtmp1, vtmp8, 1);
            vtmp3 = _mm256_slli_epi16(vtmp3, 1);
            vtmp3 = _mm256_add_epi16(vtmp3, vtmp1);
        }
        _mm256_store_si256((__m256i *)(r + tmp), vtmp3);
    }
    
    for(tmp = 32; tmp >= 0; tmp -= 32)
    {
        vtmp3 = _mm256_setzero_si256();
        for(i = 7; i >= 0; i--)
        {
            vtmp1 = _mm256_load_si256((__m256i *)(sp + tmp * 8 + i * 32));
            vtmp1 = _mm256_packs_epi16(vtmp1, vtmp1);
            vtmp1 =_mm256_permute4x64_epi64(vtmp1, 0xd8);
            vtmp0 = _mm256_castsi256_si128(vtmp1);
            vtmp2 = _mm256_load_si256((__m256i *)(sp + tmp * 8 + i * 32 + 16));
            vtmp2 = _mm256_packs_epi16(vtmp2, vtmp2);
            vtmp2 =_mm256_permute4x64_epi64(vtmp2, 0xd8);
            vtmp8 = _mm256_castsi256_si128(vtmp2);
            vtmp1 = _mm256_insertf128_si256(vtmp1, vtmp0, 0);
            vtmp1 = _mm256_insertf128_si256(vtmp1, vtmp8, 1);
            vtmp3 = _mm256_slli_epi16(vtmp3, 1);
            vtmp3 = _mm256_add_epi16(vtmp3, vtmp1);
        }
        _mm256_store_si256((__m256i *)(m + tmp), vtmp3);
    }

    sha3_512(r_hash, r, MESSBYTES);

    for(i = 0; i < MESSBYTES; i += 32)
    {
        vtmp1 = _mm256_load_si256((__m256i *)(m + i));
        vtmp2 = _mm256_load_si256((__m256i *)(r_hash + i));
        vtmp1 = _mm256_xor_si256(vtmp1, vtmp2);
        _mm256_store_si256((__m256i *)(m + i), vtmp1);
    }

}