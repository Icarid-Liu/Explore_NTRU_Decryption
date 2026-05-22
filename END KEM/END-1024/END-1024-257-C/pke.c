#include <string.h>
#include <stdio.h>
#include "param.h"
#include "mq_ntt.h"
#include "poly.h"
#include "fips202.h"
#include "code.h"

static  int16_t barrett_reduce_514(int32_t a) 
{
  int32_t t;
  t = ((int32_t)2040 * a);
  t = t >> 20;
  t *= 514;

  return a - t;
}

static  int32_t barrett_reduce_6316546(int64_t a) 
{
    int64_t t;
    t = ((int64_t)11140383 * a);
    t = t >> 46; 
    t *= 6316546;
    return a - t;
}

static  int16_t barrett_reduce_769(int32_t a) 
{
    int64_t t;
    t = ((int64_t)349070 * a);
    t = t >> 28;  
    t *= 769;
    return (int16_t)(a - t);
}

static int16_t barrett_reduce_257(int32_t a) 
{
    int64_t t;
    t = ((int64_t)16711935 * a);
    t = t >> 32;  
    t *= 257;
    return (int16_t)(a - t);
}

static int32_t barrett_reduce_197633(int32_t a) 
{
    int64_t t;
    t = ((int64_t)679 * a);
    t = t >> 27; 
    t *= 197633;
    return a - t;
}

static  int16_t div_6316546(int32_t a) 
{
    int64_t t;
    t = ((int64_t)11140383 * a);
    t = t >> 46; 
    t = (int16_t)t;
    t += (t>>15) & 1;
    return (int16_t)t;
}


static  int16_t montgomery_reduce_257(int32_t a)
{
  int32_t t;
  int16_t u;

  u=a*(-255);
  t=(int32_t)u*257;
  t=a-t;
  t>>=16;

  return t;
}

static  int16_t fqmul_257(int16_t a, int16_t b) 
{
  return montgomery_reduce_257((int32_t)a*b);
}

static int16_t montgomery_reduce_769(int32_t a)
{
    int32_t t;
    int16_t u, mask;

    u=a*(-767);
    t=(int32_t)u*769;
    t=a-t;
    t>>=16;

    return t;
}

static  int16_t fqmul_769(int16_t a, int16_t b) 
{
  return montgomery_reduce_769((int32_t)a*b);
}

static  int16_t montgomery_reduce_12289(int32_t a)
{
  int32_t t;
  int16_t u;

  u=a*(53249);
  t=(int32_t)u*12289;
  t=a-t;
  t>>=16;

  return t;
}

static  int16_t fqmul_12289(int16_t a, int16_t b) 
{
  return montgomery_reduce_12289((int32_t)a*b);
}


int END_encrypt(int16_t *h, int8_t *m, int8_t *r, int16_t *c)
{
    int8_t  m_prime[MESSBYTES];
    int16_t s[DIM_N];
    int16_t t1[DIM_N],t2[DIM_N],res[DIM_N];
    int16_t i, j, max_value;
    sha3_512(m_prime, r, MESSBYTES);

    for(i = 0; i < MESSBYTES; i++)
    {
        m_prime[i] = m[i] ^ m_prime[i];
    }

    for(i = 0; i < 8; i++)
    {
        for(j = 0; j < MESSBYTES; j++)
        {
            s[512 + i * MESSBYTES + j] = (r[j] & 1);
            r[j] = (r[j] >> 1);
        }
    }

    for(i = 0; i < 8; i++)
    {
        for(j = 0; j < MESSBYTES; j++)
        {
            s[i * MESSBYTES + j] = (m_prime[j] & 1);
            m_prime[j] = (m_prime[j] >> 1);
        }
    }

    //check e
    mq_poly_ntt_257(s);
    mq_poly_mul_ntt_257(c, h, s);
    mq_poly_intt_257(c);
    poly_round(c);

    return 0;

}

void END_decrypt(int16_t *h, int16_t *c, int8_t *f, int8_t *g, int8_t *F, int8_t *G, int16_t *w, int8_t *m)
{
    int16_t i, tmp, u1, idx0, idx1, idx2;
    uint32_t mask1, mask2, mask3;
    uint32_t norm1, norm2, norm3, norm4;
    int32_t t, num1, num2;
    int32_t C2[DIM_N], Tmp0[DIM_N], Tmp1[DIM_N];
    int16_t tmp1[DIM_N],tmp2[DIM_N], tmp3[DIM_N], tmp4[DIM_N], tmp5[DIM_N], tmp6[DIM_N], tmp7[DIM_N];
    int16_t c1[DIM_N], c2[DIM_N], sp[DIM_N], ep[DIM_N], Fp[DIM_N], Gp[DIM_N];

    int8_t r[MESSBYTES], r_hash[MESSBYTES];

    poly_deround(c);

    for(i = 0; i < DIM_N; i++)
    {
        tmp1[i] = 1; //tmp1 = vector(1)
    }

    poly_int8_to_int16(f, tmp2); //tmp2 = f
    poly_int8_to_int16(g, tmp3); //tmp3 = g
    poly_int8_to_int16(F, tmp4); //tmp4 = F
    poly_int8_to_int16(G, tmp5); //tmp5 = G
    norm2 = poly_norm2(tmp4) + poly_norm2(tmp5);

    mq_poly_ntt_257(tmp1); //tmp1 = vector(1)_257
    mq_poly_ntt_257(c); //c = c_257

    mq_poly_ntt_257(tmp2); //tmp2 = f_257
    mq_poly_ntt_257(tmp3); //tmp3 = g_257
    mq_poly_mul_ntt_257(c1, c, tmp2); //c1 = cf_257

    mq_poly_mul_ntt_257(tmp2, tmp1, tmp3); //tmp2 = g.1_257

    mask1 = 0;
    for(i = 0; i < DIM_N; i++)
    {
        c1[i] = montgomery_reduce_257((2 * c1[i]) - tmp2[i]);
        mask1 ^= (g[i] & 1);
    }
    mq_poly_intt_257_n(c1); //c1 = cf-g.1_mod257

    for(i = 0; i < DIM_N; i++)
    {
        idx0 = barrett_reduce_514(c1[i] * 258 - 257);
        idx1 = barrett_reduce_514(c1[i] * 258);
        idx2 = (idx0 & (-mask1)) | (idx1 & (~(-mask1)));
        u1 = idx2 > 256;
        c1[i] = idx2 - (u1 * 514);
    }

    mq_poly_ntt_257(tmp4); //tmp4 = F_257
    mq_poly_ntt_257(tmp5); //tmp5 = G_257
    mq_poly_mul_ntt_257(c2, c, tmp4); //c2 = cF_257
    mq_poly_mul_ntt_257(tmp4, tmp1, tmp5); //tmp4 = G.1_257
    mask1 = 0;
    for(i = 0; i < DIM_N; i++)
    {
        c2[i] = montgomery_reduce_257((2 * c2[i]) - tmp4[i]);
        mask1 ^= (G[i] & 1);
    }
    mq_poly_intt_257_n(c2); //c2 = cF-G.1_mod257

    for(i = 0; i < DIM_N; i++)
    {
        idx0 = barrett_reduce_514(c2[i] * 258 - 257);
        idx1 = barrett_reduce_514(c2[i] * 258);
        idx2 = (idx0 & (-mask1)) | (idx1 & (~(-mask1)));
        u1 = idx2 > 256;
        c2[i] = idx2 - (u1 * 514);
    }

    for(i = 0; i < DIM_N; i++)
    {
        tmp2[i] = montgomery_reduce_257(w[i]); //tmp2 = w % 257
        tmp3[i] = montgomery_reduce_12289(w[i] * 4091); //tmp3 = w
    }

    memcpy(tmp1, c1, DIM_N * sizeof(int16_t)); // tmp1 = c1
    mq_poly_ntt_257(tmp1); //tmp1 = c_257
    mq_poly_ntt_257(tmp2); //tmp2 = w_257
    mq_poly_mul_ntt_257(tmp4, tmp1, tmp2); //tmp4 = c1w_257
    mq_poly_intt_257(tmp4); // tmp4 = c1w % 257

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
    mq_poly_ntt_769(tmp1); //tmp1 = f_769
    mq_poly_ntt_769(c2); //c2 = c2_769
    mq_poly_mul_ntt_769(tmp2, tmp1, c2); //tmp2 = c2f_769

    for(i = 0; i < DIM_N; i++)
    {
        tmp3[i] = fqmul_769(w[i], 171); //tmp3 = w % 769
    }
    mq_poly_ntt_769(tmp3); //tmp3 = w_769
    mq_poly_mul_ntt_769(tmp4, tmp1, tmp3); //tmp4 = wf_769

    poly_int8_to_int16(F, tmp3); //tmp3 = F
    mq_poly_ntt_769(tmp3); //tmp3 = F_769
    for(i = 0; i < DIM_N; i++)
    {
        tmp3[i] = fqmul_769(fqmul_769(tmp3[i], 511) - fqmul_769(tmp4[i], 19), 171); //tmp3 = QF-wf_769
    }

    mq_poly_ntt_769(c1); //c1 = c1_769
    mq_poly_mul_ntt_769(tmp5, c1, tmp3); //tmp5 = c1(QF-wf)_769

    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = fqmul_769(tmp5[i] - tmp2[i], 675); //sp = sp_769
    }
    mq_poly_intt_769_n(sp); //sp = sp

    for(i = 0; i < DIM_N; i++)
    {
        sp[i] = (sp[i] + 1) >> 1;
    }

    memcpy(tmp3, sp, DIM_N * sizeof(int16_t)); // tmp3 = sp
    mq_poly_ntt_257(tmp3); //tmp3 = sp_257
    mq_poly_mul_ntt_257(tmp2, h, tmp3); // tmp2 = hsp_257
    for(i = 0; i < DIM_N; i++)
    {
        ep[i] = montgomery_reduce_257(c[i] - tmp2[i]); //ep = c - hsp_257
    }
    mq_poly_intt_257_n(ep); //ep = ep

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
    mq_poly_ntt_257(tmp6); //tmp2 = C2-2qw^idx0/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1); //tmp4 = C2-qwe1/qQ * f_257
    mq_poly_intt_257_n(tmp7); //tmp7 = C2-2qw^idx0/2qQ * f
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
    mq_poly_ntt_257(tmp6); //tmp2 = C2-2qw^idx0/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1); //tmp4 = C2-qwe1/qQ * f_257
    mq_poly_intt_257_n(tmp7); //tmp7 = C2-2qw^idx0/2qQ * f
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
    mq_poly_ntt_257(tmp6); //tmp2 = C2-2q(w^idx0+w^idx1)/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1); //tmp4 = C2-2q(w^idx0+w^idx1)/2qQ * f_257
    mq_poly_intt_257_n(tmp7); //tmp7 = C2-2q(w^idx0+w^idx1)/2qQ * f
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
    mq_poly_ntt_257(tmp6); //tmp2 = C2-2q(w^idx0+w^idx1)/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1); //tmp4 = C2-2q(w^idx0+w^idx1)/2qQ * f_257
    mq_poly_intt_257_n(tmp7); //tmp7 = C2-2q(w^idx0+w^idx1)/2qQ * f
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
    mq_poly_ntt_257(tmp6); //tmp2 = C2-2q(w^idx0+w^idx1)/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1); //tmp4 = C2-2q(w^idx0+w^idx1)/2qQ * f_257
    mq_poly_intt_257_n(tmp7); //tmp7 = C2-2q(w^idx0+w^idx1)/2qQ * f
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
    mq_poly_ntt_257(tmp6); //tmp2 = C2-2q(w^idx0+w^idx1)/2qQ_257
    mq_poly_mul_ntt_257(tmp7, tmp6, tmp1); //tmp4 = C2-2q(w^idx0+w^idx1)/2qQ * f_257
    mq_poly_intt_257_n(tmp7); //tmp7 = C2-2q(w^idx0+w^idx1)/2qQ * f
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

    for(i = 7; i >= 0; i--)
    {
        for(tmp = 0; tmp < MESSBYTES; tmp++)
        {
            r[tmp] = (r[tmp] << 1);
            r[tmp] += sp[512 + i * MESSBYTES + tmp];
        }
    }

    for(i = 7; i >= 0; i--)
    {
        for(tmp = 0; tmp < MESSBYTES; tmp++)
        {
            m[tmp] = (m[tmp] << 1);
            m[tmp] += sp[i * MESSBYTES + tmp];
            
        }
    }

    sha3_512(r_hash, r, MESSBYTES);

    for(i = 0; i < MESSBYTES; i++)
    {
        m[i] = m[i] ^ r_hash[i];
    }

}