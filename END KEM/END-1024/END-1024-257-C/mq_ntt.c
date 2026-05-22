#include <stdint.h>
#include <stdio.h>
#include "param.h"
#include "mq_ntt.h"
#include "mq_ntt_param.h"

static int16_t montgomery_reduce_257(int32_t a)
{
  int32_t t;
  int16_t u;

  u=a*(-255);
  t=(int32_t)u*257;
  t=a-t;
  t>>=16;
  return t;
}

static int16_t fqmul_257(int16_t a, int16_t b) 
{
  return montgomery_reduce_257((int32_t)a*b);
}


void mq_poly_ntt_257(int16_t *a) 
{
  unsigned int len,start,j,k;
  int16_t t,zeta;

  k=1;
  for(len=DIM_N>>1;len>=8;len>>=1) 
  {
    for(start=0;start<DIM_N;start=j+len) 
    {
      zeta=f_257[k++];
      for(j=start;j<start+len;j++) 
      {
        t=fqmul_257(zeta,a[j+len]);
        a[j+len]=a[j]-t;
        a[j]=a[j]+t;
      }
    }
  }
}

void mq_poly_intt_257(int16_t *a) 
{
  unsigned int start,len,j,k;
  int16_t t,zeta;

  k=0;
  for(len=8;len<=DIM_N>>1;len<<=1)
  {
    for(start=0;start<DIM_N;start=j+len)
    {
      zeta=fn_257[k++];
      for(j=start;j<start+len;j++) 
      {
        t=a[j];
        a[j]=(t+a[j+len]);
        a[j+len]=t-a[j+len];
        a[j+len]=fqmul_257(zeta,a[j+len]);
      }
    }
  }

  for(j=0;j<DIM_N;j++)
  {
    a[j]=fqmul_257(a[j],fn_257[127]);
    a[j]+=(a[j]>>15)&257;
  }
}

void mq_poly_intt_257_n(int16_t *a) 
{
  unsigned int start,len,j,k;
  int16_t t,zeta;

  k=0;
  for(len=8;len<=DIM_N>>1;len<<=1)
  {
    for(start=0;start<DIM_N;start=j+len)
    {
      zeta=fn_257[k++];
      for(j=start;j<start+len;j++) 
      {
        t=a[j];
        a[j]=(t+a[j+len]);
        a[j+len]=t-a[j+len];
        a[j+len]=fqmul_257(zeta,a[j+len]);
      }
    }
  }

  for(j=0;j<DIM_N;j++)
  {
    a[j]=fqmul_257(a[j],fn_257[127]);
    // a[j]+=(a[j]>>15)&257;
  }
}

void basemul_257(int16_t *r,int16_t *a,int16_t *b,int16_t zeta)
{
  int16_t k;
  r[0]=fqmul_257(a[0],b[0]);
  k=fqmul_257(a[7],b[1]);
  k+=fqmul_257(a[6],b[2]);
  k+=fqmul_257(a[5],b[3]);
  k+=fqmul_257(a[4],b[4]);
  k+=fqmul_257(a[3],b[5]);
  k+=fqmul_257(a[2],b[6]);
  k+=fqmul_257(a[1],b[7]);
  r[0]+=fqmul_257(k,zeta);

  r[1]=fqmul_257(a[1],b[0]);
  r[1]+=fqmul_257(a[0],b[1]);
  k=fqmul_257(a[7],b[2]);
  k+=fqmul_257(a[6],b[3]);
  k+=fqmul_257(a[5],b[4]);
  k+=fqmul_257(a[4],b[5]);
  k+=fqmul_257(a[3],b[6]);
  k+=fqmul_257(a[2],b[7]);
  r[1]+=fqmul_257(k,zeta);

  r[2]=fqmul_257(a[2],b[0]);
  r[2]+=fqmul_257(a[1],b[1]);
  r[2]+=fqmul_257(a[0],b[2]);
  k=fqmul_257(a[7],b[3]);
  k+=fqmul_257(a[6],b[4]);
  k+=fqmul_257(a[5],b[5]);
  k+=fqmul_257(a[4],b[6]);
  k+=fqmul_257(a[3],b[7]);
  r[2]+=fqmul_257(k,zeta);

  r[3]=fqmul_257(a[3],b[0]);
  r[3]+=fqmul_257(a[2],b[1]);
  r[3]+=fqmul_257(a[1],b[2]);
  r[3]+=fqmul_257(a[0],b[3]);
  k=fqmul_257(a[7],b[4]);
  k+=fqmul_257(a[6],b[5]);
  k+=fqmul_257(a[5],b[6]);
  k+=fqmul_257(a[4],b[7]);
  r[3]+=fqmul_257(k,zeta);

  r[4]=fqmul_257(a[4],b[0]);
  r[4]+=fqmul_257(a[3],b[1]);
  r[4]+=fqmul_257(a[2],b[2]);
  r[4]+=fqmul_257(a[1],b[3]);
  r[4]+=fqmul_257(a[0],b[4]);
  k=fqmul_257(a[7],b[5]);
  k+=fqmul_257(a[6],b[6]);
  k+=fqmul_257(a[5],b[7]);
  r[4]+=fqmul_257(k,zeta);

  r[5]=fqmul_257(a[5],b[0]);
  r[5]+=fqmul_257(a[4],b[1]);
  r[5]+=fqmul_257(a[3],b[2]);
  r[5]+=fqmul_257(a[2],b[3]);
  r[5]+=fqmul_257(a[1],b[4]);
  r[5]+=fqmul_257(a[0],b[5]);
  k=fqmul_257(a[7],b[6]);
  k+=fqmul_257(a[6],b[7]);
  r[5]+=fqmul_257(k,zeta);

  k=fqmul_257(a[7],b[7]);
  r[6]=fqmul_257(k,zeta);
  r[6]+=fqmul_257(a[6],b[0]);
  r[6]+=fqmul_257(a[5],b[1]);
  r[6]+=fqmul_257(a[4],b[2]);
  r[6]+=fqmul_257(a[3],b[3]);
  r[6]+=fqmul_257(a[2],b[4]);
  r[6]+=fqmul_257(a[1],b[5]);
  r[6]+=fqmul_257(a[0],b[6]);

  r[7]=fqmul_257(a[7],b[0]);
  r[7]+=fqmul_257(a[6],b[1]);
  r[7]+=fqmul_257(a[5],b[2]);
  r[7]+=fqmul_257(a[4],b[3]);
  r[7]+=fqmul_257(a[3],b[4]);
  r[7]+=fqmul_257(a[2],b[5]);
  r[7]+=fqmul_257(a[1],b[6]);
  r[7]+=fqmul_257(a[0],b[7]);
}

void mq_poly_mul_ntt_257(int16_t *r, int16_t *a, int16_t *b)
{
    int i;
    for(i=0;i<DIM_N/16;i++) 
    {
        basemul_257(r+16*i, a+16*i, b+16*i, f_257[64+i]);
        basemul_257(r+16*i+8, a+16*i+8, b+16*i+8, -f_257[64+i]);
    }
}

void basemul_257_mq(int16_t *r, int16_t *a, int16_t *b, int16_t zeta)
{
  int16_t k;
  r[0] = fqmul_257(a[0], b[0]);
  k = fqmul_257(a[7], b[1]);
  k += fqmul_257(a[6], b[2]);
  k += fqmul_257(a[5], b[3]);
  k += fqmul_257(a[4], b[4]);
  k += fqmul_257(a[3], b[5]);
  k += fqmul_257(a[2], b[6]);
  k += fqmul_257(a[1], b[7]);
  r[0] += fqmul_257(k, zeta);
  r[0] = fqmul_257(r[0], 1);
  r[0] += (r[0] >> 15) & Q;

  r[1] = fqmul_257(a[1], b[0]);
  r[1] += fqmul_257(a[0], b[1]);
  k = fqmul_257(a[7], b[2]);
  k += fqmul_257(a[6], b[3]);
  k += fqmul_257(a[5], b[4]);
  k += fqmul_257(a[4], b[5]);
  k += fqmul_257(a[3], b[6]);
  k += fqmul_257(a[2], b[7]);
  r[1] += fqmul_257(k, zeta);
  r[1] = fqmul_257(r[1], 1);
  r[1] += (r[1] >> 15) & Q;

  r[2] = fqmul_257(a[2], b[0]);
  r[2] += fqmul_257(a[1], b[1]);
  r[2] += fqmul_257(a[0], b[2]);
  k = fqmul_257(a[7], b[3]);
  k += fqmul_257(a[6], b[4]);
  k += fqmul_257(a[5], b[5]);
  k += fqmul_257(a[4], b[6]);
  k += fqmul_257(a[3], b[7]);
  r[2] += fqmul_257(k, zeta);
  r[2] = fqmul_257(r[2], 1);
  r[2] += (r[2] >> 15) & Q;

  r[3] = fqmul_257(a[3], b[0]);
  r[3] += fqmul_257(a[2], b[1]);
  r[3] += fqmul_257(a[1], b[2]);
  r[3] += fqmul_257(a[0], b[3]);
  k = fqmul_257(a[7], b[4]);
  k += fqmul_257(a[6], b[5]);
  k += fqmul_257(a[5], b[6]);
  k += fqmul_257(a[4], b[7]);
  r[3] += fqmul_257(k, zeta);
  r[3] = fqmul_257(r[3], 1);
  r[3] += (r[3] >> 15) & Q;

  r[4] = fqmul_257(a[4], b[0]);
  r[4] += fqmul_257(a[3], b[1]);
  r[4] += fqmul_257(a[2], b[2]);
  r[4] += fqmul_257(a[1], b[3]);
  r[4] += fqmul_257(a[0], b[4]);
  k = fqmul_257(a[7], b[5]);
  k += fqmul_257(a[6], b[6]);
  k += fqmul_257(a[5], b[7]);
  r[4] += fqmul_257(k, zeta);
  r[4] = fqmul_257(r[4], 1);
  r[4] += (r[4] >> 15) & Q;

  r[5] = fqmul_257(a[5], b[0]);
  r[5] += fqmul_257(a[4], b[1]);
  r[5] += fqmul_257(a[3], b[2]);
  r[5] += fqmul_257(a[2], b[3]);
  r[5] += fqmul_257(a[1], b[4]);
  r[5] += fqmul_257(a[0], b[5]);
  k = fqmul_257(a[7], b[6]);
  k += fqmul_257(a[6], b[7]);
  r[5] += fqmul_257(k, zeta);
  r[5] = fqmul_257(r[5], 1);
  r[5] += (r[5] >> 15) & Q;

  k = fqmul_257(a[7], b[7]);
  r[6]=fqmul_257(k, zeta);
  r[6] += fqmul_257(a[6], b[0]);
  r[6] += fqmul_257(a[5], b[1]);
  r[6] += fqmul_257(a[4], b[2]);
  r[6] += fqmul_257(a[3], b[3]);
  r[6] += fqmul_257(a[2], b[4]);
  r[6] += fqmul_257(a[1], b[5]);
  r[6] += fqmul_257(a[0], b[6]);
  r[6] = fqmul_257(r[6], 1);
  r[6] += (r[6] >> 15) & Q;

  r[7] = fqmul_257(a[7], b[0]);
  r[7] += fqmul_257(a[6], b[1]);
  r[7] += fqmul_257(a[5], b[2]);
  r[7] += fqmul_257(a[4], b[3]);
  r[7] += fqmul_257(a[3], b[4]);
  r[7] += fqmul_257(a[2], b[5]);
  r[7] += fqmul_257(a[1], b[6]);
  r[7] += fqmul_257(a[0], b[7]);
  r[7] = fqmul_257(r[7], 1);
  r[7] += (r[7] >> 15) & Q;
}

void mq_poly_mul_ntt_257_mq(int16_t *r,  int16_t *a,  int16_t *b)
{
  int i;
  for(i = 0; i < DIM_N / 16; i++) 
  {
    basemul_257_mq(r + 16 * i, a + 16 * i, b + 16 * i, f_257[64 + i]);
    basemul_257_mq(r + 16 * i + 8, a + 16 * i + 8, b + 16 * i + 8, -f_257[64 + i]);
  }
}

static void base_adjoint_257(int16_t *r,int16_t *a, int16_t zeta)
{
  r[0] = a[0];
  r[1] = fqmul_257(a[7], zeta);
  r[2] = fqmul_257(a[6], zeta);
  r[3] = fqmul_257(a[5], zeta);
  r[4] = fqmul_257(a[4], zeta);
  r[5] = fqmul_257(a[3], zeta);
  r[6] = fqmul_257(a[2], zeta);
  r[7] = fqmul_257(a[1], zeta);
}

void mq_poly_adjoint_ntt_257(int16_t *r, int16_t *a)
{
    int i;
    for(i=0;i<DIM_N/16;i++) 
    {
      base_adjoint_257(r+DIM_N-16*(i+1)+8, a+16*i, f_257[64+i]);
      base_adjoint_257(r+DIM_N-16*(i+1), a+16*i+8, -f_257[64+i]);
    }
}

int baseinv_257(int16_t *r, int16_t *a, int16_t zeta)
{
  int16_t t, k, flag;
  int16_t b[4], c[2], e, f[4]; 

  b[0] = fqmul_257(a[1], a[7]);
  t = fqmul_257(a[3], a[5]);
  b[0] += t;
  t = fqmul_257(a[2], a[6]);
  b[0] -= t;
  b[0] = fqmul_257(b[0], 2); //2
  t = fqmul_257(a[4], a[4]);
  b[0] -= t;
  b[0] = fqmul_257(b[0], zeta);
  t = fqmul_257(a[0], a[0]);
  b[0] += t;

  b[1] = fqmul_257(a[3], a[7]);
  t = fqmul_257(a[4], a[6]);
  b[1] -= t;
  b[1] = fqmul_257(b[1], 2); //2
  t = fqmul_257(a[5], a[5]);
  b[1] += t;
  b[1] = fqmul_257(b[1], zeta);
  t = fqmul_257(a[0], a[2]);
  t = fqmul_257(t, 2); //2
  b[1] += t;
  t = fqmul_257(a[1], a[1]);
  b[1] -= t;

  b[2] = fqmul_257(a[5], a[7]);
  b[2] = fqmul_257(b[2], 2); //2
  t = fqmul_257(a[6], a[6]);
  b[2] -= t;
  b[2] = fqmul_257(b[2], zeta);
  t = fqmul_257(a[0], a[4]);
  k = fqmul_257(a[1], a[3]);
  t -= k;
  t = fqmul_257(t, 2); //2
  b[2] += t;
  t = fqmul_257(a[2], a[2]);
  b[2] += t;

  b[3] = fqmul_257(a[7], a[7]);
  b[3] = fqmul_257(b[3], zeta);
  t = fqmul_257(a[0], a[6]);
  k = fqmul_257(a[2], a[4]);
  t += k;
  k = fqmul_257(a[1], a[5]);
  t -= k;
  t = fqmul_257(t, 2); //2
  b[3] += t;
  t = fqmul_257(a[3], a[3]);
  b[3] -= t;

  c[0] = fqmul_257(b[1], b[3]);
  c[0] = fqmul_257(c[0], 2); //2
  t = fqmul_257(b[2], b[2]);
  c[0] -= t;
  c[0] = fqmul_257(c[0], zeta);
  t = fqmul_257(b[0], b[0]);
  c[0] += t;

  c[1] = fqmul_257(b[3], b[3]);
  c[1] = fqmul_257(c[1], zeta);
  t = fqmul_257(b[0], b[2]);
  t = fqmul_257(t, 2); //2
  c[1] += t;
  t = fqmul_257(b[1], b[1]);
  c[1] -= t;

  e = fqmul_257(c[1], c[1]);
  e = fqmul_257(e, zeta);
  t = fqmul_257(c[0], c[0]);
  e += t;
  e += (e >> 15) & Q;
  flag = e == 0;
  e = qinv[e];

  c[0] = fqmul_257(e, c[0]);

  c[1] = fqmul_257(e, c[1]);
  c[1] = fqmul_257(c[1], -1);

  f[0] = fqmul_257(c[1], b[2]);
  f[0] = fqmul_257(f[0], zeta);
  t = fqmul_257(c[0], b[0]);
  f[0] = t- f[0];

  f[1] = fqmul_257(c[1], b[3]);
  f[1] = fqmul_257(f[1], zeta);
  t = fqmul_257(c[0], b[1]);
  f[1] -= t;

  f[2] = fqmul_257(c[0], b[2]);
  t = fqmul_257(c[1], b[0]);
  f[2] += t;

  f[3] = fqmul_257(c[0], b[3]);
  t = fqmul_257(c[1], b[1]);
  f[3] += t;
  f[3] = fqmul_257(f[3], -1);

  r[0] = fqmul_257(f[1], a[6]);
  t = fqmul_257(f[2], a[4]);
  r[0] += t;
  t = fqmul_257(f[3], a[2]);
  r[0] += t;
  r[0] = fqmul_257(r[0], zeta);
  t = fqmul_257(f[0], a[0]);
  r[0] = t - r[0];

  r[1] = fqmul_257(f[1], a[7]);
  t = fqmul_257(f[2], a[5]);
  r[1] += t;
  t  = fqmul_257(f[3], a[3]);
  r[1] += t;
  r[1] = fqmul_257(r[1], zeta);
  t = fqmul_257(f[0], a[1]);
  r[1] -= t;

  r[2] = fqmul_257(f[2], a[6]);
  t = fqmul_257(f[3], a[4]);
  r[2] += t;
  r[2] = fqmul_257(r[2], zeta);
  t = fqmul_257(f[0], a[2]);
  k = fqmul_257(f[1], a[0]);
  t += k;
  r[2] = t - r[2];

  r[3] = fqmul_257(f[2], a[7]);
  t = fqmul_257(f[3], a[5]);
  r[3] += t;
  r[3] = fqmul_257(r[3], zeta);
  t = fqmul_257(f[0], a[3]);
  k = fqmul_257(f[1], a[1]);
  t += k;
  r[3] -= t;

  r[4] = fqmul_257(f[3], a[6]);
  r[4] = fqmul_257(r[4], zeta);
  t = fqmul_257(f[0], a[4]);
  k = fqmul_257(f[1], a[2]);
  t += k;
  k = fqmul_257(f[2], a[0]);
  t += k;
  r[4] = t - r[4];

  r[5] = fqmul_257(f[3], a[7]);
  r[5] = fqmul_257(r[5], zeta);
  t = fqmul_257(f[0], a[5]);
  k = fqmul_257(f[1], a[3]);
  t += k;
  k = fqmul_257(f[2], a[1]);
  t += k;
  r[5] -= t;

  r[6] = fqmul_257(f[1], a[4]);
  t = fqmul_257(f[2], a[2]);
  r[6] += t;
  t = fqmul_257(f[3], a[0]);
  r[6] += t;
  t = fqmul_257(f[0], a[6]);
  r[6] += t;

  r[7] = fqmul_257(f[0], a[7]);
  t = fqmul_257(f[1], a[5]);
  r[7] += t;
  t = fqmul_257(f[2], a[3]);
  r[7] += t;
  t = fqmul_257(f[3], a[1]);
  r[7] += t;
  r[7] = fqmul_257(r[7], -1);

  return flag;
}

//if f is invertible return 0
int mq_poly_inv_ntt_257(int16_t *r, int16_t *a)
{
  int i, check = 0;
  for(i = 0; i < DIM_N / 16; i++) 
  {
    check += baseinv_257(r + 16 * i, a + 16 * i, -f_257[64 + i]);
    check += baseinv_257(r + 16 * i + 8, a + 16 * i + 8, f_257[64 + i]);
  }
  return check;
}

static int16_t montgomery_reduce_769(int32_t a)
{
  int32_t t;
  int16_t u;

  u = a * (-767);
  t = (int32_t)u * 769;
  t = a - t;
  t >>= 16;
  return t;
}

static int16_t fqmul_769(int16_t a, int16_t b) 
{
  return montgomery_reduce_769((int32_t)a * b);
}

void mq_poly_ntt_769(int16_t *a) 
{
  unsigned int len, start, j, k;
  int16_t t, zeta;

  k = 1;
  for(len = DIM_N >> 1; len >= 8; len >>= 1) 
  {
    for(start = 0; start < DIM_N; start = j + len) 
    {
      zeta = f_769[k++];
      for(j = start; j < start + len; j++) 
      {
        t = fqmul_769(zeta, a[j + len]);
        a[j + len] = a[j] - t;
        a[j] = a[j] + t;
      }
    }
  }
}


void mq_poly_intt_769(int16_t *a) 
{
  unsigned int start, len, j, k;
  int16_t t, zeta;

  k = 0;
  for(len = 8; len <= DIM_N >> 1; len <<= 1)
  {
    for(start = 0; start < DIM_N; start = j + len)
    {
      zeta = fn_769[k++];
      for(j = start; j < start + len; j++) 
      {
        t = a[j];
        a[j] = (t + a[j + len]);
        a[j+len] = t - a[j + len];
        a[j+len] = fqmul_769(zeta, a[j + len]);
      }
    }
  }

  for(j = 0; j < DIM_N; j++)
  {
    a[j] = fqmul_769(a[j], fn_769[127]);
    a[j] += (a[j] >> 15) & 769;
  }
}

void mq_poly_intt_769_n(int16_t *a) 
{
  unsigned int start,len,j,k;
  int16_t t,zeta;

  k=0;
  for(len=8;len<=DIM_N>>1;len<<=1)
  {
    for(start=0;start<DIM_N;start=j+len)
    {
      zeta=fn_769[k++];
      for(j=start;j<start+len;j++) 
      {
        t=a[j];
        a[j]=(t+a[j+len]);
        a[j+len]=t-a[j+len];
        a[j+len]=fqmul_769(zeta,a[j+len]);
      }
    }
  }

  for(j=0;j<DIM_N;j++)
  {
    a[j]=fqmul_769(a[j],fn_769[127]);
    // a[j]+=(a[j]>>15)&769;
  }
}

void basemul_769(int16_t *r,int16_t *a,int16_t *b,int16_t zeta)
{
  int16_t k;
  r[0]=fqmul_769(a[0],b[0]);
  k=fqmul_769(a[7],b[1]);
  k+=fqmul_769(a[6],b[2]);
  k+=fqmul_769(a[5],b[3]);
  k+=fqmul_769(a[4],b[4]);
  k+=fqmul_769(a[3],b[5]);
  k+=fqmul_769(a[2],b[6]);
  k+=fqmul_769(a[1],b[7]);
  r[0]+=fqmul_769(k,zeta);

  r[1]=fqmul_769(a[1],b[0]);
  r[1]+=fqmul_769(a[0],b[1]);
  k=fqmul_769(a[7],b[2]);
  k+=fqmul_769(a[6],b[3]);
  k+=fqmul_769(a[5],b[4]);
  k+=fqmul_769(a[4],b[5]);
  k+=fqmul_769(a[3],b[6]);
  k+=fqmul_769(a[2],b[7]);
  r[1]+=fqmul_769(k,zeta);

  r[2]=fqmul_769(a[2],b[0]);
  r[2]+=fqmul_769(a[1],b[1]);
  r[2]+=fqmul_769(a[0],b[2]);
  k=fqmul_769(a[7],b[3]);
  k+=fqmul_769(a[6],b[4]);
  k+=fqmul_769(a[5],b[5]);
  k+=fqmul_769(a[4],b[6]);
  k+=fqmul_769(a[3],b[7]);
  r[2]+=fqmul_769(k,zeta);

  r[3]=fqmul_769(a[3],b[0]);
  r[3]+=fqmul_769(a[2],b[1]);
  r[3]+=fqmul_769(a[1],b[2]);
  r[3]+=fqmul_769(a[0],b[3]);
  k=fqmul_769(a[7],b[4]);
  k+=fqmul_769(a[6],b[5]);
  k+=fqmul_769(a[5],b[6]);
  k+=fqmul_769(a[4],b[7]);
  r[3]+=fqmul_769(k,zeta);

  r[4]=fqmul_769(a[4],b[0]);
  r[4]+=fqmul_769(a[3],b[1]);
  r[4]+=fqmul_769(a[2],b[2]);
  r[4]+=fqmul_769(a[1],b[3]);
  r[4]+=fqmul_769(a[0],b[4]);
  k=fqmul_769(a[7],b[5]);
  k+=fqmul_769(a[6],b[6]);
  k+=fqmul_769(a[5],b[7]);
  r[4]+=fqmul_769(k,zeta);

  r[5]=fqmul_769(a[5],b[0]);
  r[5]+=fqmul_769(a[4],b[1]);
  r[5]+=fqmul_769(a[3],b[2]);
  r[5]+=fqmul_769(a[2],b[3]);
  r[5]+=fqmul_769(a[1],b[4]);
  r[5]+=fqmul_769(a[0],b[5]);
  k=fqmul_769(a[7],b[6]);
  k+=fqmul_769(a[6],b[7]);
  r[5]+=fqmul_769(k,zeta);

  k=fqmul_769(a[7],b[7]);
  r[6]=fqmul_769(k,zeta);
  r[6]+=fqmul_769(a[6],b[0]);
  r[6]+=fqmul_769(a[5],b[1]);
  r[6]+=fqmul_769(a[4],b[2]);
  r[6]+=fqmul_769(a[3],b[3]);
  r[6]+=fqmul_769(a[2],b[4]);
  r[6]+=fqmul_769(a[1],b[5]);
  r[6]+=fqmul_769(a[0],b[6]);

  r[7]=fqmul_769(a[7],b[0]);
  r[7]+=fqmul_769(a[6],b[1]);
  r[7]+=fqmul_769(a[5],b[2]);
  r[7]+=fqmul_769(a[4],b[3]);
  r[7]+=fqmul_769(a[3],b[4]);
  r[7]+=fqmul_769(a[2],b[5]);
  r[7]+=fqmul_769(a[1],b[6]);
  r[7]+=fqmul_769(a[0],b[7]);

}

void mq_poly_mul_ntt_769(int16_t *r, int16_t *a, int16_t *b)
{
  int i;
  for(i=0;i<DIM_N/16;i++) 
  {
    basemul_769(r+16*i, a+16*i, b+16*i, f_769[64+i]);
    basemul_769(r+16*i+8, a+16*i+8, b+16*i+8, -f_769[64+i]);
  }
}

static void base_adjoint_769(int16_t *r,int16_t *a, int16_t zeta)
{
  r[0] = a[0];
  r[1] = fqmul_769(a[7], zeta);
  r[2] = fqmul_769(a[6], zeta);
  r[3] = fqmul_769(a[5], zeta);
  r[4] = fqmul_769(a[4], zeta);
  r[5] = fqmul_769(a[3], zeta);
  r[6] = fqmul_769(a[2], zeta);
  r[7] = fqmul_769(a[1], zeta);
}

void mq_poly_adjoint_ntt_769(int16_t *r, int16_t *a)
{
    int i;
    for(i=0;i<DIM_N/16;i++) 
    {
      base_adjoint_769(r+DIM_N-16*(i+1)+8, a+16*i, f_769[64+i]);
      base_adjoint_769(r+DIM_N-16*(i+1), a+16*i+8, -f_769[64+i]);
    }
}

static int16_t montgomery_reduce_12289(int32_t a)
{
  int32_t t;
  int16_t u;

  u=a*(53249);
  t=(int32_t)u*12289;
  t=a-t;
  t>>=16;
  return t;
}

static int16_t fqmul_12289(int16_t a, int16_t b) 
{
  return montgomery_reduce_12289((int32_t)a*b);
}

void mq_poly_ntt_12289(int16_t *a) 
{
  unsigned int len,start,j,k;
  int16_t t,zeta;

  k=1;
  for(len=DIM_N>>1;len>=1;len>>=1) 
  {
    for(start=0;start<DIM_N;start=j+len) 
    {
      zeta=f_12289[k++];
      for(j=start;j<start+len;j++) 
      {
        t=fqmul_12289(zeta,a[j+len]);
        t+=(t>>15)&12289;
        a[j]+=(a[j]>>15)&12289;
        a[j+len]=a[j]-t;
        a[j]=a[j]+t;
        a[j]-=12289;
      }
    }
  }

}

void mq_poly_intt_12289(int16_t *a) 
{
  unsigned int start,len,j,k;
  int16_t t,zeta;

  k=0;
  for(len=1;len<=DIM_N>>1;len<<=1)
  {
    for(start=0;start<DIM_N;start=j+len)
    {
      zeta=fn_12289[k++];
      for(j=start;j<start+len;j++) 
      {
        t=a[j];
        a[j]=t+a[j+len];
        a[j]-=12289;
        a[j]+=(a[j]>>15)&12289;
        a[j+len]=t-a[j+len];
        a[j+len]=fqmul_12289(zeta,a[j+len]);
        a[j+len]+=(a[j+len]>>15)&12289;
      }
    }
  }

  for(j=0;j<DIM_N;j++)
  {
    a[j]=fqmul_12289(a[j],fn_12289[1023]);
    // a[j]+=(a[j]>>15)&12289;
  }
}


void mq_poly_intt_12289_n(int16_t *a) 
{
  unsigned int start,len,j,k;
  int16_t t,zeta;

  k=0;
  for(len=1;len<=DIM_N>>1;len<<=1)
  {
    for(start=0;start<DIM_N;start=j+len)
    {
      zeta=fn_12289[k++];
      for(j=start;j<start+len;j++) 
      {
        t=a[j];
        a[j]=t+a[j+len];
        a[j]-=12289;
        a[j]+=(a[j]>>15)&12289;
        a[j+len]=t-a[j+len];
        a[j+len]=fqmul_12289(zeta,a[j+len]);
        a[j+len]+=(a[j+len]>>15)&12289;
      }
    }
  }

  for(j=0;j<DIM_N;j++)
  {
    a[j]=fqmul_12289(a[j],fn_12289[1023]);
    a[j]+=(a[j]>>15)&12289;
  }
}

void mq_poly_mul_ntt_12289(int16_t *r, int16_t *a, int16_t *b)
{
  int i;
  for(i=0;i<DIM_N;i++) 
  {
    r[i]=fqmul_12289(a[i],b[i]);
    r[i]+=(r[i]>>15)&12289;
  }
}

void mq_poly_adjoint_ntt_12289(int16_t *r, int16_t *a)
{
  int i;
  for(i=0;i<DIM_N;i++) 
  {
    r[i] = a[DIM_N-i-1];
  }
}
