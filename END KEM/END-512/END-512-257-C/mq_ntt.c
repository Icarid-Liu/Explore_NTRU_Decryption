#include <stdint.h>
#include <stdio.h>
#include "param.h"
#include "mq_ntt.h"
#include "mq_ntt_param.h"

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


void mq_poly_ntt_257(int16_t *a) 
{
  unsigned int len,start,j,k;
  int16_t t,zeta;

  k=1;
  for(len=DIM_N>>1;len>=4;len>>=1) 
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
  for(len=4;len<=DIM_N>>1;len<<=1)
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
  for(len=4;len<=DIM_N>>1;len<<=1)
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
  int c0,c1,c2,c3;
  c0 = fqmul_257(a[0],b[0]);
  c1 = fqmul_257(a[1],b[1]);
  c2 = fqmul_257(a[2],b[2]);
  c3 = fqmul_257(a[3],b[3]);

  r[0] = fqmul_257((a[1]+a[3]),(b[1]+b[3]));
  r[0] -= c1;
  r[0] -= c3;
  r[0] += c2;
  r[0] = fqmul_257(r[0],zeta);
  r[0] += c0;

  r[1] = fqmul_257((a[2]+a[3]),(b[2]+b[3]));
  r[1] -= c2;
  r[1] -= c3;
  r[1] = fqmul_257(r[1],zeta);
  r[1] += fqmul_257((a[0]+a[1]),(b[0]+b[1]));
  r[1] -= c0;
  r[1] -= c1;

  r[2] = fqmul_257(c3,zeta);
  r[2] += c1;
  r[2] += fqmul_257((a[0]+a[2]),(b[0]+b[2]));
  r[2] -= c0;
  r[2] -= c2;

  r[3] = fqmul_257((a[0]+a[3]),(b[0]+b[3]));
  r[3] -= c0;
  r[3] -= c3;
  r[3] += fqmul_257((a[1]+a[2]),(b[1]+b[2]));
  r[3] -= c1;
  r[3] -= c2;
}

void mq_poly_mul_ntt_257(int16_t *r, int16_t *a, int16_t *b)
{
    int i;
    for(i = 0; i < DIM_N/8; i++) 
    {
        basemul_257(r + 8 * i, a + 8 * i, b + 8 * i, f_257[64+i]);
        basemul_257(r + 8 * i + 4, a + 8 * i + 4, b + 8 * i + 4, -f_257[64+i]);
    }
}

void basemul_257_mq(int16_t *r, int16_t *a, int16_t *b, int16_t zeta)
{
  int c0,  c1,  c2,  c3;
  c0 = fqmul_257(a[0], b[0]);
  c1 = fqmul_257(a[1], b[1]);
  c2 = fqmul_257(a[2], b[2]);
  c3 = fqmul_257(a[3], b[3]);

  r[0] = fqmul_257((a[1] + a[3]), (b[1] + b[3]));
  r[0] -= c1;
  r[0] -= c3;
  r[0] += c2;
  r[0] = fqmul_257(r[0], zeta);
  r[0] += c0;
  r[0] = fqmul_257(r[0], 1);
  r[0] += (r[0] >> 15) & Q;

  r[1] = fqmul_257((a[2] + a[3]), (b[2] + b[3]));
  r[1] -= c2;
  r[1] -= c3;
  r[1] = fqmul_257(r[1], zeta);
  r[1] += fqmul_257((a[0] + a[1]), (b[0] + b[1]));
  r[1] -= c0;
  r[1] -= c1;
  r[1] = fqmul_257(r[1], 1);
  r[1] += (r[1] >> 15) & Q;

  r[2] = fqmul_257(c3, zeta);
  r[2] += c1;
  r[2] += fqmul_257((a[0] + a[2]), (b[0] + b[2]));
  r[2] -= c0;
  r[2] -= c2;
  r[2] = fqmul_257(r[2], 1);
  r[2] += (r[2] >> 15) & Q;

  r[3] = fqmul_257((a[0] + a[3]), (b[0] + b[3]));
  r[3] -= c0;
  r[3] -= c3;
  r[3] += fqmul_257((a[1] + a[2]), (b[1] + b[2]));
  r[3] -= c1;
  r[3] -= c2;
  r[3] = fqmul_257(r[3], 1);
  r[3] += (r[3] >> 15) & Q;
}

void mq_poly_mul_ntt_257_mq(int16_t *r, int16_t *a, int16_t *b)
{
    int i;
    for(i = 0; i < DIM_N/8; i++) 
    {
        basemul_257_mq(r + 8 * i, a + 8 * i, b + 8 * i, f_257[64+i]);
        basemul_257_mq(r + 8 * i + 4, a + 8 * i + 4, b + 8 * i + 4, -f_257[64+i]);
    }
}


int baseinv_257(int16_t *r,int16_t *a,int16_t zeta)
{
  int check;
  int t,k,det;
  int zeta2=fqmul_257(zeta,zeta);

  r[0]=fqmul_257(a[2],a[2]);
  t=fqmul_257(a[1],a[3]);
  t=fqmul_257(t,2);
  r[0]+=t;
  r[0]=fqmul_257(r[0],a[0]);
  t=fqmul_257(a[1],a[1]);
  t=fqmul_257(t,a[2]);
  r[0]-=t;
  r[0]=fqmul_257(r[0],zeta);
  t=fqmul_257(a[3],a[3]);
  t=fqmul_257(t,a[2]);
  t=fqmul_257(t,zeta2);
  r[0]-=t;
  t=fqmul_257(a[0],a[0]);
  t=fqmul_257(t,a[0]);
  r[0]-=t;


  r[1]=fqmul_257(a[1],a[2]);
  t=fqmul_257(a[0],a[3]);
  t=fqmul_257(t,2);
  r[1]-=t;
  r[1]=fqmul_257(r[1],a[2]);
  t=fqmul_257(a[1],a[1]);
  t=fqmul_257(t,a[3]);
  r[1]-=t;
  r[1]=fqmul_257(r[1],zeta);
  t=fqmul_257(a[3],a[3]);
  t=fqmul_257(t,a[3]);
  t=fqmul_257(t,zeta2);
  r[1]+=t;
  t=fqmul_257(a[0],a[0]);
  t=fqmul_257(t,a[1]);
  r[1]+=t;


  r[2]=fqmul_257(a[1],a[3]);
  r[2]=fqmul_257(r[2],2);
  t=fqmul_257(a[2],a[2]);
  r[2]-=t;
  r[2]=fqmul_257(r[2],a[2]);
  t=fqmul_257(a[3],a[3]);
  t=fqmul_257(t,a[0]);
  r[2]-=t;
  r[2]=fqmul_257(r[2],zeta);
  t=fqmul_257(a[0],a[0]);
  t=fqmul_257(t,a[2]);
  r[2]+=t;
  t=fqmul_257(a[1],a[1]);
  t=fqmul_257(t,a[0]);
  r[2]-=t;

  r[3]=fqmul_257(a[2],a[2]);
  t=fqmul_257(a[1],a[3]);
  r[3]-=t;
  r[3]=fqmul_257(r[3],a[3]);
  r[3]=fqmul_257(r[3],zeta);
  t=fqmul_257(a[1],a[1]);
  t=fqmul_257(t,a[1]);
  r[3]+=t;
  t=fqmul_257(a[0],a[2]);
  t=fqmul_257(t,a[1]);
  t=fqmul_257(t,2);
  r[3]-=t;
  t=fqmul_257(a[0],a[0]);
  t=fqmul_257(t,a[3]);
  r[3]+=t;

  det=fqmul_257(a[2],a[2]);
  t=fqmul_257(a[1],a[3]);
  t=fqmul_257(t,4);
  det-=t;
  det=fqmul_257(det,a[2]);
  det=fqmul_257(det,a[2]);
  t=fqmul_257(a[0],a[2]);
  t=fqmul_257(t,2);
  t+=fqmul_257(a[1],a[1]);
  t=fqmul_257(t,a[3]);
  t=fqmul_257(t,a[3]);
  t=fqmul_257(t,2);
  det+=t;
  det=fqmul_257(det,zeta2);
  t=fqmul_257(a[3],a[3]);
  t=fqmul_257(t,a[3]);
  t=fqmul_257(t,a[3]);
  t=fqmul_257(t,zeta2);
  t=fqmul_257(t,zeta);
  det=t-det;
  t=fqmul_257(a[0],a[0]);
  t=fqmul_257(t,a[0]);
  t=fqmul_257(t,a[0]);
  det-=t;
  t=fqmul_257(a[1],a[3]);
  t=fqmul_257(t,2);
  t+=fqmul_257(a[2],a[2]);
  t=fqmul_257(t,2);
  t=fqmul_257(t,a[0]);
  t=fqmul_257(t,a[0]);
  k=fqmul_257(a[0],a[2]);
  k=fqmul_257(k,-4);
  k+=fqmul_257(a[1],a[1]);
  k=fqmul_257(k,a[1]);
  k=fqmul_257(k,a[1]);
  t+=k;
  t=fqmul_257(t,zeta);
  det+=t;
  det=fqmul_257(det,1);
  det+=(det>>15)&257;
  det=qinv[det];

  r[0]=fqmul_257(r[0],det);
  r[1]=fqmul_257(r[1],det);
  r[2]=fqmul_257(r[2],det);
  r[3]=fqmul_257(r[3],det);

  check = (uint16_t)det;
  check = (uint32_t)(-check) >> 31;
  return check - 1;
}

//if f is invertible return 0
int mq_poly_inv_ntt_257(int16_t *r, int16_t *a)
{
  int i,check=0;
  for(i=0;i<DIM_N/8;i++) 
  {
    check+=baseinv_257(r+8*i, a+8*i,f_257[64+i]);
    check+=baseinv_257(r+8*i+4, a+8*i+4,-f_257[64+i]);
  }
  return check;
}

static  int16_t montgomery_reduce_769(int32_t a)
{
  int32_t t;
  int16_t u;

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

void mq_poly_ntt_769(int16_t *a) 
{
  unsigned int len,start,j,k;
  int16_t t,zeta;

  k=1;
  for(len=DIM_N>>1;len>=4;len>>=1) 
  {
    for(start=0;start<DIM_N;start=j+len) 
    {
      zeta=f_769[k++];
      for(j=start;j<start+len;j++) 
      {
        t=fqmul_769(zeta,a[j+len]);
        a[j+len]=a[j]-t;
        a[j]=a[j]+t;
      }
    }
  }
}

void mq_poly_intt_769(int16_t *a) 
{
  unsigned int start,len,j,k;
  int16_t t,zeta;

  k=0;
  for(len=4;len<=DIM_N>>1;len<<=1)
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
  int c0,c1,c2,c3;
  c0=fqmul_769(a[0],b[0]);
  c1=fqmul_769(a[1],b[1]);
  c2=fqmul_769(a[2],b[2]);
  c3=fqmul_769(a[3],b[3]);

  r[0]=fqmul_769((a[1]+a[3]),(b[1]+b[3]));
  r[0]-=c1;
  r[0]-=c3;
  r[0]+=c2;
  r[0]=fqmul_769(r[0],zeta);
  r[0]+=c0;

  r[1]=fqmul_769((a[2]+a[3]),(b[2]+b[3]));
  r[1]-=c2;
  r[1]-=c3;
  r[1]=fqmul_769(r[1],zeta);
  r[1]+=fqmul_769((a[0]+a[1]),(b[0]+b[1]));
  r[1]-=c0;
  r[1]-=c1;

  r[2]=fqmul_769(c3,zeta);
  r[2]+=c1;
  r[2]+=fqmul_769((a[0]+a[2]),(b[0]+b[2]));
  r[2]-=c0;
  r[2]-=c2;

  r[3]=fqmul_769((a[0]+a[3]),(b[0]+b[3]));
  r[3]-=c0;
  r[3]-=c3;
  r[3]+=fqmul_769((a[1]+a[2]),(b[1]+b[2]));
  r[3]-=c1;
  r[3]-=c2;
}

void mq_poly_mul_ntt_769(int16_t *r, int16_t *a, int16_t *b)
{
    int i;
    for(i=0;i<DIM_N/8;i++) 
    {
        basemul_769(r+8*i, a+8*i, b+8*i, f_769[64+i]);
        basemul_769(r+8*i+4, a+8*i+4, b+8*i+4, -f_769[64+i]);
    }
}

static  void base_adjoint_769(int16_t *r,int16_t *a, int16_t zeta)
{
  r[0] = a[0];
  r[1] = fqmul_769(a[3], zeta);
  r[2] = fqmul_769(a[2], zeta);
  r[3] = fqmul_769(a[1], zeta);
}

void mq_poly_adjoint_ntt_769(int16_t *r, int16_t *a)
{
    int i;
    for(i=0;i<DIM_N/8;i++) 
    {
      base_adjoint_769(r+DIM_N-8*(i+1)+4, a+8*i, f_769[64+i]);
      base_adjoint_769(r+DIM_N-8*(i+1), a+8*i+4, -f_769[64+i]);
    }
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
    a[j]=fqmul_12289(a[j],fn_12289[511]);
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
    a[j]=fqmul_12289(a[j],fn_12289[511]);
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

static  int16_t montgomery_reduce_25601(int32_t a)
{
  int32_t t;
  int16_t u;

  u=a*(-25599);
  t=(int32_t)u*25601;
  t=a-t;
  t>>=16;

  return t;
}

static  int16_t fqmul_25601(int16_t a, int16_t b) 
{
  return montgomery_reduce_25601((int32_t)a*b);
}


void mq_poly_ntt_25601(int16_t *a) 
{
  unsigned int len,start,j,k;
  int16_t t,zeta;

  k=1;
  for(len=DIM_N>>1;len>=1;len>>=1) 
  {
    for(start=0;start<DIM_N;start=j+len) 
    {
      zeta=f_25601[k++];
      for(j=start;j<start+len;j++) 
      {
        t=fqmul_25601(zeta,a[j+len]);
        t+=(t>>15)&25601;
        a[j]+=(a[j]>>15)&25601;
        a[j+len]=a[j]-t;
        a[j]=a[j]+t;
        a[j]-=25601;
      }
    }
  }

}

void mq_poly_intt_25601(int16_t *a) 
{
  unsigned int start,len,j,k;
  int16_t t,zeta;

  k=0;
  for(len=1;len<=DIM_N>>1;len<<=1)
  {
    for(start=0;start<DIM_N;start=j+len)
    {
      zeta=fn_25601[k++];
      for(j=start;j<start+len;j++) 
      {
        t=a[j];
        a[j]=t+a[j+len];
        a[j]-=25601;
        a[j]+=(a[j]>>15)&25601;
        a[j+len]=t-a[j+len];
        a[j+len]=fqmul_25601(zeta,a[j+len]);
        a[j+len]+=(a[j+len]>>15)&25601;
      }
    }
  }

  for(j=0;j<DIM_N;j++)
  {
    a[j]=fqmul_25601(a[j],fn_25601[511]);
    // a[j]+=(a[j]>>15)&25601;
  }
}

void mq_poly_mul_ntt_25601(int16_t *r, int16_t *a, int16_t *b)
{
  int i;
  for(i=0;i<DIM_N;i++) 
  {
    r[i]=fqmul_25601(a[i],b[i]);
    r[i]+=(r[i]>>15)&25601;
  }
}

void mq_poly_mul_adjoint_ntt_25601(int16_t *r, int16_t *a)
{
  int i;
  for(i=0;i<DIM_N/2;i++) 
  {
    r[i]= r[DIM_N-i]=fqmul_25601(a[i],a[DIM_N-i]);
    r[i]+=(r[i]>>15)&25601;
    r[DIM_N-i]+=(r[DIM_N-i]>>15)&25601;
  }
}