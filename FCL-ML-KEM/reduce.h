#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

#define MONT -1044 // 2^16 mod q
#define QINV -3327 // q^-1 mod 2^16

#define montgomery_reduce KYBER_NAMESPACE(montgomery_reduce)
int16_t montgomery_reduce(int32_t a);

#define barrett_reduce KYBER_NAMESPACE(barrett_reduce)
int16_t barrett_reduce(int16_t a);

#define abs_16 KYBER_NAMESPACE(abs_16)
int16_t abs_16(int16_t x);

#define abs_32 KYBER_NAMESPACE(abs_32)
int32_t abs_32(int32_t x);


#endif
