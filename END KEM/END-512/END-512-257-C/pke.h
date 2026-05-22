#ifndef PKE_H
#define PKE_H

#include <stdint.h>
#include "keygen.h"

int END_encrypt(int16_t *h, int8_t *m, int8_t *r, int16_t *c);
void END_decrypt(int16_t *h, int16_t *c, int8_t *f, int8_t *g, int8_t *F, int8_t *G, int16_t *w, int8_t *m);

#endif