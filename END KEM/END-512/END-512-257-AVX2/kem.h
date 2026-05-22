#ifndef KEM_H
#define KEM_H

#include <stdint.h>

void END_cca_encrypt(int16_t *h, int8_t *c, int8_t *K);
int END_cca_decrypt(int16_t *h, int8_t *c, int8_t *K, int8_t *f, int8_t *g, int8_t *F, int8_t *G, int8_t *k, int16_t *w);
#endif