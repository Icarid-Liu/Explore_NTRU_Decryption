#ifndef KEYGEN_H
#define KEYGEN_H

#include <stdio.h>
#include <stdint.h>
#include "fips202.h"

void END_keygen(int8_t *h, int8_t *f, int8_t *g, int8_t *F, int8_t *G, int16_t *w, int8_t *k, void *tmp, size_t tmp_len);

#endif