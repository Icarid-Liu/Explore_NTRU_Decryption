#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "param.h"
#include "keygen.h"
#include "kem.h"
#include "poison.h"


static void * xmalloc(size_t len)
{
	void *buf;

	if (len == 0) {
		return NULL;
	}
	buf = malloc(len);
	if (buf == NULL) {
		fprintf(stderr, "memory allocation error\n");
		exit(EXIT_FAILURE);
	}
	return buf;
}

static void xfree(void *buf)
{
	if (buf != NULL) {
		free(buf);
	}
}

int main()
{
	int8_t *h = xmalloc(PKLEN);
	int8_t *f = xmalloc(DIM_N);
	int8_t *g = xmalloc(DIM_N);
	int8_t *F = xmalloc(DIM_N);
	int8_t *G = xmalloc(DIM_N);
	int8_t *k = xmalloc(DIM_N);
	int16_t *w = xmalloc(DIM_N * sizeof(*w));
	int16_t *wp = xmalloc(DIM_N * sizeof(*wp));
	int8_t *c = xmalloc(CIPHLEN);
	int8_t *K = xmalloc(SECRETLEN);
	void *tmpbytes = xmalloc(22 * DIM_N + 18);
	uint32_t *tmp = (void *)(((uintptr_t)tmpbytes + 7) & ~(uintptr_t)7);

	END_keygen(h, f, g, F, G, w, k, tmp, 22 * DIM_N);
	
    // sign the secret data
	poison(f, sizeof(int8_t) * DIM_N);
	poison(g, sizeof(int8_t) * DIM_N);
	poison(F, sizeof(int8_t) * DIM_N);
	poison(G, sizeof(int8_t) * DIM_N);
	poison(w, sizeof(int16_t) * DIM_N);
	poison(wp, sizeof(int16_t) * DIM_N);
	poison(k, sizeof(int8_t) * DIM_N);

	END_cca_encrypt(h, c, K);
	END_cca_decrypt(h, c, K, f, g, F, G, k, w);

	xfree(h);
	xfree(f);
	xfree(g);
	xfree(F);
	xfree(G);
	xfree(k);
	xfree(w);
	xfree(wp);
	xfree(c);
	xfree(K);
	xfree(tmpbytes);

	return 0;
}
