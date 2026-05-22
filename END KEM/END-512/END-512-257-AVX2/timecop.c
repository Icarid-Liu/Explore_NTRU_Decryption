#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "keygen.h"
#include "kem.h"
#include "param.h"
#include "poison.h"


static void * xmalloc(size_t len)
{
	void *buf;

	if (len == 0) {
		return NULL;
	}
	buf = aligned_alloc(32,len);
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
	int8_t *h = xmalloc(PKLEN * sizeof(int8_t));
	int8_t *f = xmalloc(DIM_N * sizeof(int8_t));
	int8_t *g = xmalloc(DIM_N * sizeof(int8_t));
	int8_t *F = xmalloc(DIM_N * sizeof(int8_t));
	int8_t *G = xmalloc(DIM_N * sizeof(int8_t));
	int8_t *k = xmalloc(DIM_N * sizeof(int8_t));
	int16_t *w = xmalloc(DIM_N * sizeof(int32_t));
	int16_t *wp = xmalloc(DIM_N * sizeof(int16_t));
	int8_t *c = xmalloc(CIPHLEN);
	int8_t *K = xmalloc(SECRETLEN * sizeof(int8_t));
	void *tmpbytes = xmalloc(24 * DIM_N + 18);
	uint32_t *tmp = (void *)(((uintptr_t)tmpbytes + 7) & ~(uintptr_t)7);

	END_keygen(h, f, g, F, G, w, k, tmp, 22 * DIM_N);
	
	poison(f, DIM_N * sizeof(int8_t));
	poison(g, DIM_N * sizeof(int8_t));
	poison(F, DIM_N * sizeof(int8_t));
	poison(G, DIM_N * sizeof(int8_t));
	poison(w, DIM_N * sizeof(int32_t));
	poison(wp, DIM_N * sizeof(int16_t));
	poison(k, DIM_N * sizeof(int8_t));

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
