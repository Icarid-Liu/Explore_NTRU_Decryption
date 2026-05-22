#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mq_ntt.h"
#include "param.h"
#include "keygen.h"
#include "poly.h"
#include "mq_ntt_param.h"
#include "kem.h"
#include "code.h"
#include "randombytes.h"

#define NTEST 100000

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
	int16_t i;
	int8_t *h = xmalloc(PKLEN * sizeof *h);
	int8_t *f = xmalloc(DIM_N);
	int8_t *g = xmalloc(DIM_N);
	int8_t *F = xmalloc(DIM_N);
	int8_t *G = xmalloc(DIM_N);
	int8_t *k = xmalloc(DIM_N);
	int16_t *w = xmalloc(DIM_N * sizeof *w);
	int8_t *c = xmalloc(CIPHLEN);
	int8_t *K = xmalloc(SECRETLEN);
	void *tmpbytes = xmalloc(22 * DIM_N + 18);
	uint32_t *tmp = (void *)(((uintptr_t)tmpbytes + 7) & ~(uintptr_t)7);
	uint64_t t0, t1, sum;

	// test correctness
	printf("================test correctness================\n\n");
	END_keygen(h, f, g, F, G, w, k, tmp, 22 * DIM_N);

	END_cca_encrypt(h, c, K);
	printf("--------END.KEM.Enc Shared Key--------\n");
    for(int i = 0; i < SECRETLEN; i++)
    {
        printf("%d,", K[i]);
    }
    printf("\n\n");
	END_cca_decrypt(h, c, K, f, g, F, G, k, w);
	printf("--------END.KEM.Dec Shared Key--------\n");
    for(int i = 0; i < SECRETLEN; i++)
    {
        printf("%d,", K[i]);
    }
    printf("\n\n");

	//test cpucycles
	printf("=================test cpucycles=================\n\n");
	sum = 0;
	for(i = 0; i < NTEST; i++)
	{
		t0 = __rdtsc();
		END_keygen(h, f, g, F, G, w, k, tmp, 22 * DIM_N);
		t1 = __rdtsc();
		sum += (t1-t0);
	}
	printf("END.KEM.KeyGen: %llu cycles\n",(sum / NTEST));

	END_keygen(h, f, g, F, G, w, k, tmp, 24 * DIM_N);
	sum = 0;
	for(i = 0; i < NTEST; i++)
	{
		t0 = __rdtsc();
		END_cca_encrypt(h, c, K);
		t1 = __rdtsc();
		sum += (t1-t0);
	}
	printf("\nEND.KEM.Enc:    %llu cycles\n",(sum / NTEST));

	END_keygen(h, f, g, F, G, w, k, tmp, 24 * DIM_N);
	END_cca_encrypt(h, c, K);
	sum = 0;
	for(i = 0; i < NTEST; i++)
	{
		t0 = __rdtsc();
		END_cca_decrypt(h, c, K, f, g, F, G, k, w);
		t1 = __rdtsc();
		sum += (t1-t0);
	}
	printf("\nEND.KEM.Dec:    %llu cycles\n",(sum / NTEST));

	
	xfree(h);
	xfree(f);
	xfree(g);
	xfree(F);
	xfree(G);
	xfree(k);
	xfree(w);
	xfree(c);
	xfree(K);
	xfree(tmpbytes);


	return 0;
}
