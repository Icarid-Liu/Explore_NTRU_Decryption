CFLAGS += -O3 -mavx2 -DNDEBUG -march=native -mtune=native -fomit-frame-pointer -w
SOURCES = big_integer.c fixed_point.c keygen.c kem.c code.c\
          module_p.c mq_ntt_12289_25601.c ntru.c pke.c poly.c mq_ntt_257_769.S mq_poly_inv_ntt_257.c randombytes.c fips202.c
HEADERS = big_integer.h fixed_point.h keygen.h module_p.h kem.h code.h\
          mq_ntt_param.h mq_ntt.h ntru.h param.h pke.h poly.h randombytes.h fips202.h


main: $(HEADERS) $(SOURCES)
	clang -O3 -mavx2 -DNDEBUG -march=native -mtune=native -fomit-frame-pointer -w -o $@ $(SOURCES) main.c
memory_usage: $(HEADERS) $(SOURCES)
	gcc -O3 -mavx2 -DNDEBUG -march=native -mtune=native -fomit-frame-pointer -w -o $@ $(SOURCES) main.c
timecop: $(HEADERS) $(SOURCES)
	clang -O0 -mavx2 -march=native -gdwarf-4 -Wall -w -o $@ $(SOURCES) timecop.c
clean:
	rm -f main memory_usage timecop $(objects) *.txt *.o *.exe *.out
