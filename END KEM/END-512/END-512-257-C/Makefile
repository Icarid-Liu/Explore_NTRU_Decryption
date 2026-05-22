SOURCES = big_integer.c fixed_point.c keygen.c code.c\
          module_p.c mq_ntt.c ntru.c pke.c poly.c randombytes.c fips202.c kem.c 
HEADERS = big_integer.h fixed_point.h keygen.h module_p.h kem.h code.h\
          mq_ntt_param.h mq_ntt.h ntru.h param.h pke.h poly.h randombytes.h fips202.h


main: $(HEADERS) $(SOURCES)
	clang $(CFLAGS) -O3 -DNDEBUG -march=native -mtune=native -fomit-frame-pointer -w -o $@ $(SOURCES) main.c
memory_usage: $(HEADERS) $(SOURCES)
	gcc $(CFLAGS) -O3 -DNDEBUG -march=native -mtune=native -fomit-frame-pointer -w -o $@ $(SOURCES) main.c
timecop: $(HEADERS) $(SOURCES)
	clang $(CFLAGS) -O0 -gdwarf-4 -Wall -w -o $@ $(SOURCES) timecop.c
clean:
	rm -f main memory_usage timecop $(objects) *.txt *.o *.exe *.out
