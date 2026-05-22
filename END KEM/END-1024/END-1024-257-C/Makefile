SOURCES = big_integer.c fixed_point.c keygen.c code.c pke.c kem.c\
          module_p.c mq_ntt.c ntru.c poly.c randombytes.c fips202.c
HEADERS = big_integer.h fixed_point.h keygen.h module_p.h\
          mq_ntt_param.h mq_ntt.h ntru.h param.h poly.h code.h pke.h kem.h randombytes.h fips202.h



main: $(HEADERS) $(SOURCES)
	clang -O3 -DNDEBUG -march=native -mtune=native -fomit-frame-pointer -w -o $@ $(SOURCES) main.c 
memory_usage: $(HEADERS) $(SOURCES)
	gcc -O3 -DNDEBUG -march=native -mtune=native -fomit-frame-pointer -w -o $@ $(SOURCES) main.c 
timecop: $(HEADERS) $(SOURCES)
	clang -O0 -gdwarf-4 -Wall -w -o $@ $(SOURCES) timecop.c 
clean:
	rm -f main timecop memory_usage $(objects) *.txt *.o *.exe *.out
