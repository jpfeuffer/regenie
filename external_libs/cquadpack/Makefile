AR   ?= ar
CC   ?= gcc
FLAG  = -O2 -DNDEBUG -fPIC

OBJECTS = dqags.o dqext.o dqsort.o gk21.o

libcquadpack.a: ${OBJECTS}
	${AR} rcs $@ $^

%.o: %.c
	${CC} $(FLAG) -c $<

clean:
	rm -f *.o *.a
