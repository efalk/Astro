
#CFLAGS = -g -Wall -Werror
CFLAGS = -O

LIBS = -lm

SRCS = coords.c dates.c utils.c precession.c kepler.c sun.c planets.c \
	moon.c stars.c io.c navigation.c yale.c

OBJS = $(SRCS:.c=.o)

HDRS = astro.h

PROGS = ephem dates test navigation

lib:	libastro.a

libastro.a:	$(OBJS)
	ar -ruv $@ $(OBJS)

ephem:	ephem.o libastro.a
	$(CC) $(CFLAGS) -o ephem ephem.o libastro.a $(LIBS)

test:	test.o libastro.a
	$(CC) $(CFLAGS) -o test test.o libastro.a $(LIBS)

dates:	dates.c utils.o libastro.a
	$(CC) $(CFLAGS) -DSTANDALONE -o dates dates.c libastro.a $(LIBS)

navigation:	navigation.c libastro.a
	$(CC) $(CFLAGS) -DSTANDALONE -o navigation navigation.c libastro.a $(LIBS)

yale:	yale.c libastro.a
	$(CC) $(CFLAGS) -DSTANDALONE -o yale yale.c libastro.a $(LIBS)


tags: $(SRCS) $(HDRS)
	ctags $(SRCS) $(HDRS)

clean:
	rm -f *.o

spotless: clean
	rm -f *.a ${PROGS} tags
	rm -rf *.dSYM
