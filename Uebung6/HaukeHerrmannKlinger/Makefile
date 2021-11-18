CC = mpicc
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -gdwarf-4
LFLAGS = $(CFLAGS)
TGTS = timempi

all: $(TGTS)

$(TGTS): %: %.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(TGTS)
