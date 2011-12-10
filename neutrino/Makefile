# Makefile - EP1 - MAC5742 - 2011/2
# Gabriel A. von Winckler - 3313359

CC ?= icc

ifeq ($(CC), icc)
	CFLAGS ?= -Wall -fast -openmp
else
	CFLAGS ?= -march=native -Wall -O3 -fopenmp
endif

neutrino:
	$(CC) -o neutrino $(CFLAGS) neutrino.c -lm

clean:
	rm neutrino
