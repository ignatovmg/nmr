CC=gcc
CFLAGS=-Wall -O3 -std=gnu99 -fPIC -DTEST -g -pg
SRCDIR=src
BUILDDIR=build
LIBS=`pkg-config --cflags --libs gsl` -lmol2 -lm

all: $(SRCDIR)/main.c $(SRCDIR)/spring_pair_min.c
	$(CC) $(CFLAGS) $(SRCDIR)/main.c $(LIBS) -o $(BUILDDIR)/main
	$(CC) $(CFLAGS) $(SRCDIR)/spring_pair_min.c -lmol2 -lm -o $(BUILDDIR)/spring_pair_min
    
clean:
	rm -rf ./build/*
