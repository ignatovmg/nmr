CC=gcc
CFLAGS=-Wall -O3 -std=gnu99 -fPIC -DMETHOD=1 -DTEST
SRCDIR=src
BUILDDIR=build
LIBS=`pkg-config --cflags --libs gsl`

all: $(SRCDIR)/main.c
	$(CC) $(CFLAGS) $< $(LIBS) -o $(BUILDDIR)/main
    
clean:
	rm -rf ./build/*
