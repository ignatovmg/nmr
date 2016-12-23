CC=gcc
CFLAGS=-Wall -std=gnu99 -fPIC -DTEST
SRCDIR=src
BUILDDIR=build
LIBS=

all:
	$(CC) $(CFLAGS) $(LIBS) $(SRCDIR)/main.c `pkg-config --cflags --libs gsl` -o $(BUILDDIR)/main
	chmod 700 $(BUILDDIR)/main
    
clean:
	rm -rf ./build/*
	
run:
	$(BUILDDIR)/main
