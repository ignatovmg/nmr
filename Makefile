CC=gcc
CFLAGS=-Wall -O3 -std=gnu99 -fPIC -DTEST
SRCDIR=src
BUILDDIR=build
LIBS=`pkg-config --cflags --libs gsl` -lmol2 -lm

all: $(SRCDIR)/main.c $(SRCDIR)/spring_pair_min.c $(SRCDIR)/noe_min.c 
	$(CC) $(CFLAGS) $(SRCDIR)/main.c $(SRCDIR)/noe.c $(LIBS) -o $(BUILDDIR)/main
	$(CC) $(CFLAGS) $(SRCDIR)/spring_pair_min.c $(LIBS) -o $(BUILDDIR)/spring_pair_min
	$(CC) $(CFLAGS) $(SRCDIR)/noe_min.c $(SRCDIR)/noe.c $(LIBS) -o $(BUILDDIR)/noe_min
    
clean:
	rm -rf ./build/*
