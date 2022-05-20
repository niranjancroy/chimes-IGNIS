#!/usr/bin/env bash
#EXEC   = raytrace_rgb.so

ICC  := $(shell command -v icc 2> /dev/null)
GCC9 := $(shell command -v gcc-9 2> /dev/null)
GCC8 := $(shell command -v gcc-8 2> /dev/null)
GCC7 := $(shell command -v gcc-7 2> /dev/null)
HOST := $(shell hostname)
ifeq ($(findstring pfe,$(HOST)),pfe)
	PFE = true
endif

SRCS   = main.c
EXEC   = raytrace_rgb.so
ifdef PFE
	CFLAGS = -shared -fPIC -O3 -g -qopenmp -parallel
	CC = icc
else ifdef ICC
	CFLAGS = -shared -fPIC -O3 -g -qopenmp -mkl
	CC     = icc
else ifdef GCC9
        CFLAGS = -shared -fPIC -O3 -g -fopenmp
        CC     = gcc-9
else ifdef GCC8
	CFLAGS = -shared -fPIC -O3 -g -fopenmp
	CC     = gcc-8
else ifdef GCC7
	CFLAGS = -shared -fPIC -O3 -g -fopenmp
	CC     = gcc-8
else
	CC     = gcc
	CFLAGS = -shared -fPIC -O3 -g -fopenmp
endif

OBJS   = $(SRCS:.c=.o)
INCL   =

LNKCMD =  ld -L/usr/lib -L/usr/local/lib  -shared

LIBS   =  -lm

$(EXEC): $(OBJS)
	$(CC) -L/usr/lib -L/usr/local/lib  -shared $(OBJS) $(LIBS) $(CFLAGS)  -o $(EXEC)


$(OBJS): $(INCL)

clean:
	rm $(OBJS) $(EXEC)
