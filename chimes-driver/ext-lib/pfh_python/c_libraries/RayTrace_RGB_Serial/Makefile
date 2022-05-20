#!/usr/bin/env bash

ICC  := $(shell command -v icc 2> /dev/null)

SRCS   = main.c
EXEC   = raytrace_rgb.so
CFLAGS = -shared -fPIC -O3 -g
ifdef ICC
   CC     = icc
else
   CC     = gcc
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
