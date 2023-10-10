# This is a comment line
CC=gcc
# CFLAGS will be the options passed to the compiler.
CFLAGS= -O3 -Wall -Wextra -lm

build: fdtd.c fdtd.h
	$(CC) $(CFLAGS) -o fdtd fdtd.c