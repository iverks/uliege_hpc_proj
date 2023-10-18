# This is a comment line
CC:=gcc
# CC=mpicc
# CFLAGS will be the options passed to the compiler.
CFLAGS:= -O3 -Wall -Wextra -lm 
# CFLAGS:=$(CFLAGS) -fopenmp

build: fdtd.c fdtd.h
	$(CC) $(CFLAGS) -o fdtd fdtd.c