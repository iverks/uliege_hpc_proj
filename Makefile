# This is a comment line
# CC:=gcc
CC:=mpicc
# CFLAGS will be the options passed to the compiler.
CFLAGS:= -O3 -Wall -Wextra -lm 
# CFLAGS:=$(CFLAGS) -fopenmp

build: fdtd.c fdtd.h init.c init.h sending.c sending.h utilities.c utilities.h
	$(CC) $(CFLAGS) -o fdtd fdtd.c init.c sending.c utilities.c