# This is a comment line
# CC:=gcc
CC:=mpicc
# CFLAGS will be the options passed to the compiler.
CFLAGS:= -O3 -Wall -Wextra -lm 
CFLAGS:=$(CFLAGS) -g # Add debugging symbols
# CFLAGS:=$(CFLAGS) -fopenmp

build: main.c buffer.c buffer.h data.c data.h output.c output.h parameters.c parameters.h simulation.c simulation.h utilities.c utilities.h types.h
	$(CC) $(CFLAGS) -o fdtd main.c buffer.c data.c output.c parameters.c simulation.c utilities.c