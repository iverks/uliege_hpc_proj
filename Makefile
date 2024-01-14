# This is a comment line
CC:=clang
# CC:=mpicc
# CFLAGS will be the options passed to the compiler.
CFLAGS:= -O3 -Wall -Wextra -lm 
CFLAGS:=$(CFLAGS) -g # Add debugging symbols
CFLAGS:=$(CFLAGS) -fopenmp 
CFLAGS:=$(CFLAGS) -fopenmp-targets=nvptx64-nvidia-cuda

build: main.c data.c data.h output.c output.h parameters.c parameters.h simulation.c simulation.h utilities.c utilities.h types.h
	$(CC) $(CFLAGS) -o fdtd main.c data.c output.c parameters.c simulation.c utilities.c

score: main.c data.c data.h output.c output.h parameters.c parameters.h simulation.c simulation.h utilities.c utilities.h types.h
	scorep $(CC) $(CFLAGS) -o fdtd_instr main.c data.c output.c parameters.c simulation.c utilities.c