#compiler
CC = g++
#CC = mpic++

#flags for compiler - try to optimize
CFLAGS = -c -Wall -g -I/usr/lib64/

#Linker settings (need -lgmpxx -lgmp for gmp, -fopenmp for omp)
LDFLAGS = -fopenmp -lprofiler

#where sources and objects are
SOURCES = src/main.cpp
OBJ = build/main.o
#final executable
EXE = galago

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $@ -Wl,--no-as-needed -fopenmp -lprofiler -Wl,--as-needed

build/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ 

clean:
	rm $(EXE) $(OBJ)
