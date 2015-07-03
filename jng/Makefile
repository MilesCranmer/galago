#compiler
#CC = g++
CC = mpic++

#flags for compiler - try to optimize
CFLAGS = -c -Wall -g -O 

#Linker settings
LDFLAGS = 

#where sources and objects are
SOURCES = src/main.cpp
OBJ = build/main.o
#final executable
EXE = jng

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o $@

build/%.o: src/%.cpp
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm $(EXE) $(OBJ)
