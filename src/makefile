#.PHONY: all clean

PROG=bioCC
CC=g++
#CC=icpc
COPT=-c #-D_NO_HASHTABLE #-Wall
LOPT=-lpthread -lz
SRC=$(PROG).cpp common.cpp Data.cpp TxtFile.cpp Calc.cpp
HDR=$(PROG).h def.h common.h Data.h TxtFile.h Calc.h
OBJ=$(SRC:.cpp=.o)
EXEC=$(PROG)

all: $(HDR) $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LOPT) $(OBJ) -o $@
	cp $@ ..

.cpp.o:
	$(CC) $(COPT) $< -o $@

#clean:
#	rm *o
