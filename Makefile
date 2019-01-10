CC = g++
CXXFLAGS = $(shell pkg-config --cflags jellyfish-2.0) -std=c++11 -Wall -O3
LDFLAGS = -Wl,--rpath=$(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L//g')
LDLIBS = $(shell pkg-config --libs jellyfish-2.0)
SRC=cFP.cpp kmer.cpp Util.cpp readerWriter.cpp
EXEC=main
OBJ= $(SRC:.cpp=.o)

#all: main

all:
	$(MAKE) clean
	$(MAKE) $(EXEC)

main: $(OBJ) main.cpp
	$(CXX) -o $@ $(OBJ) main.cpp $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -lz

%.o: %.cpp %.h
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)


%.o: %.c %.h
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -rf *.o
