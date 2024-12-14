CXX = g++
CXXFLAGS = -std=c++17
SRC = main.cpp

all: 
	$(CXX) $(CXXFLAGS) -o cg $(SRC)

# test without source
run: 
	./cg bcsstk13.mtx 128

clean:
	rm -f cg

.PHONY: all run clean