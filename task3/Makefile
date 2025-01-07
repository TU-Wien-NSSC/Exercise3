CXX = g++
CXXFLAGS = -std=c++17
SRC = main.cpp

all: 
	$(CXX) $(CXXFLAGS) -o mcint $(SRC)

openmp: 
	$(CXX) $(CXXFLAGS) -fopenmp -o mcint $(SRC)

parallelsinx:
	OMP_NUM_THREADS=4 time ./mcint SINX 0 6.283185 1e8

# test SINX func
sinx: 
	./mcint SINX 0 6.283185 1e9

# test COS2XINV func
cos2xinv:
	./mcint COS2XINV 0 6.283185 1e9

# test X4M5 func
x4m5:
	./mcint X4M5 0 6.283185 1e9

clean:
	rm -f mcint

.PHONY: all run clean