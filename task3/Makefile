CXX = g++
CXXFLAGS = -std=c++17
SRC = main.cpp
ITERS = 1e7
XMIN = -1.570796326794897
XMAX = 1.570796326794897
THREADS = 4

serial:
	$(CXX) $(CXXFLAGS) -O3 -o mcint $(SRC)

parallel:
	$(CXX) $(CXXFLAGS) -fopenmp -O3 -o mcint $(SRC)

# test SINX func
sinx: 
	OMP_NUM_THREADS=$(THREADS) ./mcint SINX $(XMIN) $(XMAX) $(ITERS)

# test COS2XINV func
cos2xinv:
	OMP_NUM_THREADS=$(THREADS) ./mcint COS2XINV $(XMIN) $(XMAX) $(ITERS)

# test X4M5 func
x4m5:
	OMP_NUM_THREADS=$(THREADS) ./mcint X4M5 $(XMIN) $(XMAX) $(ITERS)

clean:
	rm -f mcint

.PHONY: all run clean