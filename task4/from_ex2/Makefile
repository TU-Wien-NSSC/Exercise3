CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math -fopenmp
SRC = main_ex2.cpp


NUM_THREADS ?= 10
SCHEDULE_MODE ?= static

all: 
	$(CXX) $(CXXFLAGS) -o solver $(SRC)

run_small: solver
	./solver reference 128 1000

run_test: solver
	OMP_NUM_THREADS=$(NUM_THREADS) OMP_SCHEDULE=$(SCHEDULE_MODE) ./solver reference 128 1000

clean:
	rm -f solver *.csv

.PHONY: all run clean


