XX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math -fopenmp

TARGET = solver
SRC = main.cpp

NUM_THREADS ?= 10
SCHEDULE_MODE ?= static


all: solver

solver: $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o solver

run: solver
	./solver reference 2048 100000

run_small: solver
	./solver reference 128 1000

run_test: solver
	OMP_NUM_THREADS=$(NUM_THREADS) OMP_SCHEDULE=$(SCHEDULE_MODE) ./solver reference 128 1000

clean:
	rm -f solver *.csv