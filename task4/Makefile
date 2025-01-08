# CXX = g++
# CXXFLAGS = -std=c++17 -fopenmp
# SRC = main_ex2.cpp

# all: 
# 	$(CXX) $(CXXFLAGS) -o solver $(SRC)

# # test without source
# nosource : 
# 	./solver without_source 128 1000 0.0 1.0

# # test with source
# source:
# 	./solver with_source 128 1000 0.0 1.0 0.5 0.5 0.1

# clean:
# 	rm -f solver *.csv

# .PHONY: all run clean



CC = g++
CFLAGS = -std=c++17 -O3 -march=native -ffast-math -fopenmp

TARGET = solver
SRC = main.cpp

NUM_THREADS ?= 10
SCHEDULE_MODE ?= static


all: solver

solver: $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o solver

run: solver
	./solver reference 2048 100000


run_small: solver
	./solver reference 128 1000

run_test: solver
	OMP_NUM_THREADS=$(NUM_THREADS) OMP_SCHEDULE=$(SCHEDULE_MODE) ./solver reference 128 1000

plot: plot.py
	python3 plot.py

clean:
	rm -f solver result_*.txt result_summary.txt *.csv