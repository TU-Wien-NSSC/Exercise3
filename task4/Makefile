CXX = g++
CXXFLAGS = -std=c++17 -fopenmp
SRC = main_ex2.cpp

all: 
	$(CXX) $(CXXFLAGS) -o solver $(SRC)

# test without source
nosource : 
	./solver without_source 128 1000 0.0 1.0

# test with source
source:
	./solver with_source 128 1000 0.0 1.0 0.5 0.5 0.1

clean:
	rm -f solver *.csv

.PHONY: all run clean



CC = g++
CFLAGS = -std=c++17 -O3 -march=native -ffast-math -fopenmp

TARGET = solver
SRC = main.cpp


all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET) 2048 100000

plot: plot.py
	python3 plot.py

clean:
	rm -f $(TARGET) result_*.txt result_summary.txt