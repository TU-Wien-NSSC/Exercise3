CXX = g++
CXXFLAGS = -std=c++17
SRC = main.cpp

all: 
	$(CXX) $(CXXFLAGS) -o cg $(SRC)

run: 
	./cg bcsstk13.mtx 1000

plot: plot_Ex1.py
	python3 plot_Ex1.py

clean:
	rm -f $(EXECUTABLE) residuals.txt PlotResidual.png PlotError.png

.PHONY: all run clean
