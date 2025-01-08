CXX = g++
CXXFLAGS = -std=c++17
SRC = mainSmart.cpp

all: 
	$(CXX) $(CXXFLAGS) -o cg $(SRC)

run: 
	./cg bcsstk13.mtx 100000

plot: plot_Ex1.py
	python3 plot_Ex1.py

clean:
	rm -f $(EXECUTABLE) residuals_and_errors.txt PlotResidual.png PlotError.png

.PHONY: all run clean
