CXX = g++
CXXFLAGS = -std=c++17
SRC = main.cpp

all: 
	$(CXX) $(CXXFLAGS) -o cg $(SRC)

run: 
	./cg bcsstk11.mtx 1000

plot: plot_Ex2.py
	python3 plot_Ex2.py

clean:
	rm -f $(EXECUTABLE) residual_noprecond.txt residual_withprecond.txt PlotResPrecond.png

.PHONY: all run clean
