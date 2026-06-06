CXX      = g++
CXXFLAGS = -std=c++20 -O2 -Wall -Wextra -pedantic
THREADS  = -pthread

BASELINE_TARGET = single_link_baseline
MST_TARGET      = single_link_mst

.PHONY: all clean

all: $(BASELINE_TARGET) $(MST_TARGET)

$(BASELINE_TARGET): single-link-baseline.cpp
	$(CXX) $(CXXFLAGS) single-link-baseline.cpp -o $(BASELINE_TARGET)

$(MST_TARGET): single-link-mst.cpp
	$(CXX) $(CXXFLAGS) $(THREADS) single-link-mst.cpp -o $(MST_TARGET)

clean:
	rm -f $(BASELINE_TARGET) $(MST_TARGET)
