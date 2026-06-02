CXX      = g++
CXXFLAGS = -std=c++20 -O2 -Wall -Wextra -pedantic
THREADS  = -pthread

BASELINE_TARGET = hac_single_baseline
MST_TARGET      = hac_single_mst

.PHONY: all clean

all: $(BASELINE_TARGET) $(MST_TARGET)

$(BASELINE_TARGET): baseline.cpp
	$(CXX) $(CXXFLAGS) baseline.cpp -o $(BASELINE_TARGET)

$(MST_TARGET): parallel_mst.cpp
	$(CXX) $(CXXFLAGS) $(THREADS) parallel_mst.cpp -o $(MST_TARGET)

clean:
	rm -f $(BASELINE_TARGET) $(MST_TARGET)
