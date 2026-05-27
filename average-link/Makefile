CXX      = g++
CXXFLAGS = -std=c++20 -O2 -Wall -Wextra

TARGET = hac_seq
SRCS   = sequential/main.cpp sequential/average-link-seq.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRCS) sequential/average-link-seq.hpp
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)