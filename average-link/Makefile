CXX      = g++
CXXFLAGS = -std=c++20 -O2 -Wall -Wextra -pthread

SEQ_TARGET   = hac_seq
NAIVE_TARGET = hac_naive
PPOP_TARGET  = hac_ppop

SEQ_SRCS    = sequential/main.cpp \
              sequential/average-link-matrix-seq.cpp

NAIVE_SRCS  = naive/main.cpp \
              naive/average-link-naive.cpp \
              sequential/average-link-matrix-seq.cpp

PPOP_SRCS   = pPOP/main.cpp \
              pPOP/average-link-pPOP.cpp \
              sequential/average-link-matrix-seq.cpp

.PHONY: all seq naive pPOP clean

all: seq naive pPOP

seq: $(SEQ_TARGET)

naive: $(NAIVE_TARGET)

pPOP: $(PPOP_TARGET)

$(SEQ_TARGET): $(SEQ_SRCS) sequential/average-link-matrix-seq.hpp
	$(CXX) $(CXXFLAGS) $(SEQ_SRCS) -o $(SEQ_TARGET)

$(NAIVE_TARGET): $(NAIVE_SRCS) naive/average-link-naive.hpp sequential/average-link-matrix-seq.hpp
	$(CXX) $(CXXFLAGS) $(NAIVE_SRCS) -o $(NAIVE_TARGET)

$(PPOP_TARGET): $(PPOP_SRCS) pPOP/average-link-pPOP.hpp sequential/average-link-matrix-seq.hpp
	$(CXX) $(CXXFLAGS) $(PPOP_SRCS) -o $(PPOP_TARGET)

clean:
	rm -f $(SEQ_TARGET) $(NAIVE_TARGET) $(PPOP_TARGET)
