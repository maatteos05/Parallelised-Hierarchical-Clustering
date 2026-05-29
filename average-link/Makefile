CXX      = g++
CXXFLAGS = -std=c++20 -O2 -Wall -Wextra

# --- Target binaries ---
SEQ_TARGET  = hac_seq
PPOP_TARGET = hac_ppop

# --- Source groups ---
SEQ_SRCS    = sequential/main.cpp \
              sequential/average-link-matrix-seq.cpp

PPOP_SRCS   = pPOP/main.cpp \
              pPOP/average-link-pPOP.cpp \
              sequential/average-link-matrix-seq.cpp

.PHONY: all seq pPOP clean

all: seq pPOP

seq: $(SEQ_TARGET)

pPOP: $(PPOP_TARGET)

$(SEQ_TARGET): $(SEQ_SRCS) sequential/average-link-matrix-seq.hpp
	$(CXX) $(CXXFLAGS) $(SEQ_SRCS) -o $(SEQ_TARGET)

$(PPOP_TARGET): $(PPOP_SRCS) pPOP/average-link-pPOP.hpp sequential/average-link-matrix-seq.hpp
	$(CXX) $(CXXFLAGS) $(PPOP_SRCS) -o $(PPOP_TARGET)

clean:
	rm -f $(SEQ_TARGET) $(PPOP_TARGET)