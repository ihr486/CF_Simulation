CXXSRCS := $(wildcard *.cpp)
CXXOBJS := $(CXXSRCS:%.cpp=%.o)
EXE := sim
CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++11 -O2
LDFLAGS := -lm
.PHONY: all clean
all: $(EXE)
clean:
	-@rm -vf $(EXE) $(CXXOBJS)
$(EXE): $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
