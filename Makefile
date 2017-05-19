CXXSRCS := $(wildcard *.cpp)
CSRCS := $(wildcard *.c)
ASSRCS := $(wildcard *.S)
CXX := g++
CC := gcc
AS := as
CXXFLAGS := -Wall -Wextra -std=c++11 -O2 -mavx2
CFLAGS := -Wall -Wextra -std=c11 -O2 -mavx2
ASFLAGS :=
LDFLAGS := -lm
BINS := self1 self2 self3
CXXOBJS := $(CXXSRCS:%.cpp=%.o)
COBJS := $(CSRCS:%.c=%.o)
ASOBJS := $(ASSRCS:%.S=%.o)
.PHONY: all clean
all: $(BINS)
clean:
	-@rm -vf $(BINS) $(CXXOBJS) $(COBJS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
%.o: %.S
	$(AS) $(ASFLAGS) -o $@ $<
self1: self1.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
self2: self2.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
self3: self3.o biot_savart.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
