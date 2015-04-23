CXX?=g++

CFLAGS=-Wall -O3 -std=c++11
CFLAGS_THREADS=$(CFLAGS) -pthread

all: barrier_test diffusion2d_serial diffusion2d_threaded diffusion2d_barrier

barrier_test: barrier_test.cpp
	$(CXX) $(CFLAGS_THREADS) -o $@ $<

diffusion2d_serial: diffusion2d_serial.cpp
	$(CXX) $(CFLAGS) -o $@ $<

diffusion2d_threaded: diffusion2d_threaded.cpp
	$(CXX) $(CFLAGS_THREADS) -o $@ $<

diffusion2d_barrier: diffusion2d_barrier.cpp
	$(CXX) $(CFLAGS_THREADS) -o $@ $<

clean:
	rm -f barrier_test diffusion2d_serial diffusion2d_threaded diffusion2d_barrier
