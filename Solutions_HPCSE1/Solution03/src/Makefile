CXX?=g++

perf=0

ifeq "$(perf)" "1"
	CFLAGS += -D_PERF_
endif

# these are for running on Euler!
ifeq ($(CXX),/cluster/apps/gcc/4.8.2/bin/g++)
	CFLAGS+=-Wall -O3 -std=c++11
	CFLAGS_THREADS=$(CFLAGS) -fopenmp
endif

ifeq ($(CXX),/cluster/apps/intel/composer_xe_2013_sp1.1.106/composer_xe_2013_sp1.1.106/bin/intel64/icpc)
	CFLAGS+=-Wall -O3 -std=c++11 -gcc-name=/cluster/apps/gcc/4.8.2/bin/gcc
	CFLAGS_THREADS=$(CFLAGS) -openmp
endif


all: diffusion2d_serial diffusion2d_openmp

diffusion2d_serial: diffusion2d_serial.cpp
	$(CXX) $(CFLAGS) -o $@ $<

diffusion2d_openmp: diffusion2d_openmp.cpp
	$(CXX) $(CFLAGS_THREADS) -o $@ $<

clean:
	rm -f diffusion2d_serial diffusion2d_openmp
