# these are for running on Euler!

# ifeq ($(CXX),/cluster/apps/gcc/4.8.2/bin/g++)
#			CFLAGS+=-Wall -O3 -std=gnu++0x
#			CFLAGS_THREADS=$(CFLAGS) -fopenmp
# endif

ifeq ($(CXX),g++)
		CFLAGS+=-Wall -O3 -std=gnu++0x
		CFLAGS_THREADS=$(CFLAGS) -fopenmp
		CFLAGS_ASSEMBLER=$(CFLAGS) -S -march=native
endif

# ifeq ($(CXX),/cluster/apps/intel/composer_xe_2013_sp1.1.106/composer_xe_2013_sp1.1.106/bin/intel64/icpc)
#			CFLAGS+=-Wall -O3 -std=c++11 -gcc-name=/cluster/apps/gcc/4.8.2/bin/gcc
#			CFLAGS_THREADS=$(CFLAGS) -openmp
# endif

#running on laptop
ifeq ($(CXX),c++)
		CXX=g++-mp-4.8
		CFLAGS+= -Wall -O3 --std=c++0x
		CFLAGS_THREADS=$(CFLAGS) -fopenmp
		CFLAGS_ASSEMBLER= -S -march=nocona $(CFLAGS)
endif

all: flops

flops_assembler: flops_test.cpp
	$(CXX) $(CFLAGS_ASSEMBLER) $<

flops: flops_test.cpp
	$(CXX) $(CFLAGS_THREADS) -o $@ $<

flops_noomp: flops_test.cpp
	$(CXX) $(CFLAGS) -o $@ $<

bandwidth_assembler: bandwidth_test.cpp
	$(CXX) $(CFLAGS_ASSEMBLER) $<

bandwidth: bandwidth_test.cpp
	$(CXX) $(CFLAGS) -o $@ $<

diff_assembler: diffusion_test.cpp
	$(CXX) $(CFLAGS_ASSEMBLER) $<

diff: diffusion_test.cpp
	$(CXX) $(CFLAGS_THREADS) -o $@ $<

diff_noomp: diffusion_test.cpp
	$(CXX) $(CFLAGS) -o $@ $<

clean:
	rm -f flops flops_noomp flops_test.s bandwidth bandwidth_noomp bandwidth_test.s diff diff_noomp diffusion_test.s
