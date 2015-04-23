# Makefile

SHELL := /bin/bash

#######################
# USER SETUP
#######################
# use release or debug settings
config ?= release
#######################

CFLAGS = -std=c++11 -fopenmp -I../source/

ifeq "$(config)" "release"
  CFLAGS += -O3 -DNDEBUG
else
  CFLAGS += -g
endif

OBJECTS = main.o Diffusion.o

# define paths for source files
VPATH := ../source/
.DEFAULT_GOAL := diffusion

diffusion: $(OBJECTS)
	g++ $(CFLAGS) $^ -o $@

%.o: %.cpp
	g++ $(CFLAGS) -c $^ -o $@
	
clean:
	rm -f *.o
	rm -f diffusion
