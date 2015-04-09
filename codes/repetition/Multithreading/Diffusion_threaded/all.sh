#!/bin/bash
g++ -std=c++11 -o diffusion_threaded diffusion_threaded.cpp
./diffusion_threaded 1 2 128 0.00001
python plot_gif.py
gnuplot plot.gnu
eog animate.gif
