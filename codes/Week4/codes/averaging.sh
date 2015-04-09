#!/bin/bash
for ((i = 1; i <= 100; i++))
do 
	./diffusion2d_serial.exe 1 1 1024 >> average_rho_1024.txt
	echo "run $i"
done  

