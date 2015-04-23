To compile on Brutus:
module load gcc
g++ cputest.cpp -o cputest

To submit the executable:
bsub -W 00:10 -n 1 ./cputest

To submit the script that runs the executable:
bsub -W 00:10 -n 1 < script



