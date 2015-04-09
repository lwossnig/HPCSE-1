#include "timer.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <future>
#include <cmath>
#include <sstream> //stringstream
#include <iomanip> //setprecision
#include <semaphore.h> //Barrier


double D; //diffusion coefficient
unsigned int N; //domain resolution
double dT; //delta T
unsigned int nthreads;
double time_point;

// std::mutex io_mutex; // global

// struct sync
// {
// sync( std::ostream& os ) : os(os)
// , lock(io_mutex) {}
//   template <class T>
//   std::ostream& operator<<(T const& x)
//   {
//     return os << x;
//   }
// private:
// std::ostream& os; std::lock_guard<std::mutex> lock;
// };

//generate a matrix sizeXsize containing the initial distribution
double** initial_distribution(unsigned int size)
{
    double** res = 0;
    res = new double*[size];
    for (int i = 0; i < size; i++)
        {
            res[i]= new double[size];
            for (int j = 0; j < size; j++)
                {
                    if(std::abs(i*2.0/size-1.0)<0.5 && std::abs(j*2.0/size-1.0)<0.5)
                        res[i][j]=1;
                    else
                        res[i][j]=0;
                }
        }
    return res;
}

//update values in submatrix from row 'start' to row 'end'
void update_matrix(double** matrix, double** prev,int start,int end){
    for(int x=start; x<= end; x++) { //Dirichlet boundary conditions
        for(int y=1; y< N-1; y++) { //
            matrix[x][y]=prev[x][y]+ dT*D/4*N*N* //(2/N)^2
                (prev[x+1][y]+prev[x-1][y]+prev[x][y+1]+prev[x][y-1]-4*prev[x][y]);
        }
    }
}

// double** compute_online_update(double time,int nthreads) {
//     double** res = initial_distribution(N); //initial distribution at time 0
//     double** prev = initial_distribution(N); //previous state
//     for(double t=0+dT; t<= time; t+=dT) {
//         update_matrix(res,prev,1,N-1);
//         for(int y=1; y<N-1; y++){ //update the last row of prev
//         prev[N-1][y]=res[N-1][y];
//         }
//     }
//     return res;
// }

//compute the result
double** compute(double time,int nthreads) {
    double** res = initial_distribution(N); //initial distribution at time 0
    double** prev = initial_distribution(N); //previous state
    for(double t=0+dT; t<= time; t+=dT) {
        std::vector<std::thread> threads;
        for(int reminder=(N-2)%nthreads,
                interval = (N-2)/nthreads,
                increment = 1,
                start = 1; start<N-1; reminder--,start+=increment+interval) { //divide the space in nthreads intervals
            if(reminder==0) increment=0;
            threads.push_back(std::thread(update_matrix,res,prev,start,start+increment+interval-1));
        }
        for (std::thread& t : threads){
            t.join();}
        //set prev to be the current step
        for(int x=0; x< N; x++) {
            for(int y=0; y< N; y++) {
                prev[x][y]=res[x][y];
            }
        }
    }
    return res;
}

int main(int argc, char* argv[]) {
    switch(argc) {
    case 1:
        //defaults
        D = 1.0;
        N = 128;
        dT = 0.00001;
        nthreads = 12;
        time_point = 0.05;
        break;
    case 6:
        sscanf(argv[1],"%lf",&D);
        sscanf(argv[2],"%d",&N);
        sscanf(argv[3],"%lf",&dT);
        if(dT>1) { // convert integer to fraction
            dT=1.0/dT;
        }
        sscanf(argv[4],"%d",&nthreads);
        sscanf(argv[5],"%lf",&time_point);
        break;
    default:
        std::cout << "Usage: " << argv[0] << " <diff_coeff> <domain_radius> <dT> <nThreads> <Time>" << std::endl;
        return 0;
        break;
    }
    std::cout << "Executing program, nThreads: " << nthreads << " Coeff. D: " << D << " resolution: " << N << " Delta t: " << dT << std::endl;
    std::stringstream filename;
    filename << "./result_N" << N << "_T" << std::setprecision(1) << time_point << "_NT" << nthreads << ".csv";
    std::ofstream fs(filename.str());//,std::ios::app);
    if(!fs)
        {
            std::cerr<<"Cannot open the output file."<<std::endl;
            return 1;
        }

    double** res=0;
    timer t;
    std::cout << "params are D: " << D << " N: " << N << " dT: " << dT << " nthreads: " << nthreads << " time point: " << time_point << std::endl;
    t.start();
    res=compute(time_point,nthreads);
    t.stop();
    std::cout << t.get_timing() << std::endl;
    //report the timing
    std::stringstream filename_t;
    filename_t << "./time_N" << N << "_T" << std::setprecision(1) << time_point << "_NT" << nthreads << ".csv";
    std::ofstream rt(filename_t.str());//,std::ios::app);
    if(!rt)
        {
            std::cerr<<"Cannot open the output file."<<std::endl;
            return 1;
        }
    rt << D << "," << N << "," << dT << "," << nthreads << "," << time_point << "," << t.get_timing() << std::endl;
    rt.close();
    //fs << nthreads << "," << res << "," << t.get_timing() << std::endl;
    for(int x=0; x< N; x++) {
        for(int y=0; y< N; y++) {
            fs << res[x][y];
            if(y<N-1) fs << ",";
        }
        fs << std::endl;
    }
    fs.close();
}
