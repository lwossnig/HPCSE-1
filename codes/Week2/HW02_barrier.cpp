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

std::mutex io_mutex; // global
struct sync
{
sync( std::ostream& os ) : os(os)
, lock(io_mutex) {}
  template <class T>
  std::ostream& operator<<(T const& x)
  {
    return os << x;
  }
private:
std::ostream& os; std::lock_guard<std::mutex> lock;
};

class Barrier {
private:
    unsigned int counter;
    unsigned int init_val;
    std::mutex mutx;
    std::condition_variable cv;
    bool stop;

public:
    Barrier(unsigned int nthread){
        counter=nthread;
        init_val=nthread;
        stop=true;
    }
    ~Barrier() {}
    void wait() {
        while(!stop) //suspend until all threads have resumed
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        std::unique_lock<std::mutex> ul(mutx);
        if(stop && --counter==0) { //check stop first to avoid side effects of -- while resuming threads
            stop=false;
            counter++; //this thread will not wait at the barrier
            cv.notify_all();
        }
        else {
            while(stop){
                cv.wait(ul);}
            if(++counter==init_val && !stop){ //reincrement counter, the first thread is was not waiting at the barrier
                stop=true; //reset barrier
            }
        }
    }
};

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

//update submatrix from row 'start' to row 'end' until time_point is reached
void update_matrix(double** matrix, double** prev,int start,int end,Barrier* b){
    for(double t=0+dT; t<= time_point; t+=dT) {
        //update the values in the cells
        for(int x=start; x<= end; x++) { //Dirichlet boundary conditions
            for(int y=1; y< N-1; y++) { //
                matrix[x][y]=prev[x][y]+ dT*D/4*N*N* //(2/N)^2
                    (prev[x+1][y]+prev[x-1][y]+prev[x][y+1]+prev[x][y-1]-4*prev[x][y]);
            }
        }
        if(nthreads!=1) b->wait(); //wait until all threads finished their computations
        //set prev to be the current step
        for(int x=start; x<=end; x++) {
            for(int y=1; y< N-1; y++) {
                prev[x][y]=matrix[x][y];
            }
        }
        if(nthreads!=1) b->wait(); //wait until all threads finished updating the matrix
        // --- debug: print error if update went wrong ---
        // for(int x=0; x<N; x++) { //
        //     for(int y=0; y<N; y++) {
        //         if(prev[x][y]!=matrix[x][y])
        //             sync(std::cout) << "cell "<<x<<"-"<<y<<" differs"<<std::endl;
        //     }
        //                         }
        // if(nthreads!=1) b->wait();
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

double** compute() {
    double** res = initial_distribution(N); //initial distribution at time 0
    double** prev = initial_distribution(N); //previous state
    Barrier* b = new Barrier(nthreads);
    std::vector<std::thread> threads;
    for(int reminder=(N-2)%nthreads,
            interval = (N-2)/nthreads,
            increment = 1,
            start = 1; start<N-1; reminder--,start+=increment+interval) { //divide the space in nthreads intervals
        if(reminder==0) increment=0;
        threads.push_back(std::thread(update_matrix,res,prev,start,start+increment+interval-1,b));
    }
    for (std::thread& t : threads){
        t.join();}
    return res;
}

int main(int argc, char* argv[]) {
    switch(argc) {
    case 1:
        //defaults
        D = 1.0;
        N = 128;
        dT = 0.00001;
        nthreads = 6;
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
    filename << "./result_barrier_N" << N << "_T" << std::setprecision(1) << time_point << "_NT" << nthreads << ".csv";
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
    res=compute();
    t.stop();
    std::cout << t.get_timing() << std::endl;
    //report the timing
    std::stringstream filename_t;
    filename_t << "./time_barrier_N" << N << "_T" << std::setprecision(1) << time_point << "_NT" << nthreads << ".csv";
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
