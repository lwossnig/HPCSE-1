#include<iostream>
#include<thread>
#include<mutex>
#include<vector>
#include<condition_variable>

int main()
{

const int nthreads = 10;
std::vector<std::thread> threads(nthreads);
std::mutex m;
std::condition_variable cv;
std::vector<int> oper(nthreads,0);
std::vector<int> oper2(nthreads,0);
std::vector<int> stat(nthreads,0);


	for(int thr=0;thr<nthreads;thr++)
	{
		threads[thr]= std::thread([&,thr]()
				{
				{
				std::lock_guard<std::mutex> lk(m);
				std::cout << "This is thread number " << thr << std::endl;
				oper[thr] = thr;
				stat[thr] = 1;
				}
				cv.notify_one();
			        });
	}

	std::unique_lock<std::mutex> lk(m);
	for(int i=0;i<nthreads;i++){
		cv.wait(lk, []{return 1==stat[i];});
	}
	

	oper2.swap(oper);

	for(auto iter:oper2){
		std::cout << iter << std::endl;
	}

	return 0;

}

