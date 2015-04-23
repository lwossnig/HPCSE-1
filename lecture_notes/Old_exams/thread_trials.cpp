#include <iostream>
#include <thread>
#include <mutex>

int main ()
	{
		std::mutex m;
		int count=0;
		bool excount=true;

		std::thread t([&]{
				std::lock_guard<std::mutex> lock(m);
				std::cout << "Thread counts " << ++count << std::endl;
				excount = false;
				});
		
		while(true){
			std::lock_guard<std::mutex> lock(m);
			if(excount==false)
				break;
			m.unlock();
		}

		std::lock_guard<std::mutex> lock(m);
		std::cout << "Main counts " << ++count << std::endl;
//		m.unlock();
		t.join();
		return 0;
	}
