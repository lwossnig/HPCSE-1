#ifndef BARRIER_HPP
#define BARRIER_HPP

#include <thread>
#include <mutex>
#include <cassert>
#include <limits>

namespace diffusion{
	class barrier
	{
		public:
			barrier(unsigned int count)
				: m_total(count)
				, m_count(count)
				, m_generation(0)
				{
					assert(count !=0);
				}
			void wait()
			{
				std::unique_lock<std::mutex> lock(m_mutex);
				unsigned int gen = m_generation;

				/// decrease count each time function is called
				/// and if not last thread continues in endless loop
				if (--m_count==0){ 
					/// decreases m_count and checks if = 0
					/// if done reset to new wait-generation
					m_count = m_total;
					m_generation++;
				}
				else{
					lock.unlock(); /// m_mutex.unlock -> other mutexes may access now again
					while (true) {
						lock.lock();
						if(gen != m_generation)
							break;
						lock.unlock();
					}
				}
			}

			unsigned int num_waiting() const
			{
				std::unique_lock<std::mutex>  lock(m_mutex);
				return m_count;
			}

		private:
			mutable std::mutex m_mutex;
			///  mutable allows a const member function to change value
			unsigned long const m_total; /// number of threads
			unsigned long m_count; /// num of threads who entered in round
			unsigned long m_generation; ///	generation of wait (n-times call = n'th gen)      
	};
} // end namespace diffusion!
#endif // defined BARRIER_HPP
