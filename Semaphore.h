#include <mutex>
#include <condition_variable>
#include <thread>
#include <iostream>

using namespace std;

class Semaphore 
{
public:
    Semaphore (int count_ = 0)
        : count(count_) {}

    inline void notify()
    {
	    unique_lock<std::mutex> lock(mtx);
	    count++;
	    cv.notify_one();
	}

    inline void wait()
    {
	    unique_lock<std::mutex> lock(mtx);
	
	    while(count == 0)
		{
			cv.wait(lock);
			        
		}
	    count--;
	}

private:
    mutex mtx;
    condition_variable cv;
    int count;
};
