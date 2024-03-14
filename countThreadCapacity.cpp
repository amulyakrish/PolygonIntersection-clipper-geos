#include<iostream>
#include<thread>
int main()
{
    unsigned int nThreads = std::thread::hardware_concurrency();
std::cout << "Number of concurrent threads supported: " << nThreads << std::endl;

// Adjust based on your application's nature
// For CPU-bound tasks, consider starting with nThreads or nThreads - 1 (to leave a core for system tasks)
// For I/O-bound tasks, experiment with higher counts

    return 0;
}