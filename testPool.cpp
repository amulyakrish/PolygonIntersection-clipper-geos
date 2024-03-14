#include "ThreadPool.h"
#include<iostream>

// Example usage:
void exampleTask(int a) {
    std::cout << "Processing " << a << std::endl;
}

int main() {
    ThreadPool pool(4); // Create a ThreadPool with 4 threads

    /*
    In this specific case, the exampleTask function does not return a value, 
    so the primary use of future.get() here is to block the main thread until the task completes. 
    This is done immediately after each task is enqueued, effectively making the enqueuing process 
    synchronous for demonstration purposes. However, in typical usage, you might enqueue multiple 
    tasks without immediately waiting for each to complete, leveraging the full asynchronous 
    capabilities of the thread pool.
    */
    // Enqueue a bunch of tasks
    for(int i = 0; i < 10; ++i) {
        auto future = pool.enqueue(exampleTask, i);
        future.get(); //without this order of printing will be jumbled
    }

    return 0; // The ThreadPool will wait for all tasks to finish before exiting
}