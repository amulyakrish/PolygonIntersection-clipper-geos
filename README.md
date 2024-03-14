This is the parallel version of Sequential Polygon Itntersection Code

Execution Commands:

compile:
g++ -c buildTree.cpp
g++ -c -mavx rtreeindex.cpp
g++ -c readSCPolygon.cpp
g++ -c segseg.cpp
g++ -o tree buildTree.o rtreeindex.o readSCPolygon.o segseg.o -std=c++11 -lpthread


run:
./tree

Main changes include

1.Parallelizing spatial Join method to make the base case execution concurrent and independent. Assumption is that base level has more work to do in identifying polygon intersection points.

To achive this, threads, locks, futures concept for asynchronous execution have been used

2. Second optimization is reducing the Rtree Creation time, inorder to get overall speed up of finding Intersections
For this, the set of data points are dividend into "bundle factor" number of sub problems, and each sub problem is allowed to execute independently.

We cannot parallelize it more because the upper level in recursion waits for lower tiles to get grouped togethr. Due to this dependency more thread level optimization is not possible and accuracy can go for a toss.

This is version 1 of the parallelized code.

Follow up work and Other Improvements:
To make the code work for worst case scenarios.
