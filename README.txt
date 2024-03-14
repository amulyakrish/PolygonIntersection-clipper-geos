File definitions:

**Note any filename that ends in 2, as in rtreeindex2.cpp, is optimized for using float types rather than doubles.**
**any filename that ends in 2_512, is optimized for using float types rather than doubles as well as 512-bit vector instructions**

buildTree.cpp:
This is the driver program. Contains code to build polysketch and run a spatial join for c and s polygons from
GPC Polygon. Uncomment specified regions to run clipper and make comparison.

clipper.cpp:
Clipper library from http://www.angusj.com
Remember to add compilation and linking to the compilation instructions at the bottom to utilize clipper.

readSCPolygon.cpp:
Code to read in GPC Polygon.

rtreeindex.cpp:
This is where our important algorithms and definitions are. To enable threading, look at comments in createRtree() and make sure to also
uncomment the includes at the top. rtreeeindex.h is where parameters like FANOUT and BUNDLEFACTOR can be adjusted.

rtreeindex_512.cpp:
This code has the same functionality as rtreeindex.cpp with the use of 512-bit vector intstructions. 

rtreeindexNOSIMD.cpp:
This code has the same functionality as rtreeindex.cpp but without use of intrinsic SIMD instructions. For comparison.

segseg.cpp:
segement intersection code provided by Joseph O'Rourke

s.txt, c.txt:
Large polygons for testing purposes.

/*
g++ -c buildTree.cpp
g++ -c -mavx rtreeindex.cpp  -->for threaded use-->   gcc -c -mavx -std=c++11 -lpthread rtreeindex.cpp
g++ -c readSCPolygon.cpp
g++ -c segseg.cpp
g++ -o tree buildTree.o rtreeindex.o readSCPolygon.o segseg.o -->for threaded use-->  g++ -o tree  -std=c+11 -lpthread buildTree.o rtreeindex.o readSCPolygon.o segseg.o
//the last step is linking step 
*/

for the version optimized for floating point values utilize rtreeindex2.cpp and buildtree2.cpp
/*
gcc -c buildTree2.cpp
gcc -c -mavx rtreeindex2.cpp  -->for threaded use-->   gcc -c -mavx -std=c++11 -lpthread rtreeindex2.cpp
gcc -c readSCPolygon.cpp
gcc -c segseg.cpp
g++ -o tree buildTree2.o rtreeindex2.o readSCPolygon.o segseg.o -->for threaded use-->  g++ -o tree  -std=c+11 -lpthread buildTree2.o rtreeindex2.o readSCPolygon.o segseg.o
*/

