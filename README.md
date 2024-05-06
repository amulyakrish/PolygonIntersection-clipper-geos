This is the code for clipper library code for (integral co-ordinates) and geos code for (double precision geometric co-ordinates) to find Polygon Intersection

Execution Commands for Polygon-Clipping library:

compile:
g++ -g -c clipp.cpp -o clipp.o -DCLIPPER_STATIC_LIB -std=c++11
g++ -g -c clipper.cpp -o clipper.o -DCLIPPER_STATIC_LIB -std=c++11
g++ -g -o tree clipp.o  clipper.o -std=c++11 -lpthread

run:
./tree

Execution Commands for Geos library:

compile:
g++ -o prog geos-test.cpp -I/usr/include -L/usr/lib/x86_64-linux-gnu -lgeos_c

run:
./prog




