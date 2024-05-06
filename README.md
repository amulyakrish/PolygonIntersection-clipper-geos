This is the clipper library code to find Polygon Intersection

Execution Commands:

compile:
g++ -g -c clipp.cpp -o clipp.o -DCLIPPER_STATIC_LIB -std=c++11
g++ -g -c clipper.cpp -o clipper.o -DCLIPPER_STATIC_LIB -std=c++11
g++ -g -o tree clipp.o  clipper.o -std=c++11 -lpthread

run:
./tree


