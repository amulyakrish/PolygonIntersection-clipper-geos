CC = g++

CFLAGS = -O3 -Wall -c

LFLAGS1 = -mavx512f

LFLAGS2 = -mavx

all: rtreeindex.o buildTree.o readSCPolygon.o segseg.o
		$(CC) -O3 $(LFLAGS2) -o spjoin rtreeindex.o buildTree.o readSCPolygon.o segseg.o
     
rtreeindex.o: rtreeindex.cpp
		$(CC) $(CFLAGS) $(LFLAGS2) rtreeindex.cpp

buildTree.o: buildTree.cpp
		$(CC) $(CFLAGS) $(LFLAGS2) buildTree.cpp
	
readSCPolygon.o: readSCPolygon.cpp
		$(CC) $(CFLAGS) readSCPolygon.cpp
	
segseg.o: segseg.cpp
		$(CC) $(CFLAGS) $(LFLAGS2) segseg.cpp
	
clean:
		rm *o spjoin