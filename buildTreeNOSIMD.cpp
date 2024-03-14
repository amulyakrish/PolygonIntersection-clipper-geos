#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "rtreeindex.h"
#include "readSCPolygon.h"
//#include "clipper.hpp"
#include <iostream>
#include <sys/time.h>
using namespace std;
//using namespace ClipperLib;
/*
g++ -c buildTree.cpp
g++ -c -mavx rtreeindex.cpp     gcc -c -mavx -std=c++11 -lpthread rtreeindex.cpp
g++ -c readSCPolygon.cpp
g++ -c segseg.cpp
g++ -o tree buildTree.o rtreeindex.o readSCPolygon.o segseg.o       g++ -o tree  -std=c+11 -lpthread buildTree.o rtreeindex.o readSCPolygon.o segseg.o
*/

//used for converting GPC polygon to our polygon
vector<struct Point> * gpc_polygon_convert(gpc_polygon *poly, double xTrans, double yTrans)
{
    int i;
    vector<struct Point> *p = new vector<struct Point>[poly->contour[0].num_vertices];

    for(i = 0; i < poly->contour[0].num_vertices; i += 1)
    {   
        p->push_back(Point((poly->contour[0].vertex[i].x)+xTrans, (poly->contour[0].vertex[i].y)+yTrans));
    }
    
    return p;
}

//used for converting GPC polygon to Clipper polygon
// void gpc_polygon_to_path(gpc_polygon *poly, Paths *convPoly)
// {
//     int i;
//     for(i = 0; i < poly->contour[0].num_vertices; i += 1)
//     {
//         (*convPoly)[0] << IntPoint((poly->contour[0].vertex[i].x),(poly->contour[0].vertex[i].y));
//     }
// }

int main()
{
    gpc_polygon subject, clip;
    FILE *sfp, *cfp, *sfp1, *cfp1, *sfp2, *cfp2, *sfp3, *cfp3;

    //uncomment for clipper comparison
    // Paths SClipper(1), CClipper(1), solution;
    // Clipper c;

    // Large test case
    sfp= fopen("s.txt", "r");
    cfp= fopen("c.txt", "r");

    sfp1= fopen("lakes/lakes_851348.txt", "r");
    cfp1= fopen("lakes/lakes-synthetic_1_1.txt", "r");

    sfp2= fopen("synth/s-synthetic_1_1.txt", "r");
    cfp2= fopen("synth/c-synthetic_100_0.txt", "r");

    sfp3= fopen("synth/s-synthetic_15_10.txt", "r");
    cfp3= fopen("synth/c-synthetic_100_0.txt", "r");

    gpc_read_polygon(sfp1, &subject);
    gpc_read_polygon(cfp1, &clip);

    cout<<"Reading polygons done "<<endl;
    
    vector<struct Point> *S = gpc_polygon_convert(&subject, 0, 0);
    vector<struct Point> *C = gpc_polygon_convert(&clip, 0, 0);

    //uncomment for clipper comparison
    // gpc_polygon_to_path(&subject, &SClipper);
    // gpc_polygon_to_path(&clip, &CClipper);
    // //vector<struct Point> *C = gpc_polygon_convert(&subject, 100, -50);

    // c.AddPaths(SClipper, ptSubject, true);
    // c.AddPaths(CClipper, ptClip, true);
    
    //int candidates = 0;

    struct timeval begin, end, beginCT, endCT, beginSP, endSP;
    int iterations = 1;
    int i;
    double time_spent_total = 0;
    double time_spent = 0;
    double time_spentCT = 0;
    double time_spentSP = 0;

    for(i = 0; i < iterations; i += 1)
    {
	gettimeofday(&begin, NULL);
	gettimeofday(&beginCT, NULL);
        struct Node *Sroot = createRTree(S, 0, (S->size())-1);
        struct Node *Croot = createRTree(C, 0, (C->size())-1);
        gettimeofday(&endCT, NULL);
	
	gettimeofday(&beginSP, NULL);
        spatialJoin(Sroot, Croot, S, C);
        gettimeofday(&endSP, NULL);
    	
	gettimeofday(&end, NULL);

        time_spent = (double) (end.tv_usec - begin.tv_usec) / 1000000 +
                                     (double) (end.tv_sec - begin.tv_sec);	

        //time_spentCT = (double) (endCT.tv_usec - beginCT.tv_usec) / 1000000 +
          //                           (double) (endCT.tv_sec - beginCT.tv_sec);

	//time_spentSP = (double) (endSP.tv_usec - beginSP.tv_usec) / 1000000 +
          //                           (double) (endSP.tv_sec - beginSP.tv_sec);

        //time_spentCT = time_spentCT*1000;
	time_spent = time_spent*1000;
        time_spent_total += time_spent;
        cout<<""<<time_spent<<endl;
    }
    printf("\npolySketch spatial join: Average elapsed time over %d iterations: %f miliseconds\n", iterations,
                    (time_spent_total/(double)iterations));


    //uncomment for clipper comparison
    // for(i = 0; i < iterations; i += 1)
    // {
    //     gettimeofday(&begin, NULL);
        
    //     c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

    //     gettimeofday(&end, NULL);
    
    //     time_spent = (double) (end.tv_usec - begin.tv_usec) / 1000000 +
    //                                  (double) (end.tv_sec - begin.tv_sec);
    //     time_spent_total += time_spent;
    // }
    // printf("\nClipper: Average elapsed time over %d iterations: %f miliseconds\n", iterations, (time_spent_total/(double)iterations)*1000);

    return 0;
}