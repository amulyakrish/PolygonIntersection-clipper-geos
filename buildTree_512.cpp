#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "rtreeindex_512.h"
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
    FILE *sfp, *cfp;

    //uncomment for clipper comparison
    // Paths SClipper(1), CClipper(1), solution;
    // Clipper c;

    // Simple test case
    
    sfp= fopen("s.txt", "r");
    cfp= fopen("c.txt", "r");

    //sfp= fopen("synth/s-synthetic_1_1.txt", "r");
    //cfp= fopen("synth/c-synthetic_100_0.txt", "r");

    //sfp= fopen("synth/s-synthetic_15_10.txt", "r");
    //cfp= fopen("synth/c-synthetic_100_0.txt", "r");

    gpc_read_polygon(sfp, &subject);
    gpc_read_polygon(cfp, &clip);

    cout<<"Reading polygons done "<<endl;
    
    vector<struct Point> *S = gpc_polygon_convert(&subject, 0, 0);
    vector<struct Point> *C = gpc_polygon_convert(&clip, 0, 0);

    //uncomment for clipper comparison
    // gpc_polygon_to_path(&subject, &SClipper);
    // gpc_polygon_to_path(&clip, &CClipper);
    // //vector<struct Point> *C = gpc_polygon_convert(&subject, 100, -50);

    // c.AddPaths(SClipper, ptSubject, true);
    // c.AddPaths(CClipper, ptClip, true);
    

    struct timeval begin, end;
    int iterations = 100;
    int i;
    double time_spent_total = 0;
    double time_spent = 0;

    for(i = 0; i < iterations; i += 1)
    {
        gettimeofday(&begin, NULL);
        struct Node *Sroot = createRTree(S, 0, (S->size())-1);
        struct Node *Croot = createRTree(C, 0, (C->size())-1);
        spatialJoin(Sroot, Croot, S, C);
        gettimeofday(&end, NULL);
    
        time_spent = (double) (end.tv_usec - begin.tv_usec) / 1000000 +
                                     (double) (end.tv_sec - begin.tv_sec);
        time_spent_total += time_spent;
        cout<<"Time spent per iteration: "<<time_spent*1000<<endl;
    }
    printf("\npolySketch spatial join: Average elapsed time over %d iterations: %f miliseconds\n", iterations,
                    (time_spent_total/(double)iterations)*1000);

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


    //checkLSI(S, C);

    return 0;
}