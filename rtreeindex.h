/**
 **  This source code originally provided by:
 **     http://www.superliminal.com/sources/sources.htm
 **
 **  This code is placed in the public domain.
 **/
#ifndef __R_TREE_INDEX_HEADER__
#define __R_TREE_INDEX_HEADER__

#define BUNDLEFACTOR   96 //try changing higher
#define FANOUT    4
#define NUMDIMS    2   
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <utility>
#include <immintrin.h>
#include "segseg.h"
#include "ThreadPool.h"
#include <atomic>
extern std::vector<std::future<void>> futures;
using namespace std;
using std::vector;

typedef double RectReal;      // float vs double
extern ThreadPool pool;
extern std::mutex console_mutex;
extern std::atomic<int> count;

// typedef float RectReal;    // SIMD ps pd


/**
 ** Global definitions.
 **/

#define NUMSIDES 2*NUMDIMS

struct MBR   // Tile structure
{                // original array of vertices of a polygon
    int start;   // array index of a starting vertex 
    int end;     // array index of a ending vertex
    vector<RectReal> *boundary; /* xmin,ymin,...,xmax,ymax,... */
    
    MBR(int a, int b, vector<RectReal> *c)
    {
        start = a;
        end = b;
        boundary = c;
    }
};

struct Node       // Tile
{
    struct MBR *mbr;
    int leaf; /* 1 if leaf, 0 if not */
    vector<struct Node *> *children;
};

struct Point
{
    double x,y;
    Point(double a, double b)
    {
       x = a;
       y = b;
    }
};

double max(double num1, double num2);

double min(double num1, double num2);

int overlap(struct MBR *rect1, struct MBR *rect2);

vector<pair<struct Node *, struct Node *> > overlappingChildren(struct Node *R1, struct Node *R2);

void rectUnion(vector<RectReal> *rect1, vector<RectReal> *rect2);

void rectUnion4(vector<RectReal> *rect1, vector<RectReal> *rect2, vector<RectReal> *rect3, vector<RectReal> *rect4);

void getRange(int *range, int id, int low, int high);

struct MBR *createTile(vector<Point> *ptArr, int first, int last);

struct Node *createLeaf(vector<Point> *ptArr, int low, int high);

struct MBR *unionJoin(vector<struct Node *> *nodes, int n);

struct Node *createRTreeThreadPool(vector<Point> *ptArr, int low, int high);

struct Node *createRTreeSequential(vector<Point> *ptArr,int low, int high);

//struct Node *createRTree(vector<struct Point> *ptArr, int low, int high, int depth = 0);

int findIndex(struct Node *R1, struct Node *R2, vector<struct Point> *ptArr, int index);

void spatialJoin(struct Node *R1, struct Node *R2, vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ);

int contains(struct Point p, MBR *R1, MBR *R2);

int isOverlap(int i, int j, vector<struct Point> *P, vector<struct Point> *Q);

void lsMBRfilter(vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ);

void checkLSI(vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ);

void printNode(struct Node *root);
#endif
