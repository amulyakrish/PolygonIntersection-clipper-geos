/**
 **  This source code originally provided by:
 **     http://www.superliminal.com/sources/sources.htm
 **
 **  This code is placed in the public domain.
 **/
#ifndef __R_TREE_INDEX_HEADER__
#define __R_TREE_INDEX_HEADER__

#define BUNDLEFACTOR    96 //try changing higher
#define FANOUT    4 //need fanout of 4 for optimizing using floats
#define NUMDIMS    2   
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <utility>
#include <immintrin.h>
using namespace std;
using std::vector;

typedef float RectReal;


/**
 ** Global definitions.
 **/

#define NUMSIDES 2*NUMDIMS

struct MBR
{
    int start;
    int end;
    vector<RectReal> *boundary; /* xmin,ymin,...,xmax,ymax,... */
    MBR(int a, int b, vector<RectReal> *c)
    {
        start = a;
        end = b;
        boundary = c;
    }
};

struct Node
{
    struct MBR *mbr;
    int leaf; /* 1 if leaf, 0 if not */
    vector<struct Node *> *children;
};

struct Point
{
    float x,y;
    Point(float a, float b)
    {
       x = a;
       y = b;
    }
};

float max(float num1, float num2);

float min(float num1, float num2);

int* overlap(struct MBR *R1, struct MBR *R2, struct MBR *R3, struct MBR *R4, struct MBR *R5, struct MBR *R6, struct MBR *R7, struct MBR *R8);

vector<pair<struct Node *, struct Node *> > overlappingChildren(struct Node *R1, struct Node *R2);

void rectUnion(vector<RectReal> *rect1, vector<RectReal> *rect2);

void rectUnion4(vector<RectReal> *rect1, vector<RectReal> *rect2, vector<RectReal> *rect3, vector<RectReal> *rect4);

void rectUnion8(vector<RectReal> *rect1, vector<RectReal> *rect2, vector<RectReal> *rect3, vector<RectReal> *rect4, vector<RectReal> *rect5, vector<RectReal> *rect6, vector<RectReal> *rect7, vector<RectReal> *rect8);

void getRange(int *range, int id, int low, int high);

struct MBR *createTile(vector<Point> *ptArr, int first, int last);

struct Node *createLeaf(vector<Point> *ptArr, int low, int high);

struct MBR *unionJoin(vector<struct Node *> *nodes, int n);

struct Node *createRTree(vector<Point> *ptArr,int low, int high);

void spatialJoin(struct Node *R1, struct Node *R2, vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ);

void checkLSI(vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ);


void printNode(struct Node *root);
#endif