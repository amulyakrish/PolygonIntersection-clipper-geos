#ifndef __SEG_INT_HEADER__
#define __SEG_INT_HEADER__
#include <stdio.h>
#include <math.h>
#include <vector>
#define	EXIT_FAILURE 1
#define	X 0
#define	Y 1
#define	DIM 2               /* Dimension of points */
typedef	int tPointi[DIM];   /* Type integer point */
typedef	double tPointd[DIM];   /* Type double point */
using namespace std;
using std::vector;


/*---------------------------------------------------------------------
Function prototypes.
---------------------------------------------------------------------*/

double* dot(__m256d A, __m256d B);
char* findCrossingLines(struct Point p1, struct Point p2, 
                                vector<struct Point> *p3, 
                                vector<struct Point> *p4);
char* findCrossingLines2(struct Point p1, struct Point p2, 
                                vector<struct Point> *p3, 
                                vector<struct Point> *p4);

bool lessThanIntersectionPoint(__m256d vec1, __m256d vec2, double intersection);

bool greaterThanIntersectionPoint(__m256d vec1, __m256d vec2, double intersection);

int dot(tPointi a, tPointi b);

bool SegSeg(tPointi p1, tPointi p2, tPointi p3, tPointi p4);

char SegSegInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p );
char ParallelInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p );
bool Between( tPointi a, tPointi b, tPointi c );
void Assigndi( tPointd p, tPointi a );
int	Collinear( tPointi a, tPointi b, tPointi c );
int AreaSign( tPointi a, tPointi b, tPointi c );

#endif