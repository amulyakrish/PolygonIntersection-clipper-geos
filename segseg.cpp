/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.  It is not written to be comprehensible without the 
explanation in that book.
Compile:  gcc -o segseg segseg.c (or simply: make)
Written by Joseph O'Rourke.
Last modified: November 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely 
redistributed in its entirety provided that this copyright notice is 
not removed.
--------------------------------------------------------------------
*/
#include<immintrin.h>
#include "segseg.h"
#include "rtreeindex_512.h" 
#include<stdlib.h>
#include<iostream>
using namespace std;
/*---------------------------------------------------------------------
SegSegInt: Finds the point of intersection p between two closed
segments ab and cd.  Returns p and a char with the following meaning:
   'e': The segments collinearly overlap, sharing a point.
   'v': An endpoint (vertex) of one segment is on the other segment,
        but 'e' doesn't hold.
   '1': The segments intersect properly (i.e., they share a point and
        neither 'v' nor 'e' holds).
   '0': The segments do not intersect (i.e., they share no points).
Note that two collinear segments that share just one point, an endpoint
of each, returns 'e' rather than 'v' as one might expect.
---------------------------------------------------------------------*/
/*
double* dot(__m256d A, __m256d B)
{
    // allocating memory to store 2 doubles for 32 byte (256 bits) aligned memory
    double *ans = (double*)_mm_malloc( 2 * sizeof(double), 32);
    __m256d vect = _mm256_mul_pd(A, B);
    //cout<<"A "<<A << "  B "<< B << vect <<endl;
    
    double temp[4];
    _mm256_store_pd(temp, vect);
    
    ans[0] = temp[0] + temp[1];
    ans[1] = temp[2] + temp[3];
    // ans[0] = vect[0]+vect[1];
    // ans[1] = vect[2]+vect[3];
    return ans;
}

char* findCrossingLines(struct Point p1, struct Point p2, 
                                vector<struct Point> *p3, vector<struct Point> *p4)
{
    char *result = (char *)malloc( 2 * sizeof( char ) );
    
    //printf("Beginning of findCrossingLines\n");

    __m256d PNT1 = _mm256_set_pd(p1.y,p1.x,p1.y,p1.x);
    __m256d PNT2 = _mm256_set_pd(p2.y,p2.x,p2.y,p2.x);   
     
    __m256d a = _mm256_sub_pd(PNT2, PNT1);

    //printf("%f %f %f %f \n ", p3[0], p3[1], p3[2], p3[3]);
    //cout<<p3[0]<< " "<< p3[1]<< " "<< p3[2]<< " "<< p3[3]<<endl;
    //printf("%f %f %f %f \n ", p3->at(1).y, p3->at(1).x, p3->at(0).y, p3->at(0).x );
    //printf("%f %f %f %f \n ", p4->at(1).y, p4->at(1).x, p4->at(0).y, p4->at(0).x );
     
    __m256d PNTS3 = _mm256_set_pd(p3->at(1).y, p3->at(1).x, p3->at(0).y, p3->at(0).x);
    __m256d PNTS4 = _mm256_set_pd(p4->at(1).y, p4->at(1).x, p4->at(0).y, p4->at(0).x);
    __m256d b = _mm256_sub_pd(PNTS4, PNTS3);

    __m256d c = _mm256_sub_pd(PNTS4, PNT2);
    
    __m256d b_perp = _mm256_set_pd(b[2],-b[3],b[0],-b[1]);

    //printf("Before Vectors Created\n");

    double* numerator = dot(b_perp,c);
    double* denominator = dot(b_perp,a); 

    //printf("After Vectors Created\n");

    for (int i=0; i<2; i++)
    {
    	bool isParallel = (denominator[i] == 0);
	    if(!isParallel)
	    {
                // Not Parallel

		  double quotient = numerator[i] / denominator[i];
                  Point intersectionPoint(quotient * a[X] + PNT2[X], quotient * a[Y] + PNT2[Y]);

                //printf("Intersection Point found\n");
             
		  __m256d arr1 = _mm256_set_pd(p1.y, p1.x, p3->at(i).y, p3->at(i).x);
	          __m256d arr2 = _mm256_set_pd(p2.y, p2.x, p4->at(i).y, p4->at(i).x);
                

           if(lessThanIntersectionPoint(arr1, arr2, intersectionPoint.x)
                       && greaterThanIntersectionPoint(arr1, arr2, intersectionPoint.y))
 		   {     
			 result[i] = '0';
		   }
		   else
		   {
 			result[i] = '1';
		   }

	}
	else
	{
		result[i] = '1';
	}    
    }
    return result;
}


// compares segment end points to the calculated intersection point
bool lessThanIntersectionPoint(__m256d vec1, __m256d vec2, double intersection){
        __m256d siX = _mm256_set1_pd(intersection);

	__m256d minVec = _mm256_min_pd(vec1, vec2);

        __m256d ans = _mm256_cmp_pd(siX, minVec, _CMP_GE_OQ);	
        
        if(ans[0] && ans[1] && ans[2] && ans[3])
        {
		return true;
	}
        else
        {
		return false;
	}

}

// compares segment end points to the calculated intersection point
bool greaterThanIntersectionPoint(__m256d vec1, __m256d vec2, double intersection){
        
        __m256d siY = _mm256_set1_pd(intersection);

	__m256d maxVec = _mm256_max_pd(vec1, vec2);
        
        __m256d ans = _mm256_cmp_pd(siY, maxVec, _CMP_LE_OQ); 

        if(ans[0] && ans[1] && ans[2] && ans[3])
        {
		return true;
	}
 	else
	{
		return false;
	}
}

*/
int dot(tPointi a, tPointi b){
     int add1 = a[X] + a[Y];
     int add2 = b[X] + b[Y];
     int ans = add1*add2;
     return ans;
}

bool SegSeg(tPointi p1, tPointi p2, 
                                tPointi p3, tPointi p4)
{
    bool result; 
    tPointi A = {p2[0], p2[1]};
    tPointi a = {p2[0] - p1[0], p2[1] - p1[1]};
    tPointi B = {p4[0], p4[1]};
    tPointi b = {p4[0] - p3[0], p4[1] - p3[1]};
    tPointi c = {B[0] - A[0], B[1]-A[1]};
    tPointi b_perp = {-b[1], b[0]};
    
    int numerator = dot(b_perp, c);
    int denominator = dot(b_perp, a);
    bool isParallel = (denominator == 0);


    float quotient = numerator / denominator;
    tPointi intersectionPoint = {(int)quotient * a[0] + A[0], (int)quotient * a[1] + A[1]};

    result = (!isParallel && 
                  intersectionPoint[0] >= min(p1[0], p2[0]) && 
                  intersectionPoint[0] >= min(p3[0], p4[0]) &&
                  intersectionPoint[0] <= max(p1[0], p2[0]) && 
                  intersectionPoint[0] <= max(p3[0], p4[0]) &&
                  intersectionPoint[1] >= min(p1[1], p2[1]) && 
                  intersectionPoint[1] >= min(p3[1], p4[1]) &&
                  intersectionPoint[1] <= max(p1[1], p2[1]) && 
                  intersectionPoint[1] <= max(p3[1], p4[1]));
 
    return result;
}


char SegSegInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p )
{
   //printf("a, b: (%d, %d), (%d, %d) - c, d: (%d, %d), (%d, %d)\n\n", a[0],a[1], b[0], b[1], c[0],c[1], d[0], d[1]);

   double  s, t;       /* The two parameters of the parametric eqns. */
   double num, denom;  /* Numerator and denoninator of equations. */
   char code = '?';    /* Return char characterizing intersection. */

   denom = a[X] * (double)( d[Y] - c[Y] ) +
           b[X] * (double)( c[Y] - d[Y] ) +
           d[X] * (double)( b[Y] - a[Y] ) +
           c[X] * (double)( a[Y] - b[Y] );

   /* If denom is zero, then segments are parallel: handle separately. */
   if (denom == 0.0)
      return  ParallelInt(a, b, c, d, p);

   num =    a[X] * (double)( d[Y] - c[Y] ) +
            c[X] * (double)( a[Y] - d[Y] ) +
            d[X] * (double)( c[Y] - a[Y] );
   if ( (num == 0.0) || (num == denom) ) code = 'v';
   s = num / denom;
   //printf("num=%lf, denom=%lf, s=%lf\n", num, denom, s);

   num = -( a[X] * (double)( c[Y] - b[Y] ) +
            b[X] * (double)( a[Y] - c[Y] ) +
            c[X] * (double)( b[Y] - a[Y] ) );
   if ( (num == 0.0) || (num == denom) ) code = 'v';
   t = num / denom;
   //printf("num=%lf, denom=%lf, t=%lf\n", num, denom, t);

   if      ( (0.0 < s) && (s < 1.0) &&
             (0.0 < t) && (t < 1.0) )
     code = '1';
   else if ( (0.0 > s) || (s > 1.0) ||
             (0.0 > t) || (t > 1.0) )
     code = '0';

   p[X] = a[X] + s * ( b[X] - a[X] );
   p[Y] = a[Y] + s * ( b[Y] - a[Y] );

   return code;
}

char	ParallelInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p )
{
   if ( !Collinear( a, b, c) )
      return '0';

   if ( Between( a, b, c ) ) {
      Assigndi( p, c );
      return 'e';
   }
   if ( Between( a, b, d ) ) {
      Assigndi( p, d );
      return 'e';
   }
   if ( Between( c, d, a ) ) {
      Assigndi( p, a );
      return 'e';
   }
   if ( Between( c, d, b ) ) {
      Assigndi( p, b );
      return 'e';
   }
   return '0';
}
void	Assigndi( tPointd p, tPointi a )
{
   int i;
   for ( i = 0; i < DIM; i++ )
      p[i] = a[i];
}
/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
Assumes it is already known that abc are collinear.
---------------------------------------------------------------------*/
bool    Between( tPointi a, tPointi b, tPointi c )
{
   tPointi      ba, ca;

   /* If ab not vertical, check betweenness on x; else on y. */
   if ( a[X] != b[X] )
      return ((a[X] <= c[X]) && (c[X] <= b[X])) ||
             ((a[X] >= c[X]) && (c[X] >= b[X]));
   else
      return ((a[Y] <= c[Y]) && (c[Y] <= b[Y])) ||
             ((a[Y] >= c[Y]) && (c[Y] >= b[Y]));
}

int Collinear( tPointi a, tPointi b, tPointi c )
{
   return AreaSign( a, b, c ) == 0;
}
int     AreaSign( tPointi a, tPointi b, tPointi c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
            ( c[0] - a[0] ) * (double)( b[1] - a[1] );

    /* The area should be an integer. */
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
}

// int main()
// {
//    return 0;
// }