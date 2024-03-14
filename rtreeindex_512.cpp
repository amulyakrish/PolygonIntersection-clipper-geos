#include "rtreeindex_512.h"
#include "segseg.h"
#include <vector>
#include <stdio.h>
#include <math.h>
#include <thread>
#include <future>

using namespace std;
using std::vector;

int count;

double max(double num1, double num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

double min(double num1, double num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

int* overlap(struct MBR *R1, struct MBR *R2, struct MBR *R3, struct MBR *R4)
{
    //0: xmin, 1: ymin, 2: xmax, 3: ymax

    double placeholder = 0;
    int* n = new int[2];

    __m512d R1MBR1 = _mm512_set_pd(R3->boundary->at(3), R3->boundary->at(2), R3->boundary->at(1), R3->boundary->at(0), R1->boundary->at(3), R1->boundary->at(2), R1->boundary->at(1), R1->boundary->at(0));
    __m512d R2MBR1 = _mm512_set_pd(R4->boundary->at(1), R4->boundary->at(0), R4->boundary->at(3), R4->boundary->at(2), R2->boundary->at(1), R2->boundary->at(0), R2->boundary->at(3), R2->boundary->at(2));
    __mmask8 comparison = _mm512_cmp_pd_mask(R1MBR1, R2MBR1, _MM_CMPINT_GT); //greater than comparison
    __mmask8 mask = 00000001;

    // [TODO] : Return array for two comparison results

    if (comparison & mask || 
	comparison >> 1 & mask ||
	!(comparison >> 2 & mask) ||
	!(comparison >> 3 & mask))
    {
        //printf("Mask Based Result: 0\n");
        n[0] = 0;
    }
    else
    {
	//printf("Mask Based Result: 1\n");
        n[0] = 1;
    }

    if (comparison >> 4 & mask ||
	comparison >> 5 & mask ||
	!(comparison >> 6 & mask) ||
	!(comparison >> 7 & mask))
    {
        //printf("Mask Based Result: 0\n");
        n[1] = 0;
    }
    else
    {
	//printf("Mask Based Result: 1\n");
        n[1] = 1;
    }

    return n;
}

vector<pair<struct Node *, struct Node *> > overlappingChildren(struct Node *R1, struct Node *R2) // TEMPORARY CHANGE PLEASE CHANGE LATER - JG
{
    int i, j;
    int* n;
    vector<pair<struct Node *, struct Node *> > overlappingPairs;
    for(i = 0; i < FANOUT; i += 1)
    {
        for(j = 0; j < FANOUT; j += 2)
        {
	    n = overlap(R1->children->at(i)->mbr, R2->children->at(j)->mbr, R1->children->at(i)->mbr, R2->children->at(j+1)->mbr);
            
            if(n[0]==1)
            {
                overlappingPairs.push_back(make_pair(R1->children->at(i), R2->children->at(j)));
            }
            if(n[1]==1)
            {
                overlappingPairs.push_back(make_pair(R1->children->at(i), R2->children->at(j+1)));
	    }
        }
    }
    return overlappingPairs;
}
/*
void rectUnion(vector<RectReal> *rect1, vector<RectReal> *rect2)//parallize?
{
    __m256d r1 = _mm256_set_pd(rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3));
    __m256d r2 = _mm256_set_pd(rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3));
    __m256d comparison = _mm256_cmp_pd(r1, r2, 14); //greater than comparison
    
    rect1->at(0) = (comparison[3] == 0) ?  rect1->at(0) : rect2->at(0);
    rect1->at(1) = (comparison[2] == 0) ?  rect1->at(1) : rect2->at(1);
    rect1->at(2) = (isnan(comparison[1])) ?  rect1->at(2) : rect2->at(2);
    rect1->at(3) = (isnan(comparison[0])) ?  rect1->at(3) : rect2->at(3);

}
*/
void rectUnion(vector<RectReal> *rect1, vector<RectReal> *rect2)//parallize?
{
    __mmask8 mask = 00000001;
    
    __m512d r1 = _mm512_set_pd(0,0,0,0,rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3));
    __m512d r2 = _mm512_set_pd(0,0,0,0,rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3));
    
    __mmask8 comparison = _mm512_cmp_pd_mask(r1, r2, 5); //greater than comparison

    rect1->at(0) = !(comparison >> 3 & mask) ?  rect1->at(0) : rect2->at(0);
    rect1->at(1) = !(comparison >> 2 & mask) ?  rect1->at(1) : rect2->at(1);
    rect1->at(2) = (comparison >> 1 & mask) ?  rect1->at(2) : rect2->at(2);
    rect1->at(3) = (comparison & mask) ?  rect1->at(3) : rect2->at(3);

}

/*
void rectUnion4(vector<RectReal> *rect1, vector<RectReal> *rect2, vector<RectReal> *rect3, vector<RectReal> *rect4)//parallize?
{
    __m512d r1 = _mm512_set_pd(0,0,0,0,rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3));
    __m512d r2 = _mm512_set_pd(0,0,0,0,rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3));
    __m512d minvec = _mm512_min_pd(r1, r2);
    __m512d maxvec = _mm512_max_pd(r1, r2);

    rect1->at(0) = minvec[3];	rect1->at(1) = minvec[2];
    rect1->at(2) = maxvec[1];	rect1->at(3) = maxvec[0];

}


*/
/* J.G. New method rectUnion4 using AVX512*/

void rectUnion4(vector<RectReal> *rect1, vector<RectReal> *rect2, vector<RectReal> *rect3, vector<RectReal> *rect4)
{

    __mmask8 mask = 00000001;

    __m512d r1 = _mm512_set_pd(rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3), rect3->at(0), rect3->at(1), rect3->at(2), rect3->at(3));
    __m512d r2 = _mm512_set_pd(rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3), rect4->at(0), rect4->at(1), rect4->at(2), rect4->at(3));

    __mmask8 comparison = _mm512_cmp_pd_mask(r1, r2, 5); //greater than comparison
    
    rect1->at(0) = !(comparison >> 7 & mask) ?  rect1->at(0) : rect2->at(0);
    rect1->at(1) = !(comparison >> 6 & mask) ?  rect1->at(1) : rect2->at(1);
    rect1->at(2) = (comparison >> 5 & mask) ?  rect1->at(2) : rect2->at(2);
    rect1->at(3) = (comparison >> 4 & mask) ?  rect1->at(3) : rect2->at(3);

    rect3->at(0) = !(comparison >> 3 & mask) ?  rect3->at(0) : rect4->at(0);
    rect3->at(1) = !(comparison >> 2 & mask) ?  rect3->at(1) : rect4->at(1);
    rect3->at(2) = (comparison >> 1 & mask) ?  rect3->at(2) : rect4->at(2);
    rect3->at(3) = (comparison & mask) ?  rect3->at(3) : rect4->at(3);

    rectUnion(rect1, rect3);
}

void getRange(int *range, int id, int low, int high)// does this not last indext back to first index? for polygons
{
    int size = high - low;
    int start = low + (id*size)/FANOUT;
    range[0] = start;

    int end;
    if(id == (FANOUT-1))
    {
        end = high;
    }
    else
    {
        end = low + ((id+1)*size)/FANOUT;
    }
    range[1] = end;
}

struct MBR * createTile(vector<Point> *ptArr, int first, int last)
{

    int index;
    Point p = ptArr->at(first);
    double minX = p.x;
    double minY = p.y;

    double maxX = minX;
    double maxY = minY;

    for(index = first; index <= last; index += 1)
    {
        struct Point pt = ptArr->at(index);
        minX = min(pt.x, minX);
        minY = min(pt.y, minY);
        maxX = max(pt.x, maxX);
        maxY = max(pt.y, maxY);
    }

    struct MBR *t = (struct MBR *)malloc(sizeof(struct MBR));
    t->start = first;
    t->end = last;
    
    
    t->boundary = new vector<RectReal>[4];     
    t->boundary->push_back(minX);
    t->boundary->push_back(minY);
    t->boundary->push_back(maxX);
    t->boundary->push_back(maxY);
    return t; 
}

struct Node * createLeaf(vector<Point> *ptArr, int low, int high)
{
    if((low == high) || (low > high))
        return NULL;
    struct Node *leaf = (struct Node *)malloc(sizeof(struct Node));
    struct MBR *env = createTile(ptArr, low, high);
    leaf->mbr = env;
    leaf->leaf = 1;
    return leaf;
}

struct MBR * unionJoin(vector<struct Node *> *nodes, int n)
{
    struct MBR *newMBR = (struct MBR *)malloc(sizeof(struct MBR));
    vector<RectReal> *unionRect = new vector<RectReal>[NUMSIDES]; 
    vector<RectReal> *tempNode2 = new vector<RectReal>[NUMSIDES];
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(0));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(1));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(2));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(3));

    tempNode2->push_back(nodes->at(2)->mbr->boundary->at(0));
    tempNode2->push_back(nodes->at(2)->mbr->boundary->at(1));
    tempNode2->push_back(nodes->at(2)->mbr->boundary->at(2));
    tempNode2->push_back(nodes->at(2)->mbr->boundary->at(3));
    
    rectUnion4(unionRect, nodes->at(1)->mbr->boundary, tempNode2, nodes->at(3)->mbr->boundary);

    newMBR->start = nodes->at(0)->mbr->start;
    newMBR->end = nodes->at(n-1)->mbr->end;
 
    newMBR->boundary = unionRect;

    return newMBR;
}

struct Node * createRTree(vector<struct Point> *ptArr,int low, int high)
{
    struct Node *root;
    if((high - low) <= BUNDLEFACTOR)
    {
        root = createLeaf(ptArr, low, high);
    }
    else
    {
        int childID;
        root = (struct Node *)malloc(sizeof(struct Node));
        vector<struct Node *> *tempChildren = new vector<struct Node *>[FANOUT];


        //uncomment for threads
        // int range[2], range1[2];
        // getRange(range0, 0, low, high);
        // getRange(range1, 1, low, high);
        // future<struct Node *> future1 = async(launch::async, createRTree, ptArr, range0[0], range0[1]);

        // struct Node *child2 = createRTree(ptArr, range1[0], range1[1]);
        // struct Node *child1 = future1.get();

        // tempChildren->push_back(child1);
        // tempChildren->push_back(child2);

        //comment out for loop for threads
        for(childID = 0; childID < FANOUT; childID += 1)
        {
            int range[2];
            getRange(range, childID, low, high);
            struct Node *child = createRTree(ptArr, range[0], range[1]);
            tempChildren->push_back(child);
        }
        root->children = tempChildren;
        root->leaf = 0;
        struct MBR *parentTile = unionJoin(root->children, FANOUT); //test unionJoin
        root->mbr = parentTile;
    }
    return root;
}

void spatialJoin(struct Node *R1, struct Node *R2, vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ)
{
    int childS, childC;
    if((R1->leaf == 1) and (R2->leaf == 1)) //if C1 and C2 are leaf Nodes
    {
        //printf("R1 Leaf R2 Leaf\n");
        int i, j;
        for(i = R1->mbr->start; i < (R1->mbr->end); i += 1)
        {
            for(j = R2->mbr->start; j < (R2->mbr->end); j += 1) //possible bug in these for loops hitting out of bounds
            {
                tPointd p;
                char ans;
		int ifOverlap;

		tPointd e = {ptArrP->at(i).x, ptArrP->at(i).y};
                tPointd f = {ptArrP->at(i+1).x, ptArrP->at(i+1).y};
                tPointd g = {ptArrQ->at(j).x, ptArrQ->at(j).y};
                tPointd h = {ptArrQ->at(j+1).x, ptArrQ->at(j+1).y};            

	        __m512d r1 = _mm512_set_pd(0,0,0,0, g[1],g[0], e[1], e[0]);
    	        __m512d r2 = _mm512_set_pd(0,0,0,0, h[1],h[0], f[1], f[0]);
                __m512d minvec = _mm512_min_pd(r1, r2);
                __m512d maxvec = _mm512_max_pd(r1, r2);

                double xmin = minvec[2];
    	        double ymin = minvec[3];
    	        double xmax = maxvec[2];
       	        double ymax = maxvec[3];

                double xmin2 = minvec[0];
    	    	double ymin2 = minvec[1];
    	        double xmax2 = maxvec[0];
       	        double ymax2 = maxvec[1];

	        if( xmin > xmax2 || xmin2 > xmax || ymax < ymin2 || ymax2 < ymin)
	            	ifOverlap = 0;
           	else
 		 	ifOverlap = 1;

		if(ifOverlap == 1)
		{
                	tPointi a = {(int)ptArrP->at(i).x, (int)ptArrP->at(i).y};
                	tPointi b = {(int)ptArrP->at(i+1).x, (int)ptArrP->at(i+1).y};
                	tPointi c = {(int)ptArrQ->at(j).x, (int)ptArrQ->at(j).y};
                	tPointi d = {(int)ptArrQ->at(j+1).x, (int)ptArrQ->at(j+1).y};         

	                ans = SegSegInt(a, b, c, d, p);
        	        if (ans == '1')
                	{
                    	//intersection fonud
			//count++;
			//printf("\n%d",count);
                    	//printf("p((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i+1).x, ptArrP->at(i+1).y,
                      	  //                                                   ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j+1).x, ptArrQ->at(j+1).y);
                	}
		}
            }
        }
    }
    else if(R1->leaf == 1)
    {  
        //printf("R1 Leaf R2 Not Leaf\n");        
        for(childC = 0; childC < FANOUT; childC += 2)
        {
            struct Node *C1 = R2->children->at(childC);
            struct Node *C2 = R2->children->at(childC+1);

            int* arr = overlap(R1->mbr, C1->mbr, R1->mbr, C2->mbr);

            if(arr[0])
            {
                spatialJoin(R1, C1, ptArrP, ptArrQ);
            }
            if(arr[1]){
  	       spatialJoin(R1, C2, ptArrP, ptArrQ);
	   }
        }
    }
    else if(R2->leaf == 1)
    {

        //printf("R1 Not Leaf R2 Leaf\n");
        for(childS = 0; childS < FANOUT; childS += 2)
        {
            struct Node *C1 = R1->children->at(childS);
            struct Node *C2 = R1->children->at(childS+1);
    
            int* arr = overlap(C1->mbr, R2->mbr, C2->mbr, R2->mbr);

            if(arr[0])
            {
                spatialJoin(C1, R2, ptArrP, ptArrQ);
            }
            if(arr[1]){
  		spatialJoin(C2, R2, ptArrP, ptArrQ);
	    }
        }

    }
    else
    {
        //printf("R1 Not Leaf R2 Not Leaf\n");
        int i;
        vector<pair<struct Node *, struct Node *> > candidatePairs = overlappingChildren(R1, R2);
        for(i = 0; i < candidatePairs.size(); ++i) 
        { 
            //printf("loop %d", i); 
            spatialJoin(candidatePairs.at(i).first, candidatePairs.at(i).second, ptArrP, ptArrQ);
        }

    }
        
}

void checkLSI(vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ)
{
    printf("\n\nBrute force LSI\n");
    int i, j;
    for(i = 0; i < ptArrP->size()-1; i += 1)
    {
        for(j = 0; j < ptArrQ->size()-1; j += 1)
        {
            tPointd p;
            char ans;
            tPointi a = {(int)ptArrP->at(i).x, (int)ptArrP->at(i).y};
            tPointi b = {(int)ptArrP->at(i+1).x, (int)ptArrP->at(i+1).y};
            tPointi c = {(int)ptArrQ->at(j).x, (int)ptArrQ->at(j).y};
            tPointi d = {(int)ptArrQ->at(j+1).x, (int)ptArrQ->at(j+1).y};         

            ans = SegSegInt(a, b, c, d, p);
            if (ans == '1')
            {
                printf("p((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i+1).x, ptArrP->at(i+1).y,
                                                                            ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j+1).x, ptArrQ->at(j+1).y);
            } 
        }
    }
}


void printNode(struct Node *node)
{
    struct MBR *tile = node->mbr;
    printf("%d : %d tile coordinates: (%f, %f)-(%f, %f)\n",tile->start, tile->end, tile->boundary->at(0), tile->boundary->at(1), 
                                                         tile->boundary->at(2), tile->boundary->at(3));
}

