//this version is used to optimize with float data types

#include "rtreeindex2_512.h"
#include "segseg.h"
#include <stdio.h>
#include <math.h>
//uncomment these for utilizing threads
//#include <thread>
//#include <future>

int count;

float max(float num1, float num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

float min(float num1, float num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

/* J.G. In what way do we use AVX512 to optimize this code? Adding more placeholders seems redundant, could we compare more than two MBRs? */ 

int* overlap(struct MBR *R1, struct MBR *R2, struct MBR *R3, struct MBR *R4, struct MBR *R5, struct MBR *R6, struct MBR *R7, struct MBR *R8)
{
    //0: xmin, 1: ymin, 2: xmax, 3: ymax

    /* single precision */

    float placeholder = 0;

    int *n = new int[4];

    __m512 R1MBR1 = _mm512_set_ps(R7->boundary->at(3), R7->boundary->at(2), R7->boundary->at(1), R7->boundary->at(0), R5->boundary->at(3), R5->boundary->at(2), R5->boundary->at(1), R5->boundary->at(0), R3->boundary->at(3), R3->boundary->at(2), R3->boundary->at(1), R3->boundary->at(0), R1->boundary->at(3), R1->boundary->at(2), R1->boundary->at(1), R1->boundary->at(0));
    __m512 R2MBR1 = _mm512_set_ps(R8->boundary->at(1), R8->boundary->at(0), R8->boundary->at(3), R8->boundary->at(2), R6->boundary->at(1), R6->boundary->at(0), R6->boundary->at(3), R6->boundary->at(2), R4->boundary->at(1), R4->boundary->at(0), R4->boundary->at(3), R4->boundary->at(2), R2->boundary->at(1), R2->boundary->at(0), R2->boundary->at(3), R2->boundary->at(2));
    __mmask16 comparison = _mm512_cmp_ps_mask(R1MBR1, R2MBR1, 5); //greater than comparison
    __mmask16 mask = 0000000000000001;
    
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
    if (comparison >> 8 & mask || 
	comparison >> 9 & mask ||
	!(comparison >> 10 & mask) ||
	!(comparison >> 11 & mask))
    {
        //printf("Mask Based Result: 0\n");
        n[2] = 0;
    }
    else
    {
	//printf("Mask Based Result: 1\n");
        n[2] = 1;
    }

    if (comparison >> 12 & mask ||
	comparison >> 13 & mask ||
	!(comparison >> 14 & mask) ||
	!(comparison >> 15 & mask))
    {
        //printf("Mask Based Result: 0\n");
        n[3] = 0;
    }
    else
    {
	//printf("Mask Based Result: 1\n");
        n[3] = 1;
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
        for(j = 0; j < FANOUT; j += 4)
        {
	    n = overlap(R1->children->at(i)->mbr, R2->children->at(j)->mbr, R1->children->at(i)->mbr, R2->children->at(j+1)->mbr, R1->children->at(i)->mbr, R2->children->at(j+2)->mbr, R1->children->at(i)->mbr, R2->children->at(j+3)->mbr);
            
            if(n[0]==1)
            {
                overlappingPairs.push_back(make_pair(R1->children->at(i), R2->children->at(j)));
            }
            if(n[1]==1)
            {
                overlappingPairs.push_back(make_pair(R1->children->at(i), R2->children->at(j+1)));
	    }
            if(n[2]==1)
            {
                overlappingPairs.push_back(make_pair(R1->children->at(i), R2->children->at(j+2)));
            }
            if(n[3]==1)
            {
                overlappingPairs.push_back(make_pair(R1->children->at(i), R2->children->at(j+3)));
	    }
        }
    }
    return overlappingPairs;
}

/*
void rectUnion(vector<RectReal> *rect1, vector<RectReal> *rect2)//parallize?
{
    __m512 r1 = _mm512_set_ps(0,0,0,0,0,0,0,0,rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3));
    __m512 r2 = _mm512_set_ps(rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3));
    __m256d comparison = _mm256_cmp_pd(r1, r2, 14); //greater than comparison
    
    rect1->at(0) = (comparison[3] == 0) ?  rect1->at(0) : rect2->at(0);  // if ( r1_minx is LTE r2_minx )
    rect1->at(1) = (comparison[2] == 0) ?  rect1->at(1) : rect2->at(1);  // if ( r1_miny is LTE r2_miny )
    
    // isnan => Is Not-A-Number
    // Returns whether x is a NaN (Not-A-Number) value.
    // isnan ( x ) argument is a floating point number 
    // Returns a non-zero value (true) if x is a NaN value; and zero (false) otherwise.
    rect1->at(2) = (isnan(comparison[1])) ?  rect1->at(2) : rect2->at(2);
    rect1->at(3) = (isnan(comparison[0])) ?  rect1->at(3) : rect2->at(3);

}*/

void rectUnion(vector<RectReal> *rect1, vector<RectReal> *rect2)//parallize?
{
    rect1->at(0) = min((double)rect1->at(0), (double)rect2->at(0));
    rect1->at(1) = min((double)rect1->at(1), (double)rect2->at(1));
    rect1->at(2) = max((double)rect1->at(2), (double)rect2->at(2));
    rect1->at(3) = max((double)rect1->at(3), (double)rect2->at(3));
}


void rectUnion4(vector<RectReal> *rect1, vector<RectReal> *rect2, vector<RectReal> *rect3, vector<RectReal> *rect4)
{
    __m256 r1 = _mm256_set_ps(rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3), rect3->at(0), rect3->at(1), rect3->at(2), rect3->at(3));
    __m256 r2 = _mm256_set_ps(rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3), rect4->at(0), rect4->at(1), rect4->at(2), rect4->at(3));
    __m256 comparison = _mm256_cmp_ps(r1, r2, 14); //greater than comparison
    
    rect1->at(0) = (comparison[7] == 0) ?  rect1->at(0) : rect2->at(0);
    rect1->at(1) = (comparison[6] == 0) ?  rect1->at(1) : rect2->at(1);
    rect1->at(2) = (isnan(comparison[5])) ?  rect1->at(2) : rect2->at(2);
    rect1->at(3) = (isnan(comparison[4])) ?  rect1->at(3) : rect2->at(3);

    rect3->at(0) = (comparison[3] == 0) ?  rect3->at(0) : rect4->at(0);
    rect3->at(1) = (comparison[2] == 0) ?  rect3->at(1) : rect4->at(1);
    rect3->at(2) = (isnan(comparison[1])) ?  rect3->at(2) : rect4->at(2);
    rect3->at(3) = (isnan(comparison[0])) ?  rect3->at(3) : rect4->at(3);

    rectUnion(rect1, rect3);
}

/* J.G. New method rectUnion8 using AVX512: we can do union of 8 children mbrs instead of 4 */
void rectUnion8(vector<RectReal> *rect1, vector<RectReal> *rect2, vector<RectReal> *rect3, vector<RectReal> *rect4, vector<RectReal> *rect5, vector<RectReal> *rect6, vector<RectReal> *rect7, vector<RectReal> *rect8)
{
    __m512 r1 = _mm512_set_ps(rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3), rect3->at(0), rect3->at(1), rect3->at(2), rect3->at(3), rect5->at(0), rect5->at(1), rect5->at(2), rect5->at(3), rect7->at(0), rect7->at(1), rect7->at(2), rect7->at(3));
    __m512 r2 = _mm512_set_ps(rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3), rect4->at(0), rect4->at(1), rect4->at(2), rect4->at(3), rect6->at(0), rect6->at(1), rect6->at(2), rect6->at(3), rect8->at(0), rect8->at(1), rect8->at(2), rect8->at(3));

    __mmask16 mask = 0000000000000001;
    __mmask16 comparison = _mm512_cmp_ps_mask(r1, r2, 5); //greater than comparison
    
    rect1->at(0) = (comparison >> 15 ^ mask) ?  rect1->at(0) : rect2->at(0);
    rect1->at(1) = (comparison >> 14 ^ mask) ?  rect1->at(1) : rect2->at(1);
    rect1->at(2) = (comparison >> 13 & mask) ?  rect1->at(2) : rect2->at(2);
    rect1->at(3) = (comparison >> 12 & mask) ?  rect1->at(3) : rect2->at(3); 

    rect3->at(0) = (comparison >> 11 ^ mask) ?  rect3->at(0) : rect4->at(0);
    rect3->at(1) = (comparison >> 10 ^ mask) ?  rect3->at(1) : rect4->at(1);
    rect3->at(2) = (comparison >> 9 & mask) ?  rect3->at(2) : rect4->at(2);
    rect3->at(3) = (comparison >> 8 & mask) ?  rect3->at(3) : rect4->at(3);

    rect5->at(0) = (comparison >> 7 ^ mask) ?  rect5->at(0) : rect6->at(0);
    rect5->at(1) = (comparison >> 6 ^ mask) ?  rect5->at(1) : rect6->at(1);
    rect5->at(2) = (comparison >> 5 & mask) ?  rect5->at(2) : rect6->at(2);
    rect5->at(3) = (comparison >> 4 & mask) ?  rect5->at(3) : rect6->at(3);

    rect7->at(0) = (comparison >> 3 ^ mask) ?  rect7->at(0) : rect8->at(0);
    rect7->at(1) = (comparison >> 2 ^ mask) ?  rect7->at(1) : rect8->at(1);
    rect7->at(2) = (comparison >> 1 & mask) ?  rect7->at(2) : rect8->at(2);
    rect7->at(3) = (comparison >> 0 & mask) ?  rect7->at(3) : rect8->at(3);

    rectUnion4(rect1, rect3, rect5, rect7);
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
    float minX = p.x;
    float minY = p.y;

    float maxX = minX;
    float maxY = minY;

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
    int rect, x;
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
    
    //rectUnion8(unionRect, nodes->at(1)->mbr->boundary, tempNode2, nodes->at(3)->mbr->boundary);
    for(rect = 1; rect < n; rect += 1)
    {
        rectUnion(unionRect, (nodes->at(rect)->mbr->boundary));
    }

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

        //uncomment these for utilizing threads
        // int range0[2], range1[2];
        // getRange(range0, 0, low, high);
        // getRange(range1, 1, low, high);
        // future<struct Node *> future1 = async(launch::async, createRTree, ptArr, range0[0], range0[1]);

        // struct Node *child2 = createRTree(ptArr, range1[0], range1[1]);
        // struct Node *child1 = future1.get();

        // tempChildren->push_back(child1);
        // tempChildren->push_back(child2);

        //comment out this for loop for utilizing threads
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

	        __m256d r1 = _mm256_set_pd(g[1],g[0], e[1], e[0]);
    	        __m256d r2 = _mm256_set_pd(h[1],h[0], f[1], f[0]);
                __m256d minvec = _mm256_min_pd(r1, r2);
                __m256d maxvec = _mm256_max_pd(r1, r2);

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
			//printf("\n%d", count);
                    	//printf("p((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i+1).x, ptArrP->at(i+1).y,
                    	//                                                       ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j+1).x, ptArrQ->at(j+1).y);
                	}
                }
                
            }
        }
    }
    else if(R1->leaf == 1)
    {  
        for(childC = 0; childC < FANOUT; childC += 4)
        {

            // does not work for fanout < 4
            struct Node *C1 = R2->children->at(childC);
            struct Node *C2 = R2->children->at(childC+1);
	    struct Node *C3 = R2->children->at(childC+2);
            struct Node *C4 = R2->children->at(childC+3);

            int* arr = overlap(R1->mbr, C1->mbr, R1->mbr, C2->mbr, R1->mbr, C3->mbr, R1->mbr, C4->mbr);

            if(arr[0])
            {
                spatialJoin(R1, C1, ptArrP, ptArrQ);
            }
            if(arr[1]){
  		spatialJoin(R1, C2, ptArrP, ptArrQ);
	    }
	    if(arr[3]){
  		spatialJoin(R1, C3, ptArrP, ptArrQ);
	    }
            if(arr[4]){
  		spatialJoin(R1, C4, ptArrP, ptArrQ);
	    }

        }
    }
    else if(R2->leaf == 1)
    {
        for(childS = 0; childS < FANOUT; childS += 4)
        {
            struct Node *C1 = R1->children->at(childS);
            struct Node *C2 = R1->children->at(childS+1);
            struct Node *C3 = R1->children->at(childS+2);
            struct Node *C4 = R1->children->at(childS+3);
    
            int* arr = overlap(C1->mbr, R2->mbr, C2->mbr, R2->mbr, C3->mbr, R2->mbr, C4->mbr, R2->mbr);

            if(arr[0])
            {
                spatialJoin(C1, R2, ptArrP, ptArrQ);
            }
            if(arr[1]){
  		spatialJoin(C2, R2, ptArrP, ptArrQ);
	    }
            if(arr[2])
            {
                spatialJoin(C3, R2, ptArrP, ptArrQ);
            }
            if(arr[3]){
  		spatialJoin(C4, R2, ptArrP, ptArrQ);
	    }

        }
    }
    else
    {
        int i;
        vector<pair<struct Node *, struct Node *> > candidatePairs = overlappingChildren(R1, R2);
        for(i = 0; i < candidatePairs.size(); ++i) 
        {
            spatialJoin(candidatePairs.at(i).first, candidatePairs.at(i).second, ptArrP, ptArrQ);
        }

    }
        
}

//edge to edge comparison to verify correct results
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

