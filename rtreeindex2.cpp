//this version is used to optimize with float data types

#include "rtreeindex2.h"
#include "segseg.h"
#include <stdio.h>
#include <math.h>
//uncomment these for utilizing threads
//#include <thread>
//#include <future>

int count = 0;

float max(float num1, float num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

float min(float num1, float num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

int* overlap(struct MBR *R1, struct MBR *R2, struct MBR *R3, struct MBR *R4)
{
    //0: xmin, 1: ymin, 2: xmax, 3: ymax
    int* n = new int[2];
    
    float placeholder = 0;
    __m256 R1MBR1 = _mm256_set_ps(R3->boundary->at(3), R3->boundary->at(2), R3->boundary->at(1), R3->boundary->at(0), R1->boundary->at(3), R1->boundary->at(2), R1->boundary->at(1), R1->boundary->at(0));
    __m256 R2MBR1 = _mm256_set_ps(R4->boundary->at(1), R4->boundary->at(0), R4->boundary->at(3), R4->boundary->at(2), R2->boundary->at(1), R2->boundary->at(0), R2->boundary->at(3), R2->boundary->at(2));
    __m256 comparison = _mm256_cmp_ps(R1MBR1, R2MBR1, 14); //greater than comparison


    if ( isnan(comparison[0]) ||
         isnan(comparison[1]) ||
         (comparison[2] == 0) ||
         (comparison[3] == 0) )
    {
        n[0] = 0;
    }
    else
    {
        n[0] = 1;
    }

    if ( isnan(comparison[4]) ||
         isnan(comparison[5]) ||
         (comparison[6] == 0) ||
         (comparison[7] == 0) )
    {
        n[1] = 0;
    }
    else
    {
        n[1] = 1;
    }

    return n;
}

void rectUnion(vector<RectReal> *rect1, vector<RectReal> *rect2)//parallize?
{
    __m256 r1 = _mm256_set_ps(0,0,0,0,rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3));
    __m256 r2 = _mm256_set_ps(0,0,0,0,rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3));
    __m256 minvec = _mm256_min_ps(r1, r2);
    __m256 maxvec = _mm256_max_ps(r1, r2);

    rect1->at(0) = minvec[3];	rect1->at(1) = minvec[2];
    rect1->at(2) = maxvec[1];	rect1->at(3) = maxvec[0];

}

//union of 4 children mbrs
//using 256 bit registers and 32 bit floats, you can fit 2 mbrs. can compare and join 4 different mbrs quickly
//compare and union join child 1 and 2 first, then 3 and 4. then union join the union of those two using unionrect for just 2 mbrs like before
void rectUnion4(vector<RectReal> *rect1, vector<RectReal> *rect2, vector<RectReal> *rect3, vector<RectReal> *rect4)
{
    // xmin, ymin, xmax, ymax
    __m256 r1 = _mm256_set_ps(rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3), rect3->at(0), rect3->at(1), rect3->at(2), rect3->at(3));
    __m256 r2 = _mm256_set_ps(rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3), rect4->at(0), rect4->at(1), rect4->at(2), rect4->at(3));
    __m256 minvec = _mm256_min_ps(r1, r2);
    __m256 maxvec = _mm256_max_ps(r1, r2);

    __m256 result = _mm256_blend_ps(minvec, maxvec, 0x33);
    __m256 result2 = _mm256_permute2f128_ps(result, result, 0x3);
   
    __m256 vec1 = _mm256_min_ps(result, result2);
    __m256 vec2 = _mm256_max_ps(result, result2); 
    
    rect1->at(0) = vec1[3];	rect1->at(1) = vec1[2];
    rect1->at(2) = vec2[1];	rect1->at(3) = vec2[0];    
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

    //original
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

// output is a vector in cpp is arraylist in java
vector<pair<struct Node *, struct Node *> > overlappingChildren(struct Node *R1, struct Node *R2)
{
    int i, j;
    int* n;
    vector<pair<struct Node *, struct Node *> > overlappingPairs;
    // for binary tree, fanout is 2
    // fanout means maximum number of children of a node in a tree 
    for(i = 0; i < FANOUT; i += 1)  // P -> R1 root node of P
    {
        for(j = 0; j < FANOUT; j += 2) // Q -> R2 root node of Q
        {
            n = overlap(R1->children->at(i)->mbr, R2->children->at(j)->mbr, R1->children->at(i)->mbr, R2->children->at(j+1)->mbr);
          
            if(n[0])
            {
                // push_back => append to an array
                // make_pair will create a pair
                overlappingPairs.push_back(make_pair( R1->children->at(i), R2->children->at(j)));
            }
            if(n[1])
            {
                // push_back => append to an array
                // make_pair will create a pair
                overlappingPairs.push_back(make_pair( R1->children->at(i), R2->children->at(j+1)));
            }

        }
    }
    return overlappingPairs;
}

void spatialJoin(struct Node *R1, struct Node *R2, vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ)
{
    int childA, childB;
    if((R1->leaf == 1) and (R2->leaf == 1)) //if C1 and C2 are leaf Nodes
    {
        int i, j;
        for(i = R1->mbr->start; i < (R1->mbr->end); i += 1)
        {
            for(j = R2->mbr->start; j < (R2->mbr->end); j += 1) //possible bug in these for loops hitting out of bounds
            {
                int ifOverlap;
                
		tPointd p;
                char ans;

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
 				//printf("%d\n",count);
                    		//printf("p((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i+1).x, ptArrP->at(i+1).y,
                                  //                                         ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j+1).x, ptArrQ->at(j+1).y);
                	}
                }
                
            }
        }
    }
    else if(R1->leaf == 1)
    {  
        for(childB = 0; childB < FANOUT; childB += 2)
        {
            struct Node *C1 = R2->children->at(childB);
            struct Node *C2 = R2->children->at(childB+1);
            int* n = overlap(R1->mbr, C1->mbr, R1->mbr, C2->mbr);
            if(n[0])
            {
                spatialJoin(R1, C1, ptArrP, ptArrQ);
            }
            if(n[1])
            {
                spatialJoin(R1, C2, ptArrP, ptArrQ);
            }

        }
    }
    else if(R2->leaf == 1)
    {
        for(childA = 0; childA < FANOUT; childA += 2)
        {
            struct Node *C1 = R1->children->at(childA);
            struct Node *C2 = R1->children->at(childA+1);
            int* n = overlap(C1->mbr, R2->mbr, C2->mbr, R2->mbr);
            if(n[0])
            {
                spatialJoin(C1, R2, ptArrP, ptArrQ);
            }
            if(n[1])
            {
                spatialJoin(C2, R2, ptArrP, ptArrQ);
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

