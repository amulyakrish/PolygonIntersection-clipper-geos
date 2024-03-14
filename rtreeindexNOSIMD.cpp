#include "rtreeindexNOSIMD.h"
#include "segseg.h"
#include <stdio.h>
#include <math.h>

int count = 0;
int candCount = 0;

double max(double num1, double num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

double min(double num1, double num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

/* if rectangle overlap with rectangle 2 */
int overlap(struct MBR *rect1, struct MBR *rect2)
{
    if ((rect1->boundary->at(0) > rect2->boundary->at(2)) || 
        (rect1->boundary->at(2) < rect2->boundary->at(0)) || 
        (rect1->boundary->at(3) < rect2->boundary->at(1)) || 
        (rect1->boundary->at(1) > rect1->boundary->at(3)))
    {
        return 0;
    }
    else
    {
        return 1;
    }
    //0: xmin, 1: ymin, 2: xmax, 3: ymax
}

int contains(tPointd p, struct MBR *R1, struct MBR *R2)
{
    double minX, minY, maxX, maxY;

    minX = max(R1->boundary->at(0), R2->boundary->at(0));
    minY = max(R1->boundary->at(1), R2->boundary->at(1));
    maxX = min(R1->boundary->at(2), R2->boundary->at(2));
    maxY = min(R1->boundary->at(3), R2->boundary->at(3));

    // find maxs of mins and mins of maxs
    //__m256d minCoord = _mm256_max_pd(r1, r2);
    //__m256d maxCoord = _mm256_min_pd(r1, r2);

    if(p[X] < minX && p[X] > maxX) 
	return 0; 
    if(p[Y] < minY && p[Y] > maxY)
	return 0;

    return 1;
}

/**
* Finds index of point intersecting overlap
* @returns index
*/ 
int findIndex(struct Node *R1, struct Node *R2, vector<struct Point> *ptArr, int index)
{
	tPointd e = {ptArr->at(index).x, ptArr->at(index).y};

	if(index == R1->mbr->end)
	{
                printf("\n**True Candidate Found**\n");
		candCount++;
		return index;
	}

        if(contains(e, R1->mbr, R2->mbr) == 0) // not in inner tile
	{		
		findIndex(R1, R2, ptArr, index++);
	}

	if(index != R1->mbr->start)
	{
                printf("\n**Candidate found**\n");
		candCount++;
		return index;
	}


	//printf("Index Found: %d\nStart Index: %d\tEnd Index: %d\n",index,R1->mbr->start,R1->mbr->end);
	return index;
}





// output is a vector in cpp is arraylist in java
vector<pair<struct Node *, struct Node *> > overlappingChildren(struct Node *R1, struct Node *R2)
{
    int i, j;
    vector<pair<struct Node *, struct Node *> > overlappingPairs;
    // for binary tree, fanout is 2
    // fanout means maximum number of children of a node in a tree 
    for(i = 0; i < FANOUT; i += 1)  // P -> R1 root node of P
    {
        for(j = 0; j < FANOUT; j += 1) // Q -> R2 root node of Q
        {
            if(overlap(R1->children->at(i)->mbr, R2->children->at(j)->mbr))
            {
                // push_back => append to an array
                // make_pair will create a pair
                overlappingPairs.push_back(make_pair( R1->children->at(i), R2->children->at(j)));
            }
        }
    }
    return overlappingPairs;
}

// output is rect1 = union of rectangle rect1 and rect2
// xmin, ymin, xmax, ymax
void rectUnion(vector<RectReal> *rect1, vector<RectReal> *rect2)//parallize?
{
    rect1->at(0) = min((double)rect1->at(0), (double)rect2->at(0));
    rect1->at(1) = min((double)rect1->at(1), (double)rect2->at(1));
    rect1->at(2) = max((double)rect1->at(2), (double)rect2->at(2));
    rect1->at(3) = max((double)rect1->at(3), (double)rect2->at(3));
}

void getRange(int *range, int id, int low, int high)
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

// in the beginning, first is 0 and last is last vertex index from ptArr vector
// create a minimum bounding rectangle for vertices in the ptArr vector
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
    
    
    t->boundary = new vector<RectReal>[4];      // make it dynamic memory allocation
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

// return union of rectangles from nodes array
struct MBR * unionJoin(vector<struct Node *> *nodes, int n)
{
    int rect, x;
    struct MBR *newMBR = (struct MBR *)malloc(sizeof(struct MBR));
    vector<RectReal> *unionRect = new vector<RectReal>[NUMSIDES]; // lowerLeft and upperRight  = 4 doubles
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(0));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(1));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(2));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(3));

    for(rect = 1; rect < n; rect += 1)
    {
        rectUnion(unionRect, (nodes->at(rect)->mbr->boundary));
    }

    newMBR->start = nodes->at(0)->mbr->start;
    newMBR->end = nodes->at(n-1)->mbr->end;
 
    newMBR->boundary = unionRect;

    return newMBR;
}

// ptArr is array of vertices
struct Node * createRTree(vector<struct Point> *ptArr,int low, int high)
{
    struct Node *root;
    if((high - low) <= BUNDLEFACTOR) // how many vertices should be combined to form a tile in a leaf node
    {
        root = createLeaf(ptArr, low, high);
    }
    else
    {
        int childID;
        root = (struct Node *)malloc(sizeof(struct Node));
        // FANOUT = how many children
        vector<struct Node *> *tempChildren = new vector<struct Node *>[FANOUT];
        
        for(childID = 0; childID < FANOUT; childID += 1)
        {
            // 100 vertices and 4 children
            // Each child node should have 25 vertices
            int range[2];
            getRange(range, childID, low, high);
            // 101 = 50, 51
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

// polygon (geometry) overlap test; divide and conquer with tile overlaps
void spatialJoin(struct Node *R1, struct Node *R2, vector<struct Point> *ptArrP, 
                             vector<struct Point> *ptArrQ)
{
    int tempi = 0;
    int tempj = 0;
    int child1, child2;
    if((R1->leaf == 1) and (R2->leaf == 1)) //if C1 and C2 are leaf Nodes
    {
        int i, j;
        

        for(i = findIndex(R1, R2, ptArrP, R1->mbr->start); i < (R1->mbr->end); i += 1)
        {
	    
            for(j = findIndex(R2, R1, ptArrQ, R2->mbr->start); j < (R2->mbr->end); j += 1) //possible bug in these for loops hitting out of bounds
            {
                tPointd p;
                char ans;
		tPointi a = {(int)ptArrP->at(i).x, (int)ptArrP->at(i).y};
                tPointi b = {(int)ptArrP->at(i+1).x, (int)ptArrP->at(i+1).y};
                tPointi c = {ptArrQ->at(j).x, ptArrQ->at(j).y};
                tPointi d = {ptArrQ->at(j+1).x, ptArrQ->at(j+1).y};         
                //printf("a, b: (%d, %d), (%d, %d) - c, d: (%d, %d), (%d, %d)\n\n", a[0],a[1], b[0], b[1], c[0],c[1], d[0], d[1]);

                ans = SegSegInt(a, b, c, d, p);
                //cout << ans << "\n";
                if (ans == '1')
                {
                    //count++;
                    //printf("%d\n", count);
                    //printf("p((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i+1).x, ptArrP->at(i+1).y,
                    //                                                        ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j+1).x, ptArrQ->at(j+1).y);
                }
                //printf("LSI on line segments p(%d, %d) and q(%d, %d)\n", R1->mbr->start, R1->mbr->end, R2->mbr->start, R2->mbr->end);
                // printf("MBR((%f, %f)(%f, %f)) overlaps MBR((%f, %f)(%f, %f))\n\n", R1->mbr->boundary->at(0), R1->mbr->boundary->at(1), R1->mbr->boundary->at(2), R1->mbr->boundary->at(3),
                
            }
        }
    }
    else if(R1->leaf == 1)
    {  
        for(child2 = 0; child2 < FANOUT; child2 += 1)
        {
            struct Node *C2 = R2->children->at(child2);
            if(overlap(R1->mbr, C2->mbr))
            {
                spatialJoin(R1, C2, ptArrP, ptArrQ);
            }
        }
    }
    else if(R2->leaf == 1)
    {
        for(child1 = 0; child1 < FANOUT; child1 += 1)
        {
            struct Node *C1 = R1->children->at(child1);
            if(overlap(C1->mbr, R2->mbr))
            {
                spatialJoin(C1, R2, ptArrP, ptArrQ);
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
/*
        for(child1 = 0; child1 < FANOUT; child1+=1) //paralize these for loops
        {
            struct Node *C1 = R1->children->at(child1);
            for(child2 = 0; child2 < FANOUT; child2+=1)
            {
                struct Node *C2 = R2->children->at(child2);
                if(overlap(C1->mbr, C2->mbr))
                {
                    spatialJoin(C1, C2);
                }
            }
        }
*/

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
            tPointi a = {ptArrP->at(i).x, ptArrP->at(i).y};
            tPointi b = {ptArrP->at(i+1).x, ptArrP->at(i+1).y};
            tPointi c = {ptArrQ->at(j).x, ptArrQ->at(j).y};
            tPointi d = {ptArrQ->at(j+1).x, ptArrQ->at(j+1).y};         

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

