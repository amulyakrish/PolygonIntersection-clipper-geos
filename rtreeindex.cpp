#include "rtreeindex.h"
#include "segseg.h"
#include <stdio.h>
#include <math.h>
#include <thread>
#include <future>
#include <mutex>
std::mutex console_mutex;

std::vector<std::future<void>> futures;

double max(double num1, double num2)
{
    return (num1 > num2) ? num1 : num2;
}

double min(double num1, double num2)
{
    return (num1 > num2) ? num2 : num1;
}

/**
 * Determines whether two MBRs overlap each other
 * @param two MBRs
 * @return 0 for no overlap and 1 for overlap
 */
int overlap(struct MBR *R1, struct MBR *R2)
{
    // 0: xmin, 1: ymin, 2: xmax, 3: ymax
    __m256d R1MBR1 = _mm256_set_pd(R1->boundary->at(3), R1->boundary->at(2), R1->boundary->at(1), R1->boundary->at(0));
    __m256d R2MBR1 = _mm256_set_pd(R2->boundary->at(1), R2->boundary->at(0), R2->boundary->at(3), R2->boundary->at(2));
    __m256d comparison = _mm256_cmp_pd(R1MBR1, R2MBR1, 14); // greater than comparison - 14 is the code value for gt

    if (isnan(comparison[0]) ||
        isnan(comparison[1]) ||
        (comparison[2] == 0) ||
        (comparison[3] == 0))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

vector<pair<struct Node *, struct Node *>> overlappingChildren(struct Node *R1, struct Node *R2)
{
    int i, j;
    vector<pair<struct Node *, struct Node *>> overlappingPairs;
    for (i = 0; i < FANOUT; i += 1)
    {
        for (j = 0; j < FANOUT; j += 1)
        {
            if (overlap(R1->children->at(i)->mbr, R2->children->at(j)->mbr))
            {
                overlappingPairs.push_back(make_pair(R1->children->at(i), R2->children->at(j)));
            }
        }
    }
    return overlappingPairs;
}

// xmin, ymin, xmax, ymax
void rectUnionV2(vector<RectReal> *rect1, vector<RectReal> *rect2) // parallize?
{
    __m256d r1 = _mm256_set_pd(rect1->at(0), rect1->at(1), rect1->at(2), rect1->at(3));
    __m256d r2 = _mm256_set_pd(rect2->at(0), rect2->at(1), rect2->at(2), rect2->at(3));
    __m256d minvec = _mm256_min_pd(r1, r2);
    __m256d maxvec = _mm256_max_pd(r1, r2);

    rect1->at(0) = minvec[3];
    rect1->at(1) = minvec[2];
    rect1->at(2) = maxvec[1];
    rect1->at(3) = maxvec[0];
}

void getRange(int *range, int id, int low, int high) // does this not last indext back to first index? for polygons
{
    int size = high - low;
    int start = low + (id * size) / FANOUT;
    range[0] = start;

    int end;
    if (id == (FANOUT - 1))
    {
        end = high;
    }
    else
    {
        end = low + ((id + 1) * size) / FANOUT;
    }
    range[1] = end;
}

struct MBR *createTile(vector<Point> *ptArr, int first, int last)
{
    // original version
    int index;
    Point p = ptArr->at(first);
    double minX = p.x;
    double minY = p.y;

    double maxX = minX;
    double maxY = minY;

    for (index = first; index <= last; index += 1)
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

struct Node *createLeaf(vector<Point> *ptArr, int low, int high)
{
    if ((low == high) || (low > high))
        return NULL;
    struct Node *leaf = (struct Node *)malloc(sizeof(struct Node));
    struct MBR *env = createTile(ptArr, low, high);
    leaf->mbr = env;
    leaf->leaf = 1;
    return leaf;
}

struct MBR *unionJoin(vector<struct Node *> *nodes, int n)
{
    int rect;
    struct MBR *newMBR = (struct MBR *)malloc(sizeof(struct MBR));
    vector<RectReal> *unionRect = new vector<RectReal>[NUMSIDES];
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(0));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(1));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(2));
    unionRect->push_back(nodes->at(0)->mbr->boundary->at(3));

    for (rect = 1; rect < n; rect += 2)
    {
        rectUnionV2(unionRect, (nodes->at(rect)->mbr->boundary));
    }

    newMBR->start = nodes->at(0)->mbr->start;
    newMBR->end = nodes->at(n - 1)->mbr->end;

    newMBR->boundary = unionRect;

    return newMBR;
}

struct Node *createRTreeThreadPool(vector<Point> *ptArr, int low, int high) {
    struct Node *root;
    if ((high - low) <= BUNDLEFACTOR) {
        root = createLeaf(ptArr, low, high);
    } else {
        root = (struct Node *)malloc(sizeof(struct Node));
        vector<struct Node *> *tempChildren = new vector<struct Node *>[FANOUT];

        for (int childID = 0; childID < FANOUT; ++childID) {
            int range[2];
            getRange(range, childID, low, high);
            auto future = pool.enqueue(createRTreeSequential, ptArr, range[0], range[1]);
            tempChildren->push_back(future.get());
        }

        root->children = tempChildren;
        root->leaf = 0;
        struct MBR *parentTile = unionJoin(root->children, FANOUT); // test unionJoin
        root->mbr = parentTile;
    }
    return root;
}

struct Node *createRTreeSequential(vector<Point> *ptArr, int low, int high)
{
    struct Node *root;
    if ((high - low) <= BUNDLEFACTOR)
    {
        root = createLeaf(ptArr, low, high);
    }
    else
    {
        int childID;
        root = (struct Node *)malloc(sizeof(struct Node));
        vector<struct Node *> *tempChildren = new vector<struct Node *>[FANOUT];

        // uncomment for threads
        //  int range0[2], range1[2];
        //  getRange(range0, 0, low, high);
        //  getRange(range1, 1, low, high);
        //  future<struct Node *> future1 = async(launch::async, createRTree, ptArr, range0[0], range0[1]);

        // struct Node *child2 = createRTree(ptArr, range1[0], range1[1]);
        // struct Node *child1 = future1.get();

        // tempChildren->push_back(child1);
        // tempChildren->push_back(child2);

        // comment out for loop for threads
        for (childID = 0; childID < FANOUT; childID += 1)
        {
            int range[2];
            getRange(range, childID, low, high);
            struct Node *child = createRTreeSequential(ptArr, range[0], range[1]);
            tempChildren->push_back(child);
        }
        root->children = tempChildren;
        root->leaf = 0;
        struct MBR *parentTile = unionJoin(root->children, FANOUT); // test unionJoin
        root->mbr = parentTile;
    }
    return root;
}

int contains(struct Point p, int *CMBR)
{
    if (p.x < CMBR[0] || p.x > CMBR[2])
        return 0;
    if (p.y < CMBR[1] || p.y > CMBR[3])
        return 0;

    return 1;
}

/*void processIntersectionCheck(
    const Node *R1, const Node *R2,
    const std::vector<struct Point> *ptArrP, const std::vector<struct Point> *ptArrQ,
    std::atomic<int> &count)
{
    int i, j;

    for (i = R1->mbr->start; i < (R1->mbr->end); i += 1)
    {
        tPointd e = {ptArrP->at(i).x, ptArrP->at(i).y};
        tPointd f = {ptArrP->at(i + 1).x, ptArrP->at(i + 1).y};

        for (j = R2->mbr->start; j < (R2->mbr->end); j += 1) // possible bug in these for loops hitting out of bounds
        {
            int ifOverlap;

            tPointd p;
            char ans;

            tPointd g = {ptArrQ->at(j).x, ptArrQ->at(j).y};
            tPointd h = {ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y};

            __m256d r1 = _mm256_set_pd(g[1], g[0], e[1], e[0]);
            __m256d r2 = _mm256_set_pd(h[1], h[0], f[1], f[0]);
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

            if (xmin > xmax2 || xmin2 > xmax || ymax < ymin2 || ymax2 < ymin)
                ifOverlap = 0;
            else
                ifOverlap = 1;

            if (ifOverlap == 1)
            {
                tPointi a = {(int)ptArrP->at(i).x, (int)ptArrP->at(i).y};
                tPointi b = {(int)ptArrP->at(i + 1).x, (int)ptArrP->at(i + 1).y};
                tPointi c = {(int)ptArrQ->at(j).x, (int)ptArrQ->at(j).y};
                tPointi d = {(int)ptArrQ->at(j + 1).x, (int)ptArrQ->at(j + 1).y};

                ans = SegSegInt(a, b, c, d, p);
                if (ans == '1')
                {
                    // intersection found
                    count++;
                    std::cout << "\nOverlap #" << count.load() << std::endl;
                    {
                        std::lock_guard<std::mutex> guard(print_mutex);
                        printf("\np((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i + 1).x, ptArrP->at(i + 1).y,
                               ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y);
                    }
                } // end if
            }     // end LSMBR
        }         // end j loop
    }
}*/
void processIntersectionCheck(
    const Node *R1, const Node *R2,
    const std::vector<struct Point> *ptArrP, const std::vector<struct Point> *ptArrQ,
    std::atomic<int> &count)
{
    int i, j;

    for (i = R1->mbr->start; i < (R1->mbr->end); i += 1)
    {
        tPointd e = {ptArrP->at(i).x, ptArrP->at(i).y};
        tPointd f = {ptArrP->at(i + 1).x, ptArrP->at(i + 1).y};

        for (j = R2->mbr->start; j < (R2->mbr->end); j += 1) // possible bug in these for loops hitting out of bounds
        {
            int ifOverlap;

            tPointd p;
            char ans;

            tPointd g = {ptArrQ->at(j).x, ptArrQ->at(j).y};
            tPointd h = {ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y};

            __m256d r1 = _mm256_set_pd(g[1], g[0], e[1], e[0]);
            __m256d r2 = _mm256_set_pd(h[1], h[0], f[1], f[0]);
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

            if (xmin > xmax2 || xmin2 > xmax || ymax < ymin2 || ymax2 < ymin)
                ifOverlap = 0;
            else
                ifOverlap = 1;

            if (ifOverlap == 1)
            {
                tPointi a = {(int)ptArrP->at(i).x, (int)ptArrP->at(i).y};
                tPointi b = {(int)ptArrP->at(i + 1).x, (int)ptArrP->at(i + 1).y};
                tPointi c = {(int)ptArrQ->at(j).x, (int)ptArrQ->at(j).y};
                tPointi d = {(int)ptArrQ->at(j + 1).x, (int)ptArrQ->at(j + 1).y};

                ans = SegSegInt(a, b, c, d, p);
                int currentCount = 0;
                if (ans == '1')
                {
                    // intersection found
                    /*{
    //std::lock_guard<std::mutex> guard(count_mutex);
    currentCount = ++count; // If count is std::atomic, else just use count++
    // Print the rest of your message here
}*/
 currentCount = ++count; // If count is std::atomic, else just use count++
                    //count++;
                    /*std::cout << "\nOverlap #" << count.load() << std::endl;
                    {
                        std::lock_guard<std::mutex> guard(print_mutex);
                        printf("\np((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i + 1).x, ptArrP->at(i + 1).y,
                               ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y);
                    }*/
                    #ifndef NDEBUG // If NDEBUG is not defined (debug mode)
                    {
    std::lock_guard<std::mutex> guard(console_mutex);
   // std::cout << "Overlap #" << count.load() << std::endl;
    std::cout << "Overlap #" << currentCount << "\n";
    std::cout << "p((" << ptArrP->at(i).x << ", " << ptArrP->at(i).y << "), (" 
              << ptArrP->at(i + 1).x << ", " << ptArrP->at(i + 1).y << ")) intersects q((" 
              << ptArrQ->at(j).x << ", " << ptArrQ->at(j).y << ")(" 
              << ptArrQ->at(j + 1).x << ", " << ptArrQ->at(j + 1).y << "))\n";
}
#endif
                } // end if
            }     // end LSMBR
        }         // end j loop
    }
}

void spatialJoin(struct Node *R1, struct Node *R2, vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ)
{
    int child1, child2;
    if ((R1->leaf == 1) and (R2->leaf == 1)) // if C1 and C2 are leaf Nodes
    {
       
       auto future = pool.enqueue([=] { 
        processIntersectionCheck(R1, R2, ptArrP, ptArrQ,count);
    });
    futures.push_back(std::move(future));
    return;

    /*int child1, child2;
    if ((R1->leaf == 1) and (R2->leaf == 1)) // if C1 and C2 are leaf Nodes
    {
        // parallelize the checking intersection part
       
       /* 
        

        pool.enqueue([&]()
                     { processIntersectionCheck(R1, R2, ptArrP, ptArrQ,count); }); */

        /*sequential*/
        
  /*      int i, j;

    for (i = R1->mbr->start; i < (R1->mbr->end); i += 1)
    {
        tPointd e = {ptArrP->at(i).x, ptArrP->at(i).y};
        tPointd f = {ptArrP->at(i + 1).x, ptArrP->at(i + 1).y};

        for (j = R2->mbr->start; j < (R2->mbr->end); j += 1) // possible bug in these for loops hitting out of bounds
        {
            int ifOverlap;

            tPointd p;
            char ans;

            tPointd g = {ptArrQ->at(j).x, ptArrQ->at(j).y};
            tPointd h = {ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y};

            __m256d r1 = _mm256_set_pd(g[1], g[0], e[1], e[0]);
            __m256d r2 = _mm256_set_pd(h[1], h[0], f[1], f[0]);
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

            if (xmin > xmax2 || xmin2 > xmax || ymax < ymin2 || ymax2 < ymin)
                ifOverlap = 0;
            else
                ifOverlap = 1;

            if (ifOverlap == 1)
            {
                tPointi a = {(int)ptArrP->at(i).x, (int)ptArrP->at(i).y};
                tPointi b = {(int)ptArrP->at(i + 1).x, (int)ptArrP->at(i + 1).y};
                tPointi c = {(int)ptArrQ->at(j).x, (int)ptArrQ->at(j).y};
                tPointi d = {(int)ptArrQ->at(j + 1).x, (int)ptArrQ->at(j + 1).y};

                ans = SegSegInt(a, b, c, d, p);
                if (ans == '1')
                {
                    // intersection found
                    count++;
                    std::cout << "\nOverlap #" << count << std::endl;
                    printf("\np((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i + 1).x, ptArrP->at(i + 1).y,
                               ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y);
                    
                } // end if
            }     // end LSMBR
        }         // end j loop
    }*/
    }
    else if (R1->leaf == 1)
    {
        // printf("R1 Leaf R2 Not Leaf\n"); //comment

        for (child2 = 0; child2 < FANOUT; child2 += 1)
        {
            struct Node *C2 = R2->children->at(child2);
            if (overlap(R1->mbr, C2->mbr))
            {
                spatialJoin(R1, C2, ptArrP, ptArrQ);
            }
        }
    }
    else if (R2->leaf == 1)
    {

        // printf("R1 Not Leaf R2 Leaf\n"); //comment
        for (child1 = 0; child1 < FANOUT; child1 += 1)
        {
            struct Node *C1 = R1->children->at(child1);
            if (overlap(C1->mbr, R2->mbr))
            {
                spatialJoin(C1, R2, ptArrP, ptArrQ);
            }
        }
    }
    else
    {

        // printf("R1 Not Leaf R2 Not Leaf\n"); //comment
        int i;
        vector<pair<struct Node *, struct Node *>> candidatePairs = overlappingChildren(R1, R2);
        for (i = 0; i < (int)candidatePairs.size(); ++i)
        {
            spatialJoin(candidatePairs.at(i).first, candidatePairs.at(i).second, ptArrP, ptArrQ);
        }
    }
}

/*
int contains(Point p, struct MBR *R1, struct MBR *R2)
{
    double minX = max(R1->boundary->at(0), R2->boundary->at(0));
    double minY = max(R1->boundary->at(1), R2->boundary->at(1));
    double maxX = min(R1->boundary->at(2), R2->boundary->at(2));
    double maxY = min(R1->boundary->at(3), R2->boundary->at(3));

    printf("minx: %f\tminy: %f\tmaxx: %f\tmaxy: %f\n", minX,minY,maxX,maxY);

    if(p.x < minX || p.x > maxX)
    return 0;
    if(p.y < minY || p.y > maxY)
    return 0;

    return 1;
}
*/
int isOverlap(int i, int j, vector<struct Point> *P, vector<struct Point> *Q)
{
    // Version 1: Create MBR's for the four points and pass through overlap function

    vector<RectReal> *b1 = new vector<RectReal>[4];
    vector<RectReal> *b2 = new vector<RectReal>[4];

    b1->push_back(min(P->at(i).x, P->at(i + 1).x));
    b1->push_back(min(P->at(i).y, P->at(i + 1).y));
    b1->push_back(max(P->at(i).x, P->at(i + 1).x));
    b1->push_back(max(P->at(i).y, P->at(i + 1).y));

    b2->push_back(min(Q->at(j).x, Q->at(j + 1).x));
    b2->push_back(min(Q->at(j).y, Q->at(j + 1).y));
    b2->push_back(max(Q->at(j).x, Q->at(j + 1).x));
    b2->push_back(max(Q->at(j).y, Q->at(j + 1).y));

    MBR R1(0, 3, b1);
    MBR R2(0, 3, b2);

    int result = overlap(&R1, &R2);
    return result;

    /*
        // Version 2: Find mins and maxes. Directly use overlap function logic.

        __m256d r1 = _mm256_set_pd(Q->at(j).y,Q->at(j).x, P->at(i).y, P->at(i).x);
        __m256d r2 = _mm256_set_pd(Q->at(j+1).y,Q->at(j+1).x, P->at(i+1).y, P->at(i+1).x);
        __m256d minvec = _mm256_min_pd(r1, r2);
        __m256d maxvec = _mm256_max_pd(r1, r2);

        __m256d R1MBR1 = _mm256_set_pd(maxvec[1], maxvec[0], minvec[1], minvec[0]);
        __m256d R2MBR1 = _mm256_set_pd(minvec[2], minvec[3], maxvec[3], maxvec[2]);
        __m256d comparison = _mm256_cmp_pd(R1MBR1, R2MBR1, 14); //greater than comparison

        if ( isnan(comparison[0]) ||
             isnan(comparison[1]) ||
             (comparison[2] == 0) ||
             (comparison[3] == 0) )
        {
            return 0;
        }
        else
        {
            return 1;
        }
    */
}

void lsMBRfilter(vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ)
{
    printf("\n\nBrute force LSI\n");
    int i, j, ifOverlap;
    for (i = 0; i < (int)ptArrP->size() - 1; i += 1)
    {
        for (j = 0; j < (int)ptArrQ->size() - 1; j += 1)
        {
            tPointd e = {ptArrP->at(i).x, ptArrP->at(i).y};
            tPointd f = {ptArrP->at(i + 1).x, ptArrP->at(i + 1).y};
            tPointd g = {ptArrQ->at(j).x, ptArrQ->at(j).y};
            tPointd h = {ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y};

            __m256d r1 = _mm256_set_pd(g[1], g[0], e[1], e[0]);
            __m256d r2 = _mm256_set_pd(h[1], h[0], f[1], f[0]);
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

            if (xmin > xmax2 || xmin2 > xmax || ymax < ymin2 || ymax2 < ymin)
                ifOverlap = 0;
            else
                ifOverlap = 1;

            if (ifOverlap == 1)
            {
                tPointd p;
                char ans;
                tPointi a = {(int)ptArrP->at(i).x, (int)ptArrP->at(i).y};
                tPointi b = {(int)ptArrP->at(i + 1).x, (int)ptArrP->at(i + 1).y};
                tPointi c = {(int)ptArrQ->at(j).x, (int)ptArrQ->at(j).y};
                tPointi d = {(int)ptArrQ->at(j + 1).x, (int)ptArrQ->at(j + 1).y};

                ans = SegSegInt(a, b, c, d, p);
                if (ans == '1')
                {
                    count++;
                    // printf("%d\n",count);
                    printf("p((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i + 1).x, ptArrP->at(i + 1).y,
                           ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y);
                }
            }
        }
    }
}

// used for verifying correctness
void checkLSI(vector<struct Point> *ptArrP, vector<struct Point> *ptArrQ)
{
    printf("\n\nBrute force LSI\n");
    int i, j;
    for (i = 0; i < (int)ptArrP->size() - 1; i += 1)
    {
        for (j = 0; j < (int)ptArrQ->size() - 1; j += 1)
        {
            tPointd p;
            char ans;
            tPointi a = {(int)ptArrP->at(i).x, (int)ptArrP->at(i).y};
            tPointi b = {(int)ptArrP->at(i + 1).x, (int)ptArrP->at(i + 1).y};
            tPointi c = {(int)ptArrQ->at(j).x, (int)ptArrQ->at(j).y};
            tPointi d = {(int)ptArrQ->at(j + 1).x, (int)ptArrQ->at(j + 1).y};

            ans = SegSegInt(a, b, c, d, p);
            if (ans == '1')
            {
                // count++;
                //  printf("%d\n",count);

                // printf("p((%f, %f), (%f, %f)) intersects q((%f, %f)(%f, %f))\n", ptArrP->at(i).x, ptArrP->at(i).y, ptArrP->at(i + 1).x, ptArrP->at(i + 1).y,
                //        ptArrQ->at(j).x, ptArrQ->at(j).y, ptArrQ->at(j + 1).x, ptArrQ->at(j + 1).y);
            }
        }
    }
}

void printNode(struct Node *node)
{
    struct MBR *tile = node->mbr;
    printf("%d : %d tile coordinates: (%f, %f)-(%f, %f)\n", tile->start, tile->end, tile->boundary->at(0), tile->boundary->at(1),
           tile->boundary->at(2), tile->boundary->at(3));
}
