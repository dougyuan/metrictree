#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "tree.h"

// Threshold of points to switch to sequential execution
#define POINT_THRESHOLD 1000
// Threshold of maximum simultaneous live threads
#define THREADS_MAX 16
// Development flags to switch execution mode from serial to parallel for distance calculation and subtree creation
#define PARALLELDIS true
#define PARALLELSUB true

// Function Prototypes
Tree* build(double *X, int n, int d);
Tree* getInner(Tree* T);
Tree* getOuter(Tree* T);
double getMD(Tree* T);
double * getVP(Tree* T);
int getIDX(Tree* T);
Tree* build_tree(double *points, int *ids, int n, int d);
void *euclidean(void *arg);
void swap(double *a, double *b);
int partition (double arr[], int low, int high);
double quickselect_median(double arr[], int length);
double quickselect(double arr[], int length, int idx);
void *build_tree_wrapper(void *arg);

// Mutex and counter to keep track of live threads
pthread_mutex_t threadMutex = PTHREAD_MUTEX_INITIALIZER;
int threadCount = 0;

// Struct used to pass arguments to threads for distance calculation
typedef struct
{
    int tid, n, d;
    double *point, *points, *distances;
} dtargs;

// Struct used to pass arguments to threads for subtree creation
typedef struct
{
    int tid, n, d;
    double *points;
    int *ids;
    Tree**subtree;
} stargs;

// Wrapper function to use build_tree() with pthreads
void *build_tree_wrapper(void *arg)
{
    *((stargs *)arg)->subtree = build_tree(((stargs *)arg)->points, ((stargs *)arg)->ids, ((stargs *)arg)->n, ((stargs *)arg)->d);
    return 0;
}

// Function to safely alter the live thread count
void modThreadCount(int n){
    pthread_mutex_lock( &threadMutex );
    threadCount += n;
    pthread_mutex_unlock( &threadMutex );
}

// Application entry point
Tree* build(double *X, int n, int d)
{
    // Allocate space for the index array
    int *ids = calloc(n, sizeof(int));

    // Build the initial ids array
    for (int i = 0; i < n; i++)
        ids[i] = i;

    // Call build_tree to get the pointer to the root of the tree
    return build_tree(X, ids, n, d);
}

// Function that recursively builds the binary tree and returns a pointer to its root
Tree* build_tree(double *points, int *ids, int n, int d)
{

    // Create node to be returned
    Tree* node = calloc(1, sizeof(Tree));

    // Check to end recursion. If points contains only one point we are returning a leaf
    if (n == 1)
    {
        // Build node
        node->inner = NULL;
        node->outer = NULL;
        node->idx = ids[0];
        node->md = 0;
        node->vp = calloc(d, sizeof(double));
        memcpy(node->vp, points, sizeof(double) * d);

        // Free memory for ids and points arrays
        free(ids);
        free(points);

        // Return node
        return node;
    }

    // Choose the last point in X as the vantage point
    double *point = (points + (n-1)*d);
    double id = ids[n-1];

    // Create array that holds euclidean distance of point from all other points
    double *distances = calloc(n-1, sizeof(double));

    // Dynamically decide number of threads for distance calculation
    int threads = 1;
    if(THREADS_MAX - threadCount >= 16){
        threads = 16;
    }
    else if(THREADS_MAX - threadCount >= 8){
        threads = 4;
    }
    else if(THREADS_MAX - threadCount >= 4){
        threads = 2;
    }
    else if(THREADS_MAX - threadCount >= 2){
        threads = 1;
    }

    // Declare pointers for pthreads and argument structure arrays
    pthread_t * disThread;
    dtargs * disArg;

    // Block size * threads = points number (n-1), blockSize = Points per thread
    int blockSize = floor((float)(n-1) / threads);

    // Calculate distances in parallel if true, else do it sequentially
    if((n-1 > POINT_THRESHOLD) && (PARALLELDIS == true) && (threads > 1))
    {
        // Increment live thread count
        modThreadCount(threads - 1);

        // Allocate memory for threads
        disThread = calloc(threads - 1, sizeof(pthread_t));
        disArg = calloc(threads, sizeof(dtargs));

        // Create threads and assign work to them
        for (int i = 0; i < threads - 1; i++)
        {
            disArg[i].d = d;
            disArg[i].point = point;
            disArg[i].tid = i;
            disArg[i].points = points + i * blockSize * d;
            disArg[i].distances = distances + i * blockSize;
            disArg[i].n = blockSize;
            pthread_create(&disThread[i], NULL, euclidean, (void *)&disArg[i]);
        }

        // Work for the main thread
        disArg[threads - 1].d = d;
        disArg[threads - 1].point = point;
        disArg[threads - 1].tid = threads - 1;
        disArg[threads - 1].points = points + (threads - 1) * blockSize * d;
        disArg[threads - 1].distances = distances + (threads - 1) * blockSize;
        disArg[threads - 1].n = (n-1) - blockSize * (threads - 1);
        euclidean((void *)&disArg[threads - 1]);

        // Join threads and decrement live threads count
        for (int i = 0; i < threads - 1; i++)
        {
            pthread_join(disThread[i], NULL);
            modThreadCount(-1);
        }
        free(disThread);
    }
    else
    {
        disArg = calloc(1, sizeof(dtargs));
        disArg[0].d = d;
        disArg[0].point = point;
        disArg[0].tid = 0;
        disArg[0].points = points;
        disArg[0].n = n-1;
        disArg[0].distances = distances;
        euclidean((void *)&disArg[0]);
    }
    // Free memory
    free(disArg);

    // At this point distances[i] indicates the distance of point i in points from the vantage point
    // Find median by creating a copy of distances and passing it to QuickSelect
    double *distancesCopy = calloc(n-1, sizeof(double));
    memcpy(distancesCopy, distances, sizeof(double) * (n-1));
    double median = quickselect_median(distancesCopy, n-1);
    free(distancesCopy);

    // Sort points into two new arrays
    // Calculate array sizes for subtrees. Values up to and equal to the median go on the inner tree
    int innerLength = 0;
    for (int i = 0; i < n-1; i++)
    {
        if(distances[i] <= median)
        {
            innerLength++;
        }
    }
    int outerLength = n - 1 - innerLength;

    // Pointers to keep track of inner and outer arrays content while sorting points
    int innerPointer = 0;
    int outerPointer = 0;

    // Create arrays for inner and outer points. Arrays contain the points and a list of ids (one for each point)
    double *innerPoints = calloc(innerLength * d, sizeof(double));
    double *outerPoints = calloc(outerLength * d, sizeof(double));
    int *innerIDs = calloc(innerLength, sizeof(int));
    int *outerIDs = calloc(outerLength, sizeof(int));

    // Sort points to inner and outer subtree
    for (int i = 0; i < n-1; i++)
    {
        if(distances[i] <= median)
        {
            memcpy(innerPoints + innerPointer * d, points + i*d, sizeof(double) * d);
            innerIDs[innerPointer] = ids[i];
            innerPointer++;
        }
        else
        {
            memcpy(outerPoints + outerPointer * d, points + i*d, sizeof(double) * d);
            outerIDs[outerPointer] = ids[i];
            outerPointer++;
        }
    }

    // Set node fields
    // Copy the point into vp because we will call free(points) that will also free(point)
    node->vp = calloc(d, sizeof(double));
    node->md = median;
    memcpy(node->vp, point, sizeof(double) * d);
    node->idx = id;

    // De-allocate unused memory
    free(points);
    free(distances);
    free(ids);

    // Thread and stargs arrays for parallel subtree threads
    pthread_t subThread;
    stargs* subArg;

    // Build subtrees in parallel or sequentially
    if((PARALLELSUB == true) && (THREADS_MAX - threadCount >= 2))
    {

        // Create threads
        if((innerLength > 0) && (outerLength > 0))
        {
            // Start inner tree creation on a thread
            modThreadCount(1);
            subArg = malloc(sizeof(stargs));
            subArg->d = d;
            subArg->n = innerLength;
            subArg->subtree = &node->inner;
            subArg->tid = 0;
            subArg->points = innerPoints;
            subArg->ids = innerIDs;
            pthread_create(&subThread, NULL, build_tree_wrapper, (void *)subArg);

            // Run outer tree creation in the main thread
            node->outer = build_tree(outerPoints, outerIDs, outerLength, d);

            // Join thread
            pthread_join(subThread, NULL);
            modThreadCount(-1);
        }
        else if(innerLength > 0)
        {
        	node->inner = build_tree(innerPoints, innerIDs, innerLength, d);
        }
        else if(outerLength > 0)
        {
        	node->outer = build_tree(outerPoints, outerIDs, outerLength, d);
        }
    }
    else
    {
        if(innerLength > 0)
        {
            node->inner = build_tree(innerPoints, innerIDs, innerLength, d);
        }
        if(outerLength > 0)
        {
            node->outer = build_tree(outerPoints, outerIDs, outerLength, d);
        }
    }

    // Check if a subtree is empty and assign null to its pointer
    if(innerLength < 1)
    {
        node->inner = NULL;
    }
    if(outerLength < 1)
    {
        node->outer= NULL;
    }

    return node;
}

// Return vantage-point subtree with points inside radius
Tree* getInner(Tree* T)
{
    return T->inner;
}

// Return vantage-point subtree with points outside radius
Tree* getOuter(Tree* T)
{
    return T->outer;
}

// Return median of distances to vantage point
double getMD(Tree* T)
{
    return T->md;
}

// Return the coordinates of the vantage point
double * getVP(Tree* T)
{
    return T->vp;
}

// Return the index of the vantage point
int getIDX(Tree* T)
{
    return T->idx;
}

// Calculates the distances of all points from point and writes them to an array
void *euclidean(void *arg)
{
    // Retrieve variables by casting the argument to a dtarg struct pointer
    double *point = ((dtargs *) arg) -> point;
    double *points = ((dtargs *) arg) -> points;
    double *distances = ((dtargs *) arg) -> distances;
    int n = ((dtargs *) arg) -> n;
    int d = ((dtargs *) arg) -> d;

    double accumulator = 0;

    for (int i = 0; i < n; i++)
    {
        accumulator = 0;
        for (int j = 0; j < d; j++)
        {
            accumulator += (point[j] - *(points + i * d + j)) * (point[j] - *(points + i * d + j));
        }
        distances[i] = sqrt(accumulator);
    }
    return 0;
}

// A utility function to swap two elements
void swap(double *a, double *b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

// QuickSort Partition function. low and high define the range of indexes in arr where partition should work
int partition (double arr[], int low, int high)
{

    // Select a pivot and initialize flag to position of smallest element before pivot
    double pivot = arr[high];
    int i = (low - 1);

    // Go through the array examining each element
    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot, increment i and swap it out with the one currently pointed by i
        if (arr[j] < pivot)
        {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }

    // Finally place pivot in its correct position in the array and return the position as the middle point
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

// Returns the median using the QuickSelect algorithm
double quickselect_median(double arr[], int length)
{

    if (length % 2 == 1)
    {
        return quickselect(arr, length, (length+1)/2);
    }
    else
    {
        return 0.5 * (quickselect(arr, length, length/2) + quickselect(arr, length, length/2 + 1));
    }
}

// Returns the idx-th element of arr when arr is sorted
// idx is the index (starting from 1) of the point we want to find when the array is sorted. For the median idx should be the middle one (i.e (length+1)/2 for odd lengths etc)
double quickselect(double arr[], int length, int idx)
{

    // Check to end recursion
    if (length == 1)
    {
        return arr[0];
    }

    // Select last array element as pivot
    double pivot = arr[length - 1];
    // Get index of pivot after we partition the array
    int pivotIndex = partition(arr, 0, length - 1);

    // Create the higher and lower arrays that occur after partitioning in QuickSort fashion
    int lowerLength = pivotIndex;
    pivotIndex++;
    int higherLength = (length - (lowerLength + 1));

    // At this point pivotIndex, lowerLength and higherLength all start from 1 not 0
    double *lower = calloc(lowerLength, sizeof(double));
    double *higher = calloc(higherLength, sizeof(double));
    memcpy(lower, arr, sizeof(double) * lowerLength);
    memcpy(higher, arr + pivotIndex, sizeof(double) * higherLength);

    // Variable to store result of following recursive calls
    double result = 0;

    // This means that the point we're looking (median in our case) is in the lower partition
    if (idx <= lowerLength)
    {
        result = quickselect(lower, lowerLength, idx);
    }
    // This means that the median is our pivot point
    else if(idx == pivotIndex)
    {
        result = pivot;
    }
    // This means that the median is in the higher partition
    else
    {
        result = quickselect(higher, higherLength, idx - pivotIndex);
    }

    // Free memory allocated to lower and higher
    free(lower);
    free(higher);

    // Return result
    return result;
}

int main(int argc, char *argv[]) {

    if(argc != 3)   {printf("Incorrect Usage. Quitting...\n"); exit(EXIT_FAILURE);}
    int testPoints = atoi(argv[1]);
    int testDimensions = atoi(argv[2]);  
    // Intialize random number generator
    time_t t;
    srand((unsigned) time(&t));

    // benchmark();

    // Create a random X array
    double *X = calloc(testPoints* testDimensions, sizeof(double));
    for(int i = 0; i < testPoints* testDimensions; i++)
        X[i] = rand();

    clock_t start, end;
    start = clock();
    Tree* tree = build(X, testPoints, testDimensions);
    end = clock();
    printf("------------ Benchmarking Results -------------\n");
    printf("Start point: %f\n", *(tree->vp));
    printf("Found median: %f\n", tree->md);
    double elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time elapsed: %fs\n", elapsed);

    return 0;
}

