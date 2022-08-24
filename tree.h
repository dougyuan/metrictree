// tree.h

#ifndef METRICTREE_H
#define METRICTREE_H

typedef struct MetricTree Tree;

struct MetricTree{
  // Array of size d containing the point's coordinates
  double* vp;
  // Median distance of point from all other points
  double md;
  // Index of point in original n*d matrix
  int idx;
  // Inner and outer subtrees
  Tree* inner;
  Tree* outer;
};

// Getter functions
// Build vantage-point tree given input dataset X
/*
    Input: X Input data points, stored as [n-by-d] array
    Input: n Number of data points (rows of X)
    Input: d Number of dimensions (columns of X)
    Return: The vantage-point tree
*/
Tree* buildvp(double *X, int n, int d);
// Return vantage-point subtree with points inside radius
/*
    Input: node A vantage-point tree
    Return: The vantage-point subtree
*/
Tree* getInner(Tree* T);
// Return vantage-point subtree with points outside radius
/*
    Input: node A vantage-point tree
    Return: The vantage-point subtree
*/
Tree* getOuter(Tree* T);
// Return median of distances to vantage point
/*
    Input: node A vantage-point tree
    Return: The median distance
*/
double getMD(Tree* T);
// Return the coordinates of the vantage point
/*
    Input: node A vantage-point tree
    Return: The coordinates [d-dimensional vector]
*/
double * getVP(Tree* T);
// Return the index of the vantage point
/*
    Input: node A vantage-point tree
    Return: The index to the input vector of data points
*/
int getIDX(Tree* T);

#endif
