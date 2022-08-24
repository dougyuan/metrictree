// tree.h
#ifndef METRICTREE_H
#define METRICTREE_H

typedef struct MetricTree Tree;

// MetricTree draws inspiration from Steve Hanov's
// implementation of a partition tree. Citation in the report.
struct MetricTree{
  // Starting (vantage) point
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
    Return: The tree
*/
Tree* build(double *X, int n, int d);
// Return vantage-point subtree with points inside radius
/*
    Input: A metric tree node
    Return: The subtree
*/
Tree* getInner(Tree* T);
// Return vantage-point subtree with points outside radius
/*
    Input: A metric tree node
    Return: The subtree
*/
Tree* getOuter(Tree* T);
// Return median of distances to vantage point
/*
    Input: A metric tree node
    Return: The subtree
*/
double getMD(Tree* T);
// Return the coordinates of the vantage point
/*
    Input: A metric tree node
    Return: The coordinates [d-dimensional vector]
*/
double * getVP(Tree* T);
// Return the index of the starting point
/*
    Input: A metric tree node
    Return: The index to the input vector of data points
*/
int getIDX(Tree* T);

#endif
