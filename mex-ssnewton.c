/*
 * mex-ssnewton.c: .MEX file that computes the semismooth Newton iteration
 *                   for L_{1,inf} projection of a Matrix A.
 * 
 * Syntax: 
 *      [B] = myssnewton(A, C)
 * 
 * Inputs:
 *      A: d x m matrix to project int the L_{1,inf} ball.
 *      C: ball bound
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

#define A_IN prhs[0]
#define C_IN prhs[1]
#define B_OUT plhs[0]

#include "semismoothnewton.c"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[])
{
    double *B;
    double *A, *c;
    int nRows, nCols;
    
    nRows = mxGetM(A_IN);
    nCols = mxGetN(A_IN);
    
    B_OUT = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
    if(B_OUT==NULL)
        printf("ERROR: NOT ENOUGH MEMORY.");
    
    B = mxGetPr(B_OUT);
    A = mxGetPr(A_IN);
    c = mxGetPr(C_IN);
//     printf("Rows: %d, Cols: %d\n", nRows, nCols);
//     printf("Data A:\n");
//     for(int i=0; i<nRows*nCols; i++)
//         printf("%f ", A[i]);
    
    mySemismoothNewton(B, *c, A, nRows, nCols);
    return;
}        
        