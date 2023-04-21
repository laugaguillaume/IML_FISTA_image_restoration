/*=================================================================
 *
 * mex-projL1Inf.C   .MEX file that computes the l1Inf projection of a Matrix W
 *
 * The calling syntax is:
 *
 *  [B] = projL1Inf(A, C, w) 
 *
 * Inputs:
 *
 *  A: a d-times-m matrix to project into the l1Inf ball
 *
 *  w: a d-dimensional vector where the ith entry corresponds to the weight of the ith row of A
 *  
 *  C: ball bound
 *
 * Outputs:
 *
 *  B : a d-times-m matrix such that sum_i=1:d [w_i * max_k |B(i,k)| ] = c
 *
 *=================================================================*/

/* $Revision: 2$ */

#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	A_IN prhs[0]
#define	C_IN prhs[1]
#define W_IN prhs[2]


/* Output Arguments */

#define	B_OUT plhs[0]

#include "projL1Inf.c"

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray*prhs[]) { 
  double *B; 
  double *A,*c,*w; 
  int nRows,nCols; 
  
  nRows = mxGetM(A_IN); 
  nCols = mxGetN(A_IN);
  
  /* Create a matrix for the return argument */ 
  B_OUT = mxCreateDoubleMatrix(nRows, nCols, mxREAL); 
  if(B_OUT==NULL)
    printf("ERROR: NO MORE HEAP SPACE");
  
  /* Assign pointers to the various parameters */ 
  B = mxGetPr(B_OUT);
  A = mxGetPr(A_IN); 
  c = mxGetPr(C_IN);
  w = mxGetPr(W_IN);
  /* Do the actual computations in a subroutine */
  projL1Inf(B,*c,A,w,nRows,nCols); 
  return;
}


