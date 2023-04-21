/*=================================================================
 *
 * projL1Inf.c	
 *
 *
 * first released on June 2009
 * revised on September 2011
 *=================================================================*/

#include <stdlib.h>
#include <math.h>


struct residual_t{
    int rowInd;
    double value;
};
        
static void printMatrix(double *m, int nRows,int nCols) {
  int i,j,index;
  double value;
  
  for(i=0;i<nRows;i++) {
    for (j=0;j<nCols;j++) {
      index=(j)*nRows + i;       
      value=m[index];
      printf("%f ",value);                
    }   
    printf("\n");
  } 
}

static double normL1Inf(double *m, int nRows, int nCols)
{
    double norm=0, value;
    int i, j, index;
    
    for(i=0;i<nRows;i++) 
    {
        double maxRow=0;
        for (j=0;j<nCols;j++)
        {
            index = j*nRows + i;
            value = fabs(m[index]);
            if(value>maxRow)	{
                maxRow=value;
            }
        }
        norm= norm + maxRow;
    }
    return norm;
}

int compare(const void *_a, const void *_b) {  
  double *a, *b;      
  a = (double *) _a;
  b = (double *) _b;              
  if(*a==*b)
    return 0;
  if(*a > *b)
    return -1;
  if(*a < *b)
    return 1;                          
}

int compareR(const void *_a, const void *_b) {
  struct residual_t *a, *b;      
  a = (struct residual_t *) _a;
  b = (struct residual_t *) _b;
  
  if((*a).value==(*b).value)
    return 0;
  if((*a).value < (*b).value)
    return -1;
  if((*a).value > (*b).value)
    return 1;
}


static void projL1Inf(double B[], double C, double A[], double w[], int nRows, int nCols) {
    
    int i, j, index;
    int rowIndex, colIndex, nonZero;
    double *maxes;
    int *sparseMap;
    double maxRow,value,l1InfNorm;
    double *sparseW;
    double *S;
    int indx;
    struct residual_t *R;
    double residual,si;
    int idxR; double Gradient; int *Ks;
    double Norm; double nextNorm;
    double theta; double *mu; double sum;
    
    maxes = (double*) malloc(sizeof(double)*nRows);
    
    for(i=0;i<nCols*nRows;i++) {
        B[i] = A[i];
    }
    
    /* 1- Compute L1InfNorm and maxes, if <= C copy matrix A to output B and return. Else move to 2 */
    
    l1InfNorm=0;
    nonZero=0;
    for(i=0;i<nRows;i++) {
        maxRow=0;
        for (j=0;j<nCols;j++) {
            index=j*nRows + i;
            value=fabs(A[index]);
            if(value>maxRow)	{
                maxRow=value;
            }
        }
        maxes[i]=maxRow;
        l1InfNorm= l1InfNorm + w[i] * maxRow;
        
        if(maxRow>0) {
            nonZero++;
        }
    }
    
    if(l1InfNorm<=C) {
        free(maxes);
        return;
    }
    
    /* 2- Copy to new matrix AP of absolute values, do this for rows with maximums larger than 0, remember the mappings between the rows of A and the rows of AP */
    
    sparseMap = (int*) malloc(sizeof(int)*nonZero);
    sparseW = (double*) malloc(sizeof(double)*nonZero);
    index=0;
    
    for(i=0;i<nRows;i++) {
        if(maxes[i]>0) {
            sparseMap[index]=i;
            sparseW[index]=w[i];
            index++;
        }
    }
    
    
    /* create nonZero rows of nCols+1 size */
    
    S = (double*) calloc((nCols +1) *nonZero, sizeof(double));
    
    /*
    for(i=0;i<(nCols+1)*nonZero;i++) {
       S[i]=0;
    }
     */
    
    for(i=0;i<nonZero;i++) {
        for (j=0; j<nCols; j++) {
            index=j*nRows + sparseMap[i];
            value=fabs(A[index]);
            index=i*(nCols+1)+j;
            S[index]=value;
        }
        /*
         * printf("S[%d]:", sparseMap[i]);
         * for (j=0;j<nCols;++j) printf(" %f", S[i*(nCols+1)+j]);
         * printf("\n");
         */
    }
    
    
    /* 3- Sort Colums of S */
    
    indx=0;
    for(i=0;i<nonZero;i++) {
        qsort(&S[indx], nCols+1, sizeof(double), &compare);
        indx=indx + nCols + 1;
    }
    
    
    /*
    printf("S = \n");
    printMatrix(S,nCols+1,nonZero);
     */
    
    /* 5- Compute Residual */
    
    R = (struct residual_t*) malloc(sizeof(struct residual_t)*nonZero*(nCols+1));
    
    indx=0;
    for (i=0;i<nonZero;i++) {
        R[indx].rowInd=i;
        R[indx].value=0;
        residual=S[indx];
        indx++;
        
        for (j=0;j<nCols;j++) {
            R[indx].rowInd=i;
            R[indx].value=residual - S[indx] * (j+1);
            residual=residual + S[indx];
            indx++;
        }
    }
    
    
    /* 6- Merge Residual */
    qsort(R, (nCols+1)*nonZero, sizeof(struct residual_t), &compareR);
    
    /*
    printf("Residuals \n");
    for (i=0;i<(nCols+1)*nonZero;i++)
    {
    printf("row index %d value %f \n", R[i].rowInd,R[i].value);
    } */
    
    
    /* 7- Enter Loop to compute theta */
    
    /* points to the current position in R */
    idxR=0;
    /* gradient of the N function, at interval [ R[idxR].value , R[idxR+1].value ] */
    Gradient=0;
    /* Ks[j] counts how many points of row j we did visit in R, this includes the first point of the current interval;
     -1/Ks[j] is the gradient of h_j in the current interval */
    Ks = (int*) malloc(sizeof(int)*nonZero);
    
    for (i=0;i<nonZero;i++) {
        Ks[i]=0;
    }
    
    /* to initialize we skip all the entries in R corresponding to 0 residual,
     for each row i there will be n_i entries, where n_i is the number of maximums of row i */
    while(R[idxR].value<=0.0) {
        Ks[R[idxR].rowInd]++;
        idxR++;
    }
    idxR=idxR-1;
    
    for (i=0;i<nonZero;i++) {
        Gradient += -sparseW[i]/(double) Ks[i];
    }
    
    
    /* norm at R[idxR].value */
    Norm=l1InfNorm;
    /* norm at next point, i.e. R[idxR+1].value */
    nextNorm = Gradient * (R[idxR+1].value - R[idxR].value ) + Norm;
    /*
    printf("INI iteration %d  row=--  norm=%f  gradient=%f  nextNorm=%f  K=[", idxR, Norm, Gradient, nextNorm);
    { int k; for(k=0;k<nonZero;++k) { printf("%d ", Ks[k]); } printf("]\n"); }*/
    
    /* keep crossing points in R until we hit the right interval */
    while (nextNorm>C) {
        int r;
        
        idxR++;
        Norm=nextNorm;
        r=R[idxR].rowInd;
        
        
        if (1) {
            // remove r-th term from the Gradient
            Gradient += sparseW[r]/(double) Ks[r];
            // increment the number of points we visited of the r-th row
            Ks[r]++;
            // update the r-th component of the gradient
            if (Ks[r]<=nCols) {
                Gradient += -sparseW[r]/(double) Ks[r];
            }
        }
        else {
            Ks[r]++;
            /* the r'th gradient becomes 0 when we have crossed nCols+1 points */
            Gradient = 0;
            for (i=0; i<nonZero; ++i)
            {
                if (Ks[i]<=nCols)
                {
                    Gradient += -sparseW[i]/(double) Ks[i];
                }
            }
        }
        
        nextNorm = Gradient * (R[idxR+1].value - R[idxR].value ) + Norm;
        
        /*
      printf("iteration %d  row=%d  norm=%f  gradient=%f  nextNorm=%f  K=[", idxR, r, Norm, Gradient, nextNorm);
      { int k; for(k=0;k<nonZero;++k) { printf("%d ", Ks[k]); } printf("]\n"); }*/
    }
    
    theta = (C - Norm + Gradient*R[idxR].value)/Gradient;
    
    
    /* 8- Compute mu_i for all i */
    
    mu = (double*) malloc(nonZero*sizeof(double));
    sum = 0;
    for (i=0; i<nonZero; ++i) {
        if (Ks[i]<=nCols) {
            int j=idxR;
            while (R[j].rowInd!=i) j--;
            mu[i] = S[i*(nCols+1)+Ks[i]-1] - (1/(double) Ks[i])*(theta-R[j].value);
        }
        else mu[i] = 0;
        sum += sparseW[i] * mu[i];
    }
    
    
    
    /*
    printf("theta: %f;  new norm: %f; \n", theta, sum);
     */
    
    /* 9- Create output matrix */
    
    /*
    printf("Ms = \n");
    printMatrix(mu,1,nonZero);
     */
    
    for(i=0;i<nRows*nCols;i++) {
        B[i]=A[i];
    }
    
    /*
     printf("A= \n");
     printMatrix(B,nRows,nCols);
     */
    
    
    for(i=0;i<nonZero;i++) {
        for(j=0;j<nCols;j++) {
            index=(j)*nRows + sparseMap[i];
            if(B[index] > mu[i] || -B[index] > mu[i]) {
                if(A[index]>=0.0)
                    B[index] = mu[i];
                else
                    B[index]=-mu[i];
            }
        }
    }
    
    /*
     printf("B = \n");
     printMatrix(B,nRows,nCols);
     */
    
    free(maxes);
    free(sparseMap);
    free(S);
    free(R);
    free(Ks);
    free(mu);
    
    return;
}

