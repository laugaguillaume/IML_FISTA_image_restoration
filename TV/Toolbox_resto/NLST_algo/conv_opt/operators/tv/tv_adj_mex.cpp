#include "mx_image.hpp"
#include "vect_image.hpp"

#ifdef _OPENMP
#include <omp.h>
#define NUM_OF_THREADS 3
#endif


template <typename T>
void tv_adj(Image<T> &y, const Image<T> &x);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check the params
    if (nrhs != 1)
        mexErrMsgTxt("This function takes 1 input");
    if (nlhs != 1)
        mexErrMsgTxt("This function must give 1 output");
    
    // pixel type
    typedef double T;
    
    // read the inputs
    MxImage<T> x( prhs[0] );
    
    // check data
    if( x.rows() != 2 )
        mexErrMsgTxt("The function expects 2 bands.");
    
    // init. the output
    int Nb = x.size(1);
    int Nr = x.size(2);
    int Nc = x.size(3);
    MxImage<T> y( plhs[0], Nr, Nc, Nb );
    
    // write the output
    tv_adj(y, x);
}


template <typename T>
void tv_adj(Image<T> &y, const Image<T> &x)
{
    // get the size
    int Nb = x.size(1);
    int Nr = x.size(2);
    int Nc = x.size(3);
    
    // scroll the spectral bands
    #ifdef _OPENMP
    omp_set_num_threads(NUM_OF_THREADS);
    #pragma omp parallel for default(none) shared(y, x, Nr, Nc, Nb)
    #endif
    for(int b=0; b < Nb; b++) 
    {
        // scroll the image
        for(int c=0; c < Nc-1; c++) 
        {
            T* p_y = y.ptr(0,c,b);
            for(int r=0; r < Nr-1; r++) 
            {
                const T* p_x = x.ptr(0,b,r,c);
        
                // compute the adjoint
                p_y[r]    -= p_x[0] + p_x[1];
                p_y[r+Nr] += p_x[0];          // horizontal
                p_y[r+1]  += p_x[1];          // vertical
            }
        
            // complete the horizontal
            p_y[  Nr-1] -= x(0,b,Nr-1,c);
            p_y[2*Nr-1] += x(0,b,Nr-1,c);
        }
    
        // complete the vertical
        T* p_y = y.ptr(0,Nc-1,b);
        for(int r=0; r < Nr-1; r++) 
        {           
            p_y[r]   -= x(1,b,r,Nc-1);
            p_y[r+1] += x(1,b,r,Nc-1);
        }
        
    }
}