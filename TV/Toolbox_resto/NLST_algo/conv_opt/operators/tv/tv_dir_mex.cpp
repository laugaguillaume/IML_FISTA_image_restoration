#include "mx_image.hpp"
#include "vect_image.hpp"

#ifdef _OPENMP
#include <omp.h>
#define NUM_OF_THREADS 3
#endif


template <typename T>
void tv_dir(Image<T> &y, const Image<T> &x);


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
    
    // init. the output
    int Nr = x.size(0);
    int Nc = x.size(1);
    int Nb = x.size(2);
    MxImage<T> y( plhs[0], 2, Nb, Nr, Nc );
       
    // write the output
    tv_dir(y, x);
}


template <typename T>
void tv_dir(Image<T> &y, const Image<T> &x)
{
    // get the size
    int Nr = x.size(0);
    int Nc = x.size(1);
    int Nb = x.size(2);
    
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
            const T* p_x = x.ptr(0,c,b);
        
            for(int r=0; r < Nr-1; r++) 
            {
                // get the pointer to the block
                T* p_y = y.ptr(0,b,r,c);
            
                // get the central pixel
                T xt = p_x[r];
          
                // compute the differences
                p_y[0] = p_x[r+Nr] - xt;
                p_y[1] = p_x[r+1]  - xt;
            }
        
            // complete the horizontal
            y(0,b,Nr-1,c) = p_x[2*Nr-1] - p_x[Nr-1];
        }
    
        // complete the vertical
        const T* p_x = x.ptr(0,Nc-1,b);
        for(int r=0; r < Nr-1; r++) 
        {           
            // get the central pixel
            T xt = p_x[r];
          
            // compute the differences
            y(1,b,r,Nc-1) = p_x[r+1] - xt;
        }
    }
}
