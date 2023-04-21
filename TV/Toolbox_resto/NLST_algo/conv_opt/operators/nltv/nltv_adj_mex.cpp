#include "mx_image.hpp"
#include "vect_image.hpp"

#ifdef _OPENMP
#include <omp.h>
#define NUM_OF_THREADS 3
#endif


template <typename T>
void nltv_adj(Image<T> &y, const Image<T> &x, const Image<T> &w, const Image<int32_t> &idx);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check the params
    if (nrhs != 3)
        mexErrMsgTxt("This function takes 3 inputs");
    if (nlhs != 1)
        mexErrMsgTxt("This function gives 1 output");
    
    // pixel type
    typedef double T;
    
    // read the inputs
    MxImage<T>       x   ( prhs[0] );
    MxImage<T>       w   ( prhs[1] );
    MxImage<int32_t> idx ( prhs[2] );
    
    // init. the output
    int Nb = x.size(1);
    int Nr = x.size(2);
    int Nc = x.size(3);
    MxImage<T> y( plhs[0], Nr, Nc, Nb );
    
    // check data
    if( idx.cols() != Nr*Nc || idx.rows() != x.size(0) || w.total() != idx.total() )
        mexErrMsgTxt("The inputs are not compatible");
    
    // write the output
    nltv_adj(y, x, w, idx);
}


template <typename T>
void nltv_adj(Image<T> &y, const Image<T> &x, const Image<T> &w, const Image<int32_t> &idx)
{
    // get the size
    int C = x.size(1);
    int N = x.size(2)*x.size(3);
    int B = x.size(0);
    
    // scroll the image components
    #ifdef _OPENMP
    omp_set_num_threads(NUM_OF_THREADS);
    #pragma omp parallel for default(none) shared(y, x, w, idx, N, C, B, T)
    #endif
    for(int c=0; c < C; c++) 
    {
        T* p_y  = y.ptr(0,0,c);
        T* p_yt = y.ptr(0,0,c);
    
        // scroll the component pixels
        for(int t=0; t < N; t++, p_yt++)
        {
            const       T* p_x =   x.ptr(0,c,t);
            const       T* p_w =   w.ptr(0,t);
            const int32_t* p_i = idx.ptr(0,t);
            
            // compute the adjoint
            for(int b=0; b < B; b++)
            {
                T xb = p_w[b] * p_x[b];
                *p_yt += xb;
                p_y[ p_i[b] ] -= xb;
            }
        }
    }
}