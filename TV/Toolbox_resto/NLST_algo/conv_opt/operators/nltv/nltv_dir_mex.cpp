#include "mx_image.hpp"
#include "vect_image.hpp"


#ifdef _OPENMP
#include <omp.h>
#define NUM_OF_THREADS 3
#endif


template <typename T>
void nltv_dir(Image<T> &y, const Image<T> &x, const Image<T> &w, const Image<int32_t> &idx);


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
    int Nr = x.size(0);
    int Nc = x.size(1);
    int Nb = x.size(2);
    int B = idx.rows();
    MxImage<T> y( plhs[0], B, Nb, Nr, Nc );
    
    // check data
    if( idx.cols() != x.size(0)*x.size(1) || w.total() != idx.total() )
        mexErrMsgTxt("The inputs are not compatible");
    
    // write the output
    nltv_dir(y, x, w, idx);
}


template <typename T>
void nltv_dir(Image<T> &y, const Image<T> &x, const Image<T> &w, const Image<int32_t> &idx)
{
    // get the size
    int N = x.size(0) * x.size(1);
    int C = x.size(2);
    int B = idx.rows();
    
    // scroll the image components
    #ifdef _OPENMP
    omp_set_num_threads(NUM_OF_THREADS);
    #pragma omp parallel for default(none) shared(y, x, w, idx, N, C, B, T)
    #endif
    for(int c=0; c < C; c++) 
    {
        const T* p_x = x.ptr(0,0,c);
        
        // scroll the component pixels
        for(int t=0; t < N; t++) 
        {
            // get the central pixel
            T xt = p_x[t];

            // get the pointers
                        T* p_y =   y.ptr(0,c,t);
            const       T* p_w =   w.ptr(0,t);
            const int32_t* p_i = idx.ptr(0,t);

            // compute the differences
            for(int b=0; b < B; b++) {
                p_y[b] = p_w[b] * (xt - p_x[ p_i[b] ]);
            }
        }
    }
}