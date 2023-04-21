// (2015)
//
// Author :
// Giovanni Chierchia (chierchi@telecom-paristech.fr)
//
// Contributors :
// Nelly Pustelnik (nelly.pustelnik@ens-lyon.fr)
// Jean-Christophe Pesquet (jean-christophe.pesquet@univ-paris-est.fr)
// Béatrice Pesquet (beatrice.pesquet@telecom-paristech.fr)
//
// This software contains some image processing algorithms whose purpose is to be
// used primarily for research.
//
// This software is governed by the CeCILL B license under French law and
// abiding by the rules of distribution of free software. You can use,
// modify and/ or redistribute the software under the terms of the CeCILL
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty and the software's author, the holder of the
// economic rights, and the successive licensors have only limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading, using, modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean that it is complicated to manipulate, and that also
// therefore means that it is reserved for developers and experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and, more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL B license and that you accept its terms.
// January 2015 Giovanni Chierchia
// Non-local Total-Variation 
// for multicomponent (color, multispectral, hyperspectral,...)
// image restoration
//
// This toolbox implements the algorithm presented in the paper: 
// G. Chierchia, N. Pustelnik, B. Pesquet-Popescu, J.-C. Pesquet, "A Non-Local 
// Structure Tensor Based Approach for Multicomponent Image Recovery Problems",
// IEEE Trans. on Image Process., 2014#include "src/project_L2_epi.hpp"

#include "src/prox_max.hpp"
#include "mx_image.hpp"

#ifdef _OPENMP
#include <omp.h>
#define NUM_OF_THREADS 4
#endif


template <typename T>
void prox_max(Image<T> &p, const Image<T> &y, const Image<T> &v);
template <typename T>
void prox_max(Image<T> &p, const Image<T> &y, const Image<T> &v, const Image<T> &t);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check the params
    if (nrhs < 2 || nrhs > 3)
        mexErrMsgTxt("This function takes 2 or 3 inputs");
    if (nlhs != 1)
        mexErrMsgTxt("This function gives 1 output");
 
    // pixel type
    typedef double T;
    
    // read the inputs
    MxImage<T> y( prhs[0] );
    MxImage<T> v( prhs[1] );
    
    // check the inputs
    int B = v.rows();
    int N = y.total();
    if( N*B != v.total() )
        mexErrMsgTxt("The inputs are not compatible");
    
    // init. the output
    MxImage<T> p( plhs[0], y.size(0), y.size(1), y.size(2), y.size(3) );
    
    // write the output
    if(nrhs == 2) {
        prox_max(p, y, v);
    }
    else 
    {
        MxImage<T> t( prhs[2] );    // read the 3rd input
        
        if( t.total() != v.total() )
            mexErrMsgTxt("The inputs are not compatible");
    
        prox_max(p, y, v, t);
    }
}


template <typename T>
void prox_max(Image<T> &p, const Image<T> &y, const Image<T> &v)
{
    // get the size
    int B = v.rows();
    int N = y.total();
    
    // scroll the image pixels
    #ifdef _OPENMP
    omp_set_num_threads(NUM_OF_THREADS);
    #pragma omp parallel for default(none) shared(p, y, v, N, B, T)
    #endif
    for(int i=0; i < N; i++) 
    {
        const T* v_col = v.ptr(0,i);
        
        // project the block
        p(i) = prox_max( y(i), v_col, B );
    }
}


template <typename T>
void prox_max(Image<T> &p, const Image<T> &y, const Image<T> &v, const Image<T> &t)
{
    // get the size
    int B = v.rows();
    int N = y.total();
    
    // scroll the image pixels
    #ifdef _OPENMP
    omp_set_num_threads(NUM_OF_THREADS);
    #pragma omp parallel for default(none) shared(p, y, v, t, N, B, T)
    #endif
    for(int i=0; i < N; i++) 
    {
        const T* v_col = v.ptr(0,i);
        const T* t_col = t.ptr(0,i);
        
        // project the block
        p(i) = prox_max( y(i), v_col, t_col, B );
    }
}