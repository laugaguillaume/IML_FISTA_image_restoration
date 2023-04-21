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

#include "src/project_Nuclear_epi.hpp"

#ifdef _OPENMP
#include <omp.h>
#define NUM_OF_THREADS 6
#endif


template <typename T>
void proj_epi_Nuclear_mex(Image<T> &p, Image<T> &t, const Image<T> &y, const Image<T> &xi);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check the params
    if (nrhs != 2)
        mexErrMsgTxt("This function takes 2 inputs");
    if (nlhs != 2)
        mexErrMsgTxt("This function gives 2 output");
    
    // pixel type
    typedef double T;
    
    // read the inputs
    MxImage<T> y ( prhs[0] );
    MxImage<T> xi( prhs[1] );
    
    // check the inputs
    int N1 = y.size(0);
    int N2 = y.size(1);
    int Nr = y.size(2);
    int Nc = y.size(3);
    int B = N1;//std::min(N1,N2);
    if( N1 > N2 )
        mexErrMsgTxt("This function is designed for 'fat' matrices.");
    if( B*Nr*Nc != xi.total() )
        mexErrMsgTxt("The inputs are not compatible");
    
    // init. the output
    MxImage<T> p(plhs[0], N1, N2, Nr, Nc);
    MxImage<T> t(plhs[1], B, Nr, Nc);
    
    // write the output
    proj_epi_Nuclear_mex(p, t, y, xi);
}


template <typename T>
void proj_epi_Nuclear_mex(Image<T> &p, Image<T> &t, const Image<T> &y, const Image<T> &xi)
{
    // get the size
    int N1 = y.size(0);
    int N2 = y.size(1);
    int B = N1;//std::min(N1,N2);
    int N = y.size(2)*y.size(3);
 
    // scroll the image pixels
    #ifdef _OPENMP
    Eigen::initParallel();                  // disactivate OPEN_MP in Eigen
    omp_set_num_threads(NUM_OF_THREADS);
    #pragma omp parallel for default(none) shared(p, t, y, xi, N1, N2, B, N, T)
    #endif
    for(int c=0; c < N; c++) 
    {
              T*  p_mat =  p.ptr(0,0,c);
              T*  t_mat =  t.ptr(0,c);
        const T*  y_mat =  y.ptr(0,0,c);
        const T* xi_mat = xi.ptr(0,c);
        
        // project the block
        project_Nuclear_epi(p_mat, t_mat, y_mat, xi_mat, N1, N2, B);
    }
}