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

#ifndef SVD_HPP
#define SVD_HPP

#include <mx_image.hpp>
#include <Eigen>
using namespace Eigen;

#ifdef _OPENMP
#include <omp.h>
#define NUM_OF_THREADS 4
#endif


template <typename T>
void sv_dec(Image<T> &u, Image<T> &s, Image<T> &v, const Image<T> &x)
{
    // get the size
    int N1 = x.size(0);
    int N2 = x.size(1);
    int B = std::min(N1,N2);
    int N = x.size(2)*x.size(3);
    
    // scroll the field
//     #ifdef _OPENMP
//     Eigen::initParallel();                  // disactivate OpenMP in Eigen
//     omp_set_num_threads(NUM_OF_THREADS);
//     #pragma omp parallel for default(none) shared(u, s, v, x, N1, N2, B, N, T)
//     #endif
    for(int b=0; b < N; b++) 
    {
        // pointers
        T* p_u = u.ptr(0,0,b);
        T* p_s = s.ptr(0,b);
        T* p_v = v.ptr(0,0,b);
        const T* p_x = x.ptr(0,0,b);
        
            // wrapping for Eigen
        Map<MatrixXd> uu(p_u, N1, B);
        Map<VectorXd> ss(p_s, B);
        Map<MatrixXd> vv(p_v, N2, B);
        Map<const MatrixXd> xx(p_x, N1, N2);
        
        // s.v. decomposition
        JacobiSVD<MatrixXd,HouseholderQRPreconditioner> svd(xx, ComputeThinU | ComputeThinV);
        
        // write the output
        uu = svd.matrixU();
        ss = svd.singularValues();
        vv = svd.matrixV();
    }
}



template <typename T>
void sv_rec(Image<T> &y, const Image<T> &u, const Image<T> &s, const Image<T> &v)
{
    // get the size
    int N1 = y.size(0);
    int N2 = y.size(1);
    int B = std::min(N1,N2);
    int N = y.size(2)*y.size(3);
    
    // scroll the field
//     #ifdef _OPENMP
//     Eigen::initParallel();                  // disactivate OpenMP in Eigen
//     omp_set_num_threads(NUM_OF_THREADS);
//     #pragma omp parallel for default(none) shared(y, u, s, v, N1, N2, B, N, T)
//     #endif
    for(int b=0; b < N; b++) 
    {
        // pointers
        T* p_y = y.ptr(0,0,b);
        const T* p_u = u.ptr(0,0,b);
        const T* p_s = s.ptr(0,b);
        const T* p_v = v.ptr(0,0,b);
        
        // wrapping for Eigen
        Map<MatrixXd> yy(p_y, N1, N2);
        Map<const MatrixXd> uu(p_u, N1, B);
        Map<const VectorXd> ss(p_s, B);
        Map<const MatrixXd> vv(p_v, N2, B);
    
        // s.v. reconstruction
        yy = uu * ss.asDiagonal() * vv.transpose();
    }
}


#endif