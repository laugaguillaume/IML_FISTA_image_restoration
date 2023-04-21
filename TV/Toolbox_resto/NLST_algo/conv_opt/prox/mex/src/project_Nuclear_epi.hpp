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

#ifndef PROJECT_NUCLEAR_EPI_HPP
#define PROJECT_NUCLEAR_EPI_HPP


#include "mx_image.hpp"
#include <Eigen>
using namespace Eigen;

template <typename T>
void project_Nuclear_epi(T* p_mat, T* t_mat, const T* y_mat, const T* xi_mat, int N1, int N2, int B)
{        
    /* ATTENTION: this code is optimized for "fat" matrices, hence B = N1 */
    B = N1;
    
    // wrapping for Eigen
    Map<      MatrixXd> pp( p_mat, N1, N2);
    Map<      VectorXd> tt( t_mat, B);
    Map<const MatrixXd> yy( y_mat, N1, N2);
    Map<const VectorXd> xx(xi_mat, B);
    
    // s.v. decomposition
    MatrixXd yt = yy.transpose();
    SelfAdjointEigenSolver<MatrixXd> eig(yy * yt);
    
    // compute matrix U and vector S (singular values)
    MatrixXd uu = eig.eigenvectors();
    
    // compute vector S (singular values)
    VectorXd ss = eig.eigenvalues();
    ss = (ss.array() > 0).select(ss.cwiseSqrt(), 0); // avoid negative sqrt (due to numerical errors) !!!
    
    // invert elements of S
    VectorXd s_inv = (ss.array() > 1e-13).select(ss.cwiseInverse(), 0); // don't invert zeros!!!
    
    // compute matrix V
    MatrixXd vv = yt * uu;
    vv *= s_inv.asDiagonal();
    
    // compute the epi. proj. onto the L1-norm
    Array<bool,Dynamic,Dynamic> mask = ss.array() < xx.array();
    VectorXd gg = (ss + xx) / 2;
    tt = gg.cwiseMax(0);
    tt = mask.select(xx, tt);   // t_i = (s_i < x_i) ? x_i : max(0, s_i + x_i) / 2
    ss = mask.select(ss, tt);   // s_i = (s_i < x_i) ? s_i : max(0, s_i + x_i) / 2
    
//     for(int i=0; i<B; i++)
//     {
//         tt(i) = xx(i);
//         if( ss(i) > xx(i) )
//         {
//             ss(i) = std::max( T(0), ss(i) + xx(i) ) / 2;
//             tt(i) = ss(i);
//         }
//     }
    
    // s.v. reconstruction
    pp = uu * ss.asDiagonal() * vv.transpose();
}

#endif