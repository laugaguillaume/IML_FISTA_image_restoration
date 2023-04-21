#ifndef PROJECT_SPECTRAL_EPI_HPP
#define PROJECT_SPECTRAL_EPI_HPP


#include "mx_image.hpp"
#include <Eigen>
using namespace Eigen;


template <typename T>
double find_thresh(const VectorXd &s, const T xi);


template <typename T>
void project_Spectral_epi(T* p_mat, T &t, const T* y_mat, const T xi, int N1, int N2, int B)
{        
    // wrapping for Eigen
    Map<      MatrixXd> pp(p_mat, N1, N2);
    Map<const MatrixXd> yy(y_mat, N1, N2);
    
    // s.v. decomposition
    JacobiSVD<MatrixXd,HouseholderQRPreconditioner> svd(yy, ComputeThinU | ComputeThinV);
    VectorXd ss = svd.singularValues();     // the sing. values are arranged in decreasing order
    
    // Linf projection
    T a = find_thresh( ss, xi );
    for(int r=0; r < B; r++)
        ss(r) = std::min(ss(r), a);
    t = a;
    
    // s.v. reconstruction
    pp = svd.matrixU() * ss.asDiagonal() * svd.matrixV().transpose();
}



template <typename T>
double find_thresh(const VectorXd &s, const T xi)
{
    T a = T();
    if( xi > s(0) )
        a = xi;
    else
    {
        T sum = xi;
        int B = s.size();
        for(int r = 0; r < B-1; r++)
        {
            T v_curr = s(r);
            T v_prev = s(r+1);
            
            sum += v_curr;
            T b = sum / (r+2);
            
            if(v_prev < b && b <= v_curr) {
                a = b;
            }
        }
        if(a == 0)
        	a = (sum + s(B-1)) / (B+1);
    }
    
    return std::max(a, T());
}

#endif