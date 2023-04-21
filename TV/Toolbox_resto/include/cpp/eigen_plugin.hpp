#ifndef EIGEN_PLUGIN_HPP
#define EIGEN_PLUGIN_HPP


#include <Eigen/Core>
using Eigen::MatrixXd;
using Eigen::Map;
using Eigen;    // to be commented in case of conflicts

//
// NOTE: Eigen uses a column-major order, like matlab !!!
//

template <typename T>
void convert_to_eigen(MatrixXd &y, const Image<T> &x)
{
    // map the input
    Map<const MatrixXd> xx( x.ptr(0), x.rows(), x.cols() );
    
    // copy (the resizing is automatic)
    y = xx.cast<double>();
}

template <typename T>
void convert_from_eigen(Image<T> &y, const MatrixXd &x)
{
    if( x.rows() != y.rows() || x.cols() != y.cols() )
        mexErrMsgTxt("conversion error");
    
    // map the output
    Map<MatrixXd> yy( y.ptr(0), x.rows(), x.cols() );
    
    // copy (the resizing is NOT automatic)
    yy = x.cast<T>();
}

#endif