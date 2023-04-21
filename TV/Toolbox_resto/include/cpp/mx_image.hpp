#ifndef MX_MEX_HPP
#define MX_MEX_HPP

#include "image.hpp"
#include "mex_traits.hpp"

template <typename T>
struct MxImage : public Image<T>
{

    MxImage(const mxArray* mx) {
        super(mx);
    }

    MxImage(mxArray* &mx, int N_rows, int N_cols, int N_bands = 1, int N_planes = 1) {
        const mwSize sz[4] = {N_rows, N_cols, N_bands, N_planes};
        mx = mxCreateNumericArray(4, sz, mxClassIDTrait<T>::type, mxREAL);
        super(mx);
    }

private:

    void super(const mxArray *mx)
    {        
        // complex checking
        if( mxIsComplex(mx) )
            mexErrMsgTxt("Complex values not supported");
        
        // type checking
        if( mxClassIDTrait<T>::type != mxGetClassID(mx) )
            mexErrMsgTxt("The template type does not match the matrix typed");
            
        // size checking
        mwSize Ndim = mxGetNumberOfDimensions(mx);
        if(Ndim > 4)
            mexErrMsgTxt("The array cannot have more than 4 dimensions");

        // data pointer
        T* data = (T*) mxGetData(mx);

        // matrix size
        const mwSize* sz  = mxGetDimensions(mx);
        int size[4];
        size[0] = (Ndim>0) ? sz[0] : 1;
        size[1] = (Ndim>1) ? sz[1] : 1;
        size[2] = (Ndim>2) ? sz[2] : 1;
        size[3] = (Ndim>3) ? sz[3] : 1;

        // roi handling
        size_t step[3];
        step[0] = size[0];
        step[1] = size[0]*size[1];
        step[2] = size[0]*size[1]*size[2];
        
        // configure the base class
        Image<T>::init(data, size, step);
    }
};


template <>
struct MxImage<mxLogical> : public Image<mxLogical>
{

    MxImage(const mxArray* mx) {
        super(mx);
    }

    MxImage(mxArray* &mx, int N_rows, int N_cols, int N_bands = 1, int N_planes = 1) {
        const mwSize sz[4] = {N_rows, N_cols, N_bands, N_planes};
        mx = mxCreateLogicalArray(4, sz);
        super(mx);
    }

private:

    void super(const mxArray *mx)
    {                
        // type checking
        if( mxLOGICAL_CLASS != mxGetClassID(mx) )
            mexErrMsgTxt("The matrix is not logical");
            
        // size checking
        mwSize Ndim = mxGetNumberOfDimensions(mx);
        if(Ndim > 4)
            mexErrMsgTxt("The array cannot have more than 4 dimensions");

        // data pointer
        mxLogical* data = mxGetLogicals(mx);

        // matrix size
        const mwSize* sz  = mxGetDimensions(mx);
        int size[4];
        size[0] = (Ndim>0) ? sz[0] : 1;
        size[1] = (Ndim>1) ? sz[1] : 1;
        size[2] = (Ndim>2) ? sz[2] : 1;
        size[3] = (Ndim>3) ? sz[3] : 1;

        // roi handling
        size_t step[3];
        step[0] = size[0];
        step[1] = size[0]*size[1];
        step[2] = size[0]*size[1]*size[2];
        
        // configure the base class
        Image<mxLogical>::init(data, size, step);
    }
};


#endif