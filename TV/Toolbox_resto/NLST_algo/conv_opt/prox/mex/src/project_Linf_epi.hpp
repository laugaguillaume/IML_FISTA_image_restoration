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

#ifndef PROJECT_LINF_EPI_HPP
#define PROJECT_LINF_EPI_HPP

#include "mx_image.hpp"
#include <algorithm>
#include <vector>
#include <queue>


template <typename T>
struct couple
{
	T one, two;

    couple() {}
	couple(T one_, T two_) : one(one_), two(two_) {}
    
	bool operator<(const couple<T> &that) const {
		return one < that.one;
	}
};


template <typename T>
void proj_Linf_epi(T* p_col, T &t, const T* y_col, T xi, const T* w_col, int B)
{
    // compute the weighted abs
    std::vector< couple<T> > y_abs;
    y_abs.reserve(B);
    for(int r=0; r < B; r++) 
    {
        T w = w_col[r];
        if(w != 0) {
            T w_abs = w * std::abs( y_col[r] );
            couple<T> obj(w_abs, w);
            y_abs.push_back(obj);
        }
    }
    
    // build the heap
    std::priority_queue< couple<T> > heap( y_abs.begin(), y_abs.end() );
    int H = heap.size();
    
    // compute the scaling factor 'a'
    T a = T();
    if( xi > heap.top().one )
        a = xi;
    else
    {
        T sum = xi;
        T den = 1;
        for(int r = 0; r < H-1; r++)
        {
            // next value
            couple<T> v = heap.top();
            heap.pop();
            
            // scaling factor
            T t2 = 1 / (v.two * v.two);
            sum += v.one * t2;
            den += t2;
            T b = sum / den;
            
            // stop criterion
            T v_prev = heap.top().one;
            if(v_prev < b) {
                a = b;
                break;
            }
        }
        if(a == 0) {
            couple<T> v = heap.top();
            T t2 = 1 / (v.two * v.two);
            sum += v.one * t2;
            den += t2;
            a = sum / den;
        }
    }
    
    // compute the projection
    a = std::max(a, T());
    for(int r=0; r < B; r++)
    {
        T w = w_col[r];
        if( w == 0 )
            p_col[r] = y_col[r];
        else {
            T min = std::min( std::abs(y_col[r]), a / w );
            p_col[r] = (y_col[r]<0) ? -min : min;
        }
    }
    t = a;
}


template <typename T>
void proj_Linf_epi(T* p_col, T &t, const T* y_col, T xi, int B)
{
    // compute the abs
    std::vector<T> y_abs(B);
    for(int r=0; r < B; r++) {
        y_abs[r] = std::abs( y_col[r] );
    }
    
    // build the heap
    std::priority_queue<T> heap( y_abs.begin(), y_abs.end() );
    
    // compute the scaling factor 'a'
    T a = T();
    if( xi > heap.top() )
        a = xi;
    else
    {
        T sum = xi;
        for(int r = 2; r <= B; r++)
        {
            T v_curr = heap.top();
            heap.pop();
            T v_prev = heap.top();
            
            sum += v_curr;
            T b = sum / r;
            
            if(v_prev < b && b <= v_curr) {
                a = b;
                break;
            }
        }
        if(a==0) {
            a = (sum + heap.top()) / (B+1);
        }
    }
    
    // compute the projection
    a = std::max(a, T());
    for(int r=0; r < B; r++)
    {
        T min = std::min(y_abs[r], a);
        p_col[r] = (y_col[r]<0) ? -min : min;
    }
    t = a;
}


#endif