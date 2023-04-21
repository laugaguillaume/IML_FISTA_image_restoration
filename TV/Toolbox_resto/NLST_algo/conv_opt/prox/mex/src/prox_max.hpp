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

#ifndef PROX_MAX_HPP
#define PROX_MAX_HPP

#include <algorithm>
#include <queue>



template <typename T>
T prox_max(T y, const T* v, int B)
{    
    // build the heap
    std::priority_queue<T> heap(v, v+B);
    
    // compute the projection
    T p = T();
    if( y > heap.top() )
        p = y;
    else
    {
        T sum = y;
        for(int r = 2; r <= B; r++)
        {
            // next value
            T v_curr = heap.top();
            heap.pop();
            T v_prev = heap.top();
            
            sum += v_curr;
            T b = sum / r;
            
            if(v_prev < b) {
                p = b;
                break;
            }
        }
        if(p==0) {
            p = (sum + heap.top()) / (B+1);
        }
    }
    
    return p;
}


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
T prox_max(T y, const T* v, const T* t, int B)
{    
    // store the weights
    std::vector< couple<T> > vt(B);
    for(int i=0; i < B; i++) {
        vt[i] = couple<T>( v[i], t[i] );
    }
    
    // build the heap
    std::priority_queue< couple<T> > heap( vt.begin(), vt.end() );
    
    // compute the projection
    T p = T();
    if( y > heap.top().one )
        p = y;
    else
    {
        // initialize
        T sum = y;
        T den = 1;
        
        // scroll in decreasing order
        for(int i = 0; i < B-1; i++)
        {
            // get current/next values
            couple<T> v = heap.top();
            heap.pop();
            T v_prev = heap.top().one;
            
            // accumulate the current value/weight
            T t2 = v.two * v.two;
            sum += v.one * t2;
            den += t2;
            
            // compute the candidate
            T b = sum / den;
            
            // check the stopping criterion
            if(v_prev < b) {
                p = b;
                break;
            }
        }
        
        // if it didn't stop, the projection needs the last element
        if(p == 0) 
        {
            // accumulate the last element
            couple<T> v = heap.top();
            //T t2 = 1 / (v.two * v.two);
            T t2 = v.two * v.two;
            sum += v.one * t2;
            den += t2;
            
            // compute the projected value
            p = sum / den;
        }
    }
    
    return p;
}


#endif