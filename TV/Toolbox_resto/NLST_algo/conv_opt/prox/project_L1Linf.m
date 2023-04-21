% (2015)
%
% Author :
% Giovanni Chierchia (chierchi@telecom-paristech.fr)
%
% Contributors :
% Nelly Pustelnik (nelly.pustelnik@ens-lyon.fr)
% Jean-Christophe Pesquet (jean-christophe.pesquet@univ-paris-est.fr)
% Béatrice Pesquet (beatrice.pesquet@telecom-paristech.fr)
%
% This software contains some image processing algorithms whose purpose is to be
% used primarily for research.
%
% This software is governed by the CeCILL B license under French law and
% abiding by the rules of distribution of free software. You can use,
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info".
%
% As a counterpart to the access to the source code and rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty and the software's author, the holder of the
% economic rights, and the successive licensors have only limited
% liability.
%
% In this respect, the user's attention is drawn to the risks associated
% with loading, using, modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean that it is complicated to manipulate, and that also
% therefore means that it is reserved for developers and experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or
% data to be ensured and, more generally, to use and operate it in the
% same conditions as regards security.
%
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL B license and that you accept its terms.
% January 2015 Giovanni Chierchia
% Non-local Total-Variation 
% for multicomponent (color, multispectral, hyperspectral,...)
% image restoration
%
% This toolbox implements the algorithm presented in the paper: 
% G. Chierchia, N. Pustelnik, B. Pesquet-Popescu, J.-C. Pesquet, "A Non-Local 
% Structure Tensor Based Approach for Multicomponent Image Recovery Problems",
% IEEE Trans. on Image Process., 2014

function p = project_L1Linf(x, eta, w)
%function p = project_L1Linf(x, eta, w)
%
%  Created on: 01/07/2012 - Giovanni Chierchia
%
%
% This function performs the projection onto the L_{1,\infty}-norm ball:
%
% C = { y \in \RR^{L} : \sum_{\ell=1}^L ||y^{(\ell}||_\infty \le \eta }
% 
% The function expects the input to be a P*B*N*M 4-D matrix, and it works
% by considering that L = B*N*M is the number of blocks and P is the length
% of a single block. 
%
% NOTE:
% When used on the output matrix of tv_dir or nltv_dir, this function 
% projects onto the "channel-by-channel" Linf-TV constraint:
%
%                \sum_{i=1}^B  Linf-TV(x_i) \le \eta
%
% where x = (x_1,...,x_B) is a multicomponent N*M image.



% check the inputs
[b n] = size(x);
if nargin < 3 || isempty(w)
    w = ones([n 1]);
end
if numel(w) ~= n
    error('The size of w is not compatible with the size of x');
end
%-------------------%


% get the size
B = size(x,1);
N = size(x,2)*size(x,3)*size(x,4);

% linearize
x_lin = reshape(x, [B N])';

% project
p_lin = projL1Inf( x_lin, eta, w(:) );

% restore
p = reshape(p_lin', size(x));