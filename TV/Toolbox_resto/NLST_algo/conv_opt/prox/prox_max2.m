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

function p = prox_max2(y, v, t)
%function p = prox_max2(y, v, t)
%
%  Created on: 23/03/2013 - Giovanni Chierchia
%
% The function computes the proximity operator of the function:
%
%   h(y) = 0.5 * \sum_i \max{0, t(i) * (v(i) - y)}^2
%
% where t(i) >= 0. The function works in a block-wise manner. It expects 
% the blocks to be stored on the 1st dimension of v and t. A separate 
% projection is computed on each block.

% BUGFIX:
% if the elements in v(i,:) are "almost" all equal, replace them with their average
v_orig  = sqrt(sum(v.^2,1));
v_round = sqrt(sum(round(v).^2,1));
j = abs(v_orig-v_round) < 10^-10 * abs(v_orig);
v(:,j) = repmat( mean( v(:,j), 1 ), [size(v,1) 1]);

% compute the projection
if nargin < 3 || isempty(t) || all( t(:) == 1 )
    p = prox_max_matlab(y, v);
else
    p = prox_max_mex(y, v, t);
end






 function p = prox_max_matlab(y, v)

% cumpute the candidates for p
v = sort(v, 1);
s = cumsum( flipdim(v,1) , 1 );
yy = shiftdim(y, -1);
s = bsxfun(@plus, yy, s);
s = bsxfun(@rdivide, [flipdim(s,1); yy], 2 + size(v,1) - (1:size(v,1)+1)' );

% locate the "right" candidate
infty = inf([1 size(v,2)]);
v_low  = [-infty; v];
v_high = [v; +infty];
mask = v_low < s & s <= v_high;

% select the candidate
p = reshape( s(mask), size(y) );