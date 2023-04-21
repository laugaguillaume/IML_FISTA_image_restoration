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

function [p t] = proj_epi_Linf(y, xi, w)
%function [p t] = proj_epi_Linf(y, xi, w)
%
%  Created on: 17/06/11 - Giovanni Chierchia
%
%
% The function computes the projection onto the epigraph
%
%   epi_Linf = { (y,xi) \in \RR^B \times \RR  :  ||w*y||_\infty \le xi },
%
% The function works in a block-wise manner. It expects the blocks to be 
% stored on the 1st dimension of y. A separate projection is computed on 
% each block.


% check inputs
if nargin < 3 || isempty(w)
    w = 1;
end
[B N] = size(y);
if N ~= numel(xi) 
    error('Inputs not compatible');
end
if any( w(:) <= 0 )
    error('The weights must be positive.');
end
%-----%


% project t
y_abs = w .* abs(y);
t = prox_max2(xi, y_abs, 1./w);
t = max(t, 0);

% replicate t
sz = size(y);
tt = reshape( t(:), [1 sz(2:end)] );
tt = repmat(tt, [sz(1) 1]);

% project p
p = sign(y) .* min( abs(y), tt ./ w );