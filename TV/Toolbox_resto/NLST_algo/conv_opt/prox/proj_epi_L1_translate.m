% (2015)
%
% Author :
% Giovanni Chierchia (chierchi@telecom-paristech.fr)
%
% Contributors :
% Nelly Pustelnik (nelly.pustelnik@ens-lyon.fr)
% Jean-Christophe Pesquet (jean-christophe.pesquet@univ-paris-est.fr)
% B?atrice Pesquet (beatrice.pesquet@telecom-paristech.fr)
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

function [p, t] = proj_epi_L1_translate(y, xi, d, w)
%function [p t] = proj_epi_L1(y, xi, w)
%
%  Created on: 17/06/11 - Giovanni Chierchia
%
%
% The function computes the projection onto the epigraph
%
%         epi_1 = { (y,xi) \in \RR x \RR : w * |y| \le xi },
%
% giving two outputs: (p,t) = Proj(y,xi). The function is vectorialized. 
% When the inputs are vectors, it computes the projection of each element.


if nargin < 3 || isempty(w)
    w = 1;
end
if any(w(:) <= 0)
    error('The weights must be greater than zero');
end

y = y-d;

% compute the projection
y_abs = abs(y);
c = max(0, y_abs + w .* xi) ./ (1 + w.^2);
p = c .* sign(y);
t = c .* w;

% correction 
mask = w .* y_abs < xi;
p(mask) =  y(mask);
t(mask) = xi(mask);

p = p+d;