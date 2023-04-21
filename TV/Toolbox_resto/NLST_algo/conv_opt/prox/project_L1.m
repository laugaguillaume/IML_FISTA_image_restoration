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

function [y itn] = project_L1(x, eta, w)
%function [y itn] = project_L1(x, eta, w)
%
% The function computes the projection onto the L1-norm ball:
%
%            C = { x \in R^n : sum(|x_i|) < tau }
%
% The toolbox SPGL_1 is used for the projection.


if nargin < 3 || isempty(w)
    w = 1;
end

%[y itn] = oneProjector( x(:), w(:), eta );      % this is a mex function (decomment for speedup)
[y itn] = oneProjector_matlab( x(:), w(:), eta );
y = reshape( y, size(x) );





%----------%
% We give below a matlab implementation of the "oneProjector" mexfile.
%----------%
function [x,itn] = oneProjector_matlab(b,d,tau)
% ONEPROJECTOR  Projects b onto the weighted one-norm ball of radius tau
%
%    [X,ITN] = ONEPROJECTOR(B,TAU) returns the orthogonal projection
%    of the vector b onto the one-norm ball of radius tau. The return
%    vector X which solves the problem
%
%            minimize  ||b-x||_2  st  ||x||_1 <= tau.
%               x
%
%    [X,ITN] = ONEPROJECTOR(B,D,TAU) returns the orthogonal
%    projection of the vector b onto the weighted one-norm ball of
%    radius tau, which solves the problem
%
%            minimize  ||b-x||_2  st  || Dx ||_1 <= tau.
%               x
%
%    If D is empty, all weights are set to one, i.e., D = I.
%
%    In both cases, the return value ITN given the number of elements
%    of B that were thresholded.
%
% See also spgl1.
%
%--------------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected Gradient for L1).
%--------------------------------------------------------------------------


% Check arguments
if nargin < 2
  error('The oneProjector function requires at least two parameters');
end
if nargin < 3
  tau = d;
  d   = [];
end

% Check weight vector
if isempty(d), d = 1; end;

if ~isscalar(d) && ( length(b) ~= length(d) )
  error('Vectors b and d must have the same length');
end

% Quick return for the easy cases.
if isscalar(d)  &&  d == 0
   x   = b;
   itn = 0;
   return
end

% Get sign of b and set to absolute values
s = sign(b);
b = abs(b);

% Perform the projection
if isscalar(d)
  [x,itn] = oneProjectorMex_matlab(b,tau/d);
else
  d   = abs(d);
  idx = find(d > eps); % Get index of all non-zero entries of d
  x   = b;             % Ensure x_i = b_i for all i not in index set idx
  [x(idx),itn] = oneProjectorMex_matlab(b(idx),d(idx),tau);
end

% Restore signs in x
x = x.*s;




% ----------------------------------------------------------------------
function [x, itn] = oneProjectorMex_matlab(b,d,tau)
% ----------------------------------------------------------------------
%function [x, itn] = oneProjectorMex(b,d,tau) 
% Return the orthogonal projection of the vector b >=0 onto the
% (weighted) L1 ball. In case vector d is specified, matrix D is
% defined as diag(d), otherwise the identity matrix is used.
%
% On exit,
% x      solves   minimize  ||b-x||_2  st  ||Dx||_1 <= tau.
% itn    is the number of elements of b that were thresholded.
%
% See also spgl1, oneProjector.

if nargin < 3
    tau = d;
    d   = 1;
end

if isscalar(d)
    [x,itn] = oneProjectorMex_I(b,tau/abs(d));
else
    [x,itn] = oneProjectorMex_D(b,d,tau);
end




% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_I(b,tau)
% ----------------------------------------------------------------------

% Initialization
n     = length(b);
x     = zeros(n,1);
bNorm = norm(b,1);

% Check for quick exit.
if (tau >= bNorm), x = b; itn = 0; return; end
if (tau <  eps  ),        itn = 0; return; end

% Preprocessing (b is assumed to be >= 0)
[b,idx] = sort(b,'descend'); % Descending.

csb       = -tau;
alphaPrev = 0;
for j= 1:n
    csb       = csb + b(j);
    alpha     = csb / j;
    
    % We are done as soon as the constraint can be satisfied
    % without exceeding the current minimum value of b
    if alpha >= b(j)
        break;
    end
    
    alphaPrev = alpha;
end

% Set the solution by applying soft-thresholding with
% the previous value of alpha
x(idx) = max(0,b - alphaPrev);

% Set number of iterations
itn = j;




% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_D(b,d,tau)
% ----------------------------------------------------------------------

% Initialization
n = length(b);
x = zeros(n,1);

% Check for quick exit.
if (tau >= norm(d.*b,1)), x = b; itn = 0; return; end
if (tau <  eps         ),        itn = 0; return; end

% Preprocessing (b is assumed to be >= 0)
[bd,idx] = sort(b ./ d,'descend'); % Descending.
b  = b(idx);
d  = d(idx);

% Optimize
csdb = 0; csd2 = 0;
soft = 0; alpha1 = 0; i = 1;
while (i <= n)
    csdb = csdb + d(i).*b(i);
    csd2 = csd2 + d(i).*d(i);
    
    alpha1 = (csdb - tau) / csd2;
    alpha2 = bd(i);
    
    if alpha1 >= alpha2
        break;
    end
    
    soft = alpha1;  i = i + 1;
end
x(idx(1:i-1)) = b(1:i-1) - d(1:i-1) * max(0,soft);

% Set number of iterations
itn = i;