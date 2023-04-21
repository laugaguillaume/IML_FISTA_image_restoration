 function p = fun_hyperbolic_modif(x, gamma, w, dir)
%function p = fun_hyperbolic(x, gamma, w, dir)
%
% This procedure evaluates the function:
%
%                    f(x) = gamma * sqrt( ||x||_2^2 + w)
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array compatible with the blocks of 'x'
%  w     - positive, scalar or ND array with the same size as 'x'
%  dir   - integer, direction of block-wise processing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (08-07-2021)
% Author  : Nelly Pustelnik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% default inputs
if nargin < 3 || (~isempty(dir) && dir == 0)
    dir = [];
end

% check input
sz = size(x); sz(dir) = 1;
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && (isempty(dir) || any(size(gamma)~=sz))
    error('''gamma'' must be positive and either scalar or compatible with the blocks of ''x''')
end
%
if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
    error('''w'' must be positive and either scalar or the same size as ''x''')
end
%-----%


% linearize
if isempty(dir)
    x = x(:); 
end
    
% evaluate the function
xx = gamma .* sqrt( sum(x.^2, dir) + w );
p = xx;
%p = sum( xx(:) );