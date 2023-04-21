 function p = fun_truncatedL1(x, lambda, alpha, dir)
%function p = fun_hyperbolic(x, gamma, w, dir)
%
% This procedure evaluates the function:
%
%                    f(x) = lambda^/2*x^2                     if x^2 \leq alpha/lambda^2 
%                         = alpha/2(2-alpha/(lambda^2 x^2))  otherwise
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x     - ND array
%  lambda - positive, scalar
%  beta  - positive, scalar
%  dir   - integer, direction of block-wise processing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (05-10-2021)
% Author  : Nelly Pustelnik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % default inputs
% if nargin < 3 || (~isempty(dir) && dir == 0)
%     dir = [];
% end
% 
% % check input
% sz = size(x); sz(dir) = 1;
% if any( gamma(:) <= 0 ) || ~isscalar(gamma) && (isempty(dir) || any(size(gamma)~=sz))
%     error('''gamma'' must be positive and either scalar or compatible with the blocks of ''x''')
% end
% %
% if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
%     error('''w'' must be positive and either scalar or the same size as ''x''')
% end
% %-----%


% linearize
if isempty(dir)
    x = x(:); 
end
    
% evaluate the function
ind = find(x.^2<alpha/lambda^2);
f = alpha/lambda^2*(2-alpha./(lambda^2*x.^2));
f(ind) = x(ind).^2;
p = sum( lambda^2/2*f(:) );
