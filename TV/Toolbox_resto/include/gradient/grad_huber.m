 function p = grad_huber(x, gamma, w)
%function p = grad_huber(x, gamma, w)
%
% This procedure performs the gradient of the function:
%
%           / gamma * x^2 / 2            if |x| <= w
%   f(x) = |                                                with w > 0
%           \ gamma * (w * |x| - w^2/2)  otherwise
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  w     - positive, scalar or ND array with the same size as 'x'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (08-07-2021)
% Author  : Nelly Pustelnik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check input
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
    error('''w'' must be positive and either scalar or the same size as ''x''')
end
%-----%


% evaluate the "square" branch
p = gamma .* x;

% evaluate the "abs" branch
t = gamma .* w .* sign(x);

% merge the branches
mask = abs(x) > w;
p(mask) = t(mask);

