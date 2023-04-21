 function p = prox_huber_v2(x, gamma, w)
%function p = grad_huber_v2(x, gamma, w)
%
% This procedure computes the proximity operator of the function:
%
%           / gamma * x^2 / (2*w)            if |x| <= w
%   f(x) = |                                                with w > 0
%           \ gamma * (|x| - w/2)  otherwise
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


% compute the "square" branch
p = w.*x ./ (gamma + w);

% compute the "abs" branch
t = x - gamma .* sign(x);

% select the branches
mask = abs(x) > gamma+w;
p(mask) = t(mask);
