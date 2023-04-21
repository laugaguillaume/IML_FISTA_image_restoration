 function p = fun_abs_modif(x, gamma)
%function p = fun_abs(x, gamma)
%
% This procedure evaluates the function
%
%                       f(x) = gamma * |x|
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (08-07-2021)
% Author  : Nelly Pustelnik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check input
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
%-----%


% evaluate the function
p = gamma(:) .* abs(x(:)) ;