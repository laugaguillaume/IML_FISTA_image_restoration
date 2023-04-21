 function y = get_laplace_noise(sz, mu, sigma)
%function y = get_laplace_noise(sz, mu, sigma)
%
%  Created on: 22/09/2011, Giovanni Chierchia
%
%
% The function generates 'sz' samples drawn from the laplacian pdf:
%
%                 f(x) = (1/2*b) * exp(-|x-mu| / b)
%
% where:
%  -   mean   = mu
%  - variance = 2 * b^2


% default inputs
if nargin < 1
    error('Too few inputs');
end
if nargin < 2 || isempty(mu)
    mu = 0;
end
if nargin < 3 || isempty(mu)
    sigma = 1;
end

% check the inputs
if sigma <= 0
    error('The dev. std. must be a positive value');
end

% generate the vector
x = rand(sz);

% compute the masks
low  = (x <= 0.5);
high = (x >  0.5);

% generate the laplacian noise (with mean 0 and variance 1)
y = x;
y(low)  =   log( 2 * x(low)  );
y(high) = - log( 1 - x(high) );

% enforce mean and variance
y = mu + (sigma/sqrt(2)) * y;