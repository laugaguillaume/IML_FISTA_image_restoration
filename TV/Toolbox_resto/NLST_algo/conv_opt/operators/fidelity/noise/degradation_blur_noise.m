 function deg = degradation_blur_noise(data, deg)
%function deg = degradation_blur_noise(data, deg)
%
%  Created on: 17/09/2012 - Giovanni Chierchia
%
% The function performs the 3 following operations:
%  1. blurring
%  2. noise adding
%  3. pixel decimation



%%% PRELIMINARES %%%

% Fix the noise
if deg.noiseFixed == 1
    rand('seed',10);
    randn('seed',10);
end

% add a zero if the 'psf' size is even
if isfield(deg, 'psf')
    sz = size(deg.psf);
    deg.psf = padarray(deg.psf, 1-mod(sz,2), 'pre');
end

% get positions of missing pixels
if (deg.decimation_rate~=0)
    deg.valid = get_random_indices( size(data), 1-deg.decimation_rate );
else
    deg.valid = true( size(data) );
end

% check border extension
if strcmpi(deg.border, 'zeropad')
    ext = 0;
elseif strcmpi(deg.border, 'circular')
    ext = 'circular';
else
    error('Border extension must be ''zeropad'' or ''circular''');
end



%%% BODY %%%

% blur the data
xc = imfilter(data, deg.psf, ext);

% add the noise
if strcmpi(deg.noise,'Gaussian')
    
    z = xc + deg.noise_sigma .* randn( size(xc) );
    
elseif strcmpi(deg.noise,'Poisson')
    
    z = random('poisson', deg.noise_sigma .* xc);
    
elseif strcmpi(deg.noise,'Laplace')
    
    z = xc + get_laplace_noise(size(xc), 0, deg.noise_sigma);
    
elseif strcmpi(deg.noise,'Speckle')
    
    L = deg.noise_sigma;
    z = random('Gamma', L, xc/L);   % NOTE: this is equivalent to z = x .* gamrnd(L, 1/L, size(x));
    
else
    error('Undefined noise\n');
end

% decimate the data
deg.noisy = z .* deg.valid;