 function [w idx] = squeeze_weights(w, opt)
%function [w idx] = squeeze_weights(w, opt)
%
%  Created on: 02/11/2012 - Giovanni Chierchia
%
% This functions remove the zeros and reorder the weights colum-wise.

% default input
if nargin < 2
    opt = 'zeropad';
end
%-----%


% mask of selected pixels
mask = (w ~= 0);

% equalization of the mask
mask = equalize_ones(mask);

% index computation (for C++/mex operators, so "0-based" indices)
[idx mask] = get_indices(mask, opt);
idx = int32(idx-1);

% weight reduction
[Nr Nc Nb Np] = size(w);
w = permute(w, [3 1 2 4]);
w = reshape(w, [Nb Nr*Nc Np]);
w = reshape( w(mask), size(idx) );