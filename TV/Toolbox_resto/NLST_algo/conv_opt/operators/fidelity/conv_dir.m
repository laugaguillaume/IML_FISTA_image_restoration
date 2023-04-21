 function y = conv_dir(x, psf, border)
%function y = conv_dir(x, psf, border)
%
%  Created on: 13/09/2012 - Giovanni Chierchia
%
% The function computes the convolution of vectors 'x' and 'psf', producing
% an output vector of the same size as 'x'. By default, the borders are 
% zero-padded. An optional input string allows  to set a different border 
% handling: 'replicate', 'symmetric' or 'circular'.
%
% More precisely, this function implements the linear operation:
%
%                     y = (D*A*R) * x 
% where:
%  - R is the border extension operator
%  - A is the convolution operator
%  - D is the decimation operator (which makes the output vector of the same size as the input)
% 

% default input
if nargin < 3 || isempty(border)
    border = 0;
end

% direct operator
y = imfilter(x, psf, border);