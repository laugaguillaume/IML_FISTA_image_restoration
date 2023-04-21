 function y = conv_adj(x, psf, border)
%function y = conv_adj(x, psf, border)
%
%  Created on: 13/09/2012 - Giovanni Chierchia
%
% The function computes the adjoint convolution of vectors 'x' and 'psf', 
% producing an output vector of the same size as 'x'. By default, the 
% borders are zero-padded. An optional input string allows one to set 
% a different border handling: 'replicate', 'symmetric' or 'circular'.
%
% More precisely, this function implements the linear operation:
%
%                     y = (D*A*R)' * x 
% where:
%  - R is the border extension operator
%  - A is the convolution operator
%  - D is the decimation operator (which makes the output vector of the same size as the input)
% 

% default input
if nargin < 3 || isempty(border)
    border = 0;
end

% correct "even" sizes
sz = size(psf);
psf = padarray(psf, 1-mod(sz,2), 'pre');

% compute the adjoint operator
if isscalar(border) || strcmpi(border, 'circular')  % extensions: zero padding, periodic. 
    
    y = imfilter(x, rot90(psf,2), border);
    
else                                                % extensions: symmetric, replicate.
    
    % get the radius
    rad = (size(psf) - 1) / 2;
    
    % adjoint decimation
    xx = padarray(x, rad);
    
    % adjoint (linear) convolution
    y = imfilter(xx, rot90(psf,2));
    
    % adjoint replication
    y = replicate_adj(y, rad, border);
    
end





function y = replicate_adj(y, rad, border)

if strcmpi(border, 'replicate')
    
    y(rad(1)+1,:) = y(rad(1)+1,:) + sum( y(1:rad(1),:), 1 );
    y(1:rad(1),:) = [];
    
    y(end-rad(1),:) = y(end-rad(1),:) + sum( y(end-rad(1)+1:end,:), 1 );
    y(end-rad(1)+1:end,:) = [];
    
    y(:,rad(2)+1) = y(:,rad(2)+1) + sum( y(:,1:rad(2)), 2 );
    y(:,1:rad(2)) = [];
    
    y(:,end-rad(2)) = y(:,end-rad(2)) + sum( y(:,end-rad(2)+1:end), 2 );
    y(:,end-rad(2)+1:end) = [];
    
elseif strcmpi(border, 'symmetric')
    
    y(rad(1)+1:2*rad(1),:) = y(rad(1)+1:2*rad(1),:) + y(rad(1):-1:1,:);
    y(1:rad,:) = [];
    
    y(end-2*rad(1)+1:end-rad(1),:) = y(end-2*rad(1)+1:end-rad(1),:) + y(end:-1:end-rad(1)+1,:);
    y(end-rad(1)+1:end,:) = [];
    
    y(:,rad(2)+1:2*rad(2)) = y(:,rad(2)+1:2*rad(2)) + y(:,rad(2):-1:1);
    y(:,1:rad(2)) = [];
    
    y(:,end-2*rad(2)+1:end-rad(2)) = y(:,end-2*rad(2)+1:end-rad(2)) + y(:,end:-1:end-rad(2)+1);
    y(:,end-rad(2)+1:end) = [];
    
end