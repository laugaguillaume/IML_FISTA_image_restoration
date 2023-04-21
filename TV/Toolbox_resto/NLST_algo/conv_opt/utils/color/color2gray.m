 function Y = color2gray(X)
%function Y = color2gray(X)
%
%  Created on: 30/11/2011 - Giovanni Chierchia
%
%
% The function convert a RGB image to grayscale in the following way:
%
%     I(x,y) = 0.2989 * R(x,y) + 0.5870 * G(x,y) + 0.1140 * B(x,y)
%
% REMARK:
% This function produces exactly the same result as the standard MATLAB 
% function rgb2gray when a uint8/uint16 matrix is provides as input.
% In other words, this function assumes that also single/double matrices 
% have a dinamic into the range [0, 255].
%

if size(X,3) == 1
    
    Y = X;
    
else

    % compute the RGB  weights
    T = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);
    w = T(1,:);
    
    % convert to grayscale
    Y = w(1) * X(:,:,1) + w(2) * X(:,:,2) + w(3) * X(:,:,3);
end