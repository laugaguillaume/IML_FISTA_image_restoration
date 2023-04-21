function [normfX] = blur_L2_inpainting_function(X,Z,A1,A2,A3)
%BLUR_L2_COLOR_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
[n,p,q] = size(X);
if q == 3
    X(:,:,1) = reshape(A1*reshape(X(:,:,1),n*p),n,p);
    X(:,:,2) = reshape(A2*reshape(X(:,:,2),n*p),n,p);
    X(:,:,3) = reshape(A3*reshape(X(:,:,3),n*p),n,p);
    
    normfX = 0.5*norm(X(:)-Z(:))^2;
elseif q == 1
    normfX = 0.5.*norm(A1*X(:)-Z(:))^2;
end

end

