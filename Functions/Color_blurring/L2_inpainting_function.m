function [normfX] = L2_inpainting_function(X,Z,A1,A2,A3)
%BLUR_L2_COLOR_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
[n,p,q] = size(X);
if q == 3
    X(:,:,1) = A1.*X(:,:,1);
    X(:,:,2) = A2.*X(:,:,2);
    X(:,:,3) = A3.*X(:,:,3);
    normfX = 0.5*norm(X(:)-Z(:))^2;
elseif q == 1
    normfX = 0.5.*norm(reshape(A1.*X-Z,n*p,1))^2;
end

end

