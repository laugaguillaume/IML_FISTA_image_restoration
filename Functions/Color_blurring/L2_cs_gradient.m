function [gradfx] = L2_inpainting_gradient(X,Z,A)
%L2_CS_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
[n,p,q] = size(X);
gradfx = zeros(n,p,q);
if q > 1
    for k=1:q
        gradfx(:,:,k) = A(:,:,k).*(A(:,:,k).*X(:,:,k)-Z(:,:,k));
    end
else
    gradfx = A.*(A.*X-Z);
end
end

