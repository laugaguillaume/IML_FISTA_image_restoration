function [gradfx] = L2_inpainting_gradient(X,Z,A,Ab)
%L2_CS_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
[n,p,q] = size(X);
gradfx = zeros(n,p,q);
if q > 1
    for k=1:q
        gradfx(:,:,k) = (A(:,:,k).*reshape(Ab*reshape(X(:,:,k),n*p,1),n,p)-Z(:,:,k));
        gradfx(:,:,k) = A(:,:,k).*reshape(Ab'*reshape(gradfx(:,:,k),n*p,1),n,p);
    end
else
    gradfx = (A.*reshape(Ab*reshape(X,n*p,1),n,p)-Z)
    gradfx = A.*reshape(Ab'*reshape(gradfx,n*p,1),n,p);
end
end

