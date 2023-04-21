function [normfX] = L2_cs_function(X,Z,A)
%L2_CS_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
[n,p,q] = size(X);
if q > 1
    for k=1:q
        X(:,:,k) = A(:,:,k).*X(:,:,k);
    end
else
    X = A.*X;
end
normfX = 0.5*norm(X(:)-Z(:))^2;

end

