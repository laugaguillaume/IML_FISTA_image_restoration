function [gradfx] = L2_inpainting_gradient(X,Z,A1,A2,A3)
%BLUR_L2_COLOR_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
[n,p,q] = size(X);
gradfx = zeros(n,p,q);
if q == 3
    X1 = X(:,:,1); Z1 = Z(:,:,1);
    X2 = X(:,:,2); Z2 = Z(:,:,2);
    X3 = X(:,:,3); Z3 = Z(:,:,3);
    % Forward operator
    AX1 = A1.*X1;
    AX2 = A2.*X2;
    AX3 = A3.*X3;
    % Adjoint operator
    gradfx1 = A1.* (X1 -Z1);
    gradfx2 = A2.* (X2 -Z2);
    gradfx3 = A3.* (X3 -Z3);
    % 3D gradient
    gradfx(:,:,1) = gradfx1 ;
    gradfx(:,:,2) = gradfx2;
    gradfx(:,:,3) = gradfx3 ;
elseif q == 1
    gradfx = A1.*(A1.*X-Z);
end

end

