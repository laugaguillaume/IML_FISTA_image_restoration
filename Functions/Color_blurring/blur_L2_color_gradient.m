function [gradfx] = blur_L2_color_gradient(X,Z,Ac,Ar,Acolor)
%BLUR_L2_COLOR_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
[n,p,q] = size(X);
gradfx = zeros(n,p,q);
if q == 3
    % Forward operator
    AX1 = Ac*X(:,:,1)*Ar';
    AX2 = Ac*X(:,:,2)*Ar';
    AX3 = Ac*X(:,:,3)*Ar';
    X(:,:,1) = Acolor(1)*AX1 + Acolor(2)*AX2 + Acolor(3)*AX3 ;
    X(:,:,2) = Acolor(4)*AX1 + Acolor(5)*AX2 + Acolor(6)*AX3 ;
    X(:,:,3) = Acolor(7)*AX1 + Acolor(8)*AX2 + Acolor(9)*AX3 ;
    % Adjoint operator
    gradfx1 = Ac'* (X(:,:,1) -Z(:,:,1) )*Ar;
    gradfx2 = Ac'* (X(:,:,2) -Z(:,:,2) )*Ar;
    gradfx3 = Ac'* (X(:,:,3) -Z(:,:,3) )*Ar;
    % 3D gradient
    gradfx(:,:,1) = Acolor(1)*gradfx1 + Acolor(4)*gradfx2 + Acolor(7)*gradfx3 ;
    gradfx(:,:,2) = Acolor(2)*gradfx1 + Acolor(5)*gradfx2 + Acolor(8)*gradfx3 ;
    gradfx(:,:,3) = Acolor(3)*gradfx1 + Acolor(6)*gradfx2 + Acolor(9)*gradfx3 ;
elseif q == 1
    gradfx = Ac'*(Ac*X*Ar'-Z)*Ar;
end

end

