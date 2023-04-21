function [normfX] = blur_L2_color_function(X,Z,Ac,Ar,Acolor)
%BLUR_L2_COLOR_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
[n,p,q] = size(X);
if q == 3
    AX1 = Ac*X(:,:,1)*Ar';
    AX2 = Ac*X(:,:,2)*Ar';
    AX3 = Ac*X(:,:,3)*Ar';
    
    X(:,:,1) = Acolor(1)*AX1 + Acolor(2)*AX2 + Acolor(3)*AX3 ;
    X(:,:,2) = Acolor(4)*AX1 + Acolor(5)*AX2 + Acolor(6)*AX3 ;
    X(:,:,3) = Acolor(7)*AX1 + Acolor(8)*AX2 + Acolor(9)*AX3 ;
    
    normfX = 0.5*norm(X(:)-Z(:))^2;
elseif q == 1
    normfX =0.5.*norm(reshape((Ac*X*Ar'-Z),n*p,1))^2;
end
end

