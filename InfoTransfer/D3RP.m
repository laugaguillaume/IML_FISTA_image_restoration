function [XH] = D3RP(X,R)
%3DRP Summary of this function goes here
%   R = R if fine to coarse
%   R = P if coarse to fine
[n,p] = size(R');
[~,~,q] = size(X);
XH = zeros(p,p,q);
if q == 3
    XH(:,:,1) = R*X(:,:,1)*R';
    XH(:,:,2) = R*X(:,:,2)*R';
    XH(:,:,3) = R*X(:,:,3)*R';
elseif q == 1 
    XH = R*X*R';
elseif q>3
    for k=1:q
        XH(:,:,k) = R*X(:,:,k)*R';
    end
end
end

