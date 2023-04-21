function [Y] = prox_Huber(X,tau,rho)
%PROX_HUBER Proximal operator of the standard Huber function
%   Detailed explanation goes her
[N,M] = size(X);
% N,M
Y = zeros(N,M);
% size(Y)
for k=1:N*M
    absX = abs(X(k));
    if absX > tau + rho
        Y(k) = (X(k)-tau*X(k)/absX) ;
    else
        Y(k) = rho*X(k)/(tau+rho);
    end
end

% 
% size(Y)

