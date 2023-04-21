function [Y] = prox_L21(X,tau)
%PROX_L21 Prox of the L21 mixed norm
%   Detailed explanation goes here
%   Input X is a matrix
Y = prox_L2(X,tau,1); % process along colums
end

