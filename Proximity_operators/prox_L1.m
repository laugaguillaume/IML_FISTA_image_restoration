function [Y] = prox_L1(X,tau)
Y = sign(X).*max(abs(X)-tau,0);

