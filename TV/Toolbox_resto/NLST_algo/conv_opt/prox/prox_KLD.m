function prox = prox_KLD(x,gamma,y, coeff)

% Created by N. Pustelnik
% July 8th 2021

prox = (x - coeff*gamma + sqrt(abs(x-coeff*gamma).^2 + 4*gamma*y))/2;