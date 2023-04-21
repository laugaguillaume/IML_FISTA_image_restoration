function p = prox_huberl2(x,gamma,rho)

p = x - gamma.*x./abs(x);
ind = find(abs(x)<=  (gamma + rho));
p(ind) = rho* x(ind)/(gamma+rho);

    
