function [v] = update_v_CV(g, x, v, gamma)

Ng = length(g);
t = cell(1,Ng);

for m = 1:Ng
    v_hat   = v{m}  + gamma * g(m).dir_op(x);
    t{m}    = v_hat - gamma * g(m).prox(v_hat/gamma, 1/gamma);
end