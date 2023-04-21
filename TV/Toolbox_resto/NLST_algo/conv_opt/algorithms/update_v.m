function [v t] = update_v(g, x, p, v, gamma)

Ng = length(g);
t = cell(1,Ng);

for m = 1:Ng
    v_hat   = v{m}  + gamma * g(m).dir_op(x);
    t{m}    = v_hat - gamma * g(m).prox(v_hat/gamma, 1/gamma);
    v_tilde = t{m}  + gamma * g(m).dir_op(p);
    v{m}    = v{m}  - v_hat + v_tilde;
end