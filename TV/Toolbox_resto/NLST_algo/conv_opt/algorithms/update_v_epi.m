function [v1, v2, t1, t2] = update_v_epi(g, x, x2, p, p2, v1, v2, gamma)

v_hat  = v1 + gamma * g.dir_op(x);
v2_hat = v2 + gamma * x2;

[epi_1, epi_2] = g.project_epi(v_hat/gamma, v2_hat/gamma);
t1  = v_hat - gamma * epi_1;
t2 = v2_hat - gamma * epi_2;

v_tilde  = t1 + gamma * g.dir_op(p);
v2_tilde = t2 + gamma * p2;

v1 = v1 - v_hat  + v_tilde;
v2 = v2 - v2_hat + v2_tilde;