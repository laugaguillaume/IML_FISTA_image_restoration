function x_hat = compute_adjoint(g, v)

x_hat = 0;
for m = 1:length(g)
    x_hat = x_hat + g(m).adj_op(v{m});
end