N = 9;
y = 100 * randn(N,N);
eta = 50*rand;


%%% direct

y_lin = reshape(y, [sqrt(N) sqrt(N) N]);
tic; p_dir = project_Linf(y_lin, eta); toc;
p_dir = reshape(p_dir, [N N]);


%%% Epigraph

f.prox = @(x, gamma) x;

h.grad = @(x) x - y;
h.dir_op  = @(y) y;
h.adj_op  = @(y) y;
h.beta    = 1;

g.project_epi = @(x, xi) project_Linf_epi(x, xi);
g.project_V = @(x) project_V(x, eta);
g.r = sqrt([N N]);
g.dir_op = @(x) reshape(x, [sqrt(N) sqrt(N) N]);
g.adj_op = @(x) reshape(x, [N N]);
g.beta = 1;

opt.iter = 100000;
opt.tol  = eps;
fun = @(a,b) 0;
tic; [p_epi it] = LFBF_epi(y, f, g, h, opt, fun); toc;

% test
L1inf = @(x) sum( max( abs(x), [], 2 ) );
fprintf( 'DIRECT - ||p||_inf - eta: %2.2f, ||P(y) - y||: %2.2f\n', L1inf(p_dir) - eta, norm(p_dir(:) - y(:)) );
fprintf( 'EPIGR. - ||p||_inf - eta: %2.2f, ||P(y) - y||: %2.2f\n', L1inf(p_epi) - eta, norm(p_epi(:) - y(:)) );