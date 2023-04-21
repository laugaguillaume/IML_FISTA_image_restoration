Nr = 256;
Nc = 256;
Nb = 3;

% random inputs
x = 256 * rand([Nr Nc Nb]);

% operators
f.dir_op = @(y) tv_dir_mex(y);
f.adj_op = @(y) tv_adj_mex(y);

% check adjoint
check_linop(f, [Nr Nc Nb]);

% speed test
tic; Fx = f.dir_op(x); toc;
tic; y = f.adj_op(Fx); toc;