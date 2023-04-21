Nr = 256;
Nc = 256;
Nb = 3;
B = 16;

% random inputs
x = 256 * rand([Nr Nc Nb]);
idx = int32( randi([1 Nr*Nc], [B Nr*Nc]) -1);
w = rand(size(idx));

% operators
f.dir_op = @(y) nltv_dir_mex(y, w, idx);
f.adj_op = @(y) nltv_adj_mex(y, w, idx);

% check adjoint
check_linop(f, [Nr Nc Nb]);

% speed test
tic; Fx = f.dir_op(x); toc;
tic; y = f.adj_op(Fx); toc;