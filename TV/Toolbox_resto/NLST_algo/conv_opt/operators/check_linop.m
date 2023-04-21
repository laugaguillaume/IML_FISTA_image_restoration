 function check_linop(f, sz)
%function check_linop(f, sz)
%
%  Created on: 17/09/2012 - Giovanni Chierchia
%
%
% This function checks if the adjoint operator is correct with respect 
% to his respective direct operator. To this end, the code employs the
% definition of adjoint: <x,Ay> = <A'x,y>, where A denotes a linear 
% operator and A' its adjoint.
%


if nargin < 2
    sz = 500;
end

% draw (randomly) a matrix for the direct op.
x = rand(sz);

% compute the direct operator
x_dir = f.dir_op(x);

% draw (randomly) a matrix for the adjoint op.
sz = size(x_dir);
y = rand(sz);

% compute the adjoint operator
y_adj = f.adj_op(y);

% check the correctness
a1 = y(:)' * x_dir(:);
a2 = x(:)' * y_adj(:);
d = abs(a1-a2);

% print
fprintf('|<x,Ay> - <A''x,y>| = %e\n', d);