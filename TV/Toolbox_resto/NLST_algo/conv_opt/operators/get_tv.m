 function f = get_tv(is_normalized)
%function f = get_tv(is_normalized)
%
%  Created on: 15/09/12 - Giovanni Chierchia
%
% The function creates the operator used for TV criterion. 

% default input
if nargin < 1
    is_normalized = 1;
end
%-----%

% operator norm
if is_normalized
    op_norm = sqrt(8);
    f.beta = 1;
else
    op_norm = 1;
    f.beta = 8;
end

% TV operator
f.dir_op = @(y) tv_dir_mex(y / op_norm);
f.adj_op = @(y) tv_adj_mex(y) / op_norm;

