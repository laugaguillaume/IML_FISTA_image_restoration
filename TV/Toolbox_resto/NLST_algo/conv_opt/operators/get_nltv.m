 function [f w] = get_nltv(w)
%function [f w] = get_nltv(w)
%
%  Created on: 01/05/12 - Giovanni Chierchia
%
% The function creates the operators for NLTV criterion. The input 'w' is a
% 3D-matrix of positive weights computed on all (2*rad+1)^2 neighbours.


% select & normalize
w = select_weights(w);

% image size
sz = [size(w,1) size(w,2)];

% weight reduction and reshaping
[w idx] = squeeze_weights(w);

% NLTV operators
f.dir_op = @(y) nltv_dir_mex(y, w, idx);
f.adj_op = @(y) nltv_adj_mex(y, w, idx);

% operator norm
opnorm = compute_linop_norm(sz, f);
f.beta = opnorm^2;

% normalization    
w = w / opnorm;
f.dir_op = @(y) nltv_dir_mex(y, w, idx);
f.adj_op = @(y) nltv_adj_mex(y, w, idx);
f.beta = 1;