function [U] = GetWcoef(W,j,orient) 

% [U] = GetWcoef(W,j,orient) 
% Extract the wavelet coefficients of an image decomposition
% Inputs
%   W : wavelet coefficient (e.g. obtained from FWT2_PO)
%   j : scale of the wavelet coefficients to be extracted (0<j<log2(size
%       W))
%   orient : 
%       'approx' : approximation coefficients
%       'det-H'  : detail coefficients along horizontal direction
%       'det-V'  : detail coefficients along vertical direction 
%       'det-D'  : detail coefficients along diagonal directionfor 
%       'det-ALL': detail coefficients along all 3 directions (H, V and D)

% PG : nov 2021

orient = lower(orient) ;
J = log2(min(size(W))) ;
if j >= J
    error(['extraction scale j must be smaller than ',num2str(J),...
        ' and larger than the tree depth decomposition (W)'])
end

switch orient
    case 'approx'
        U = W(1:2^(j),1:2^(j)) ;
    case 'det-v'
        U = W(2^(j)+1:2^(j+1),1:2^j) ;
    case 'det-h'
        U = W(1:2^j,2^(j)+1:2^(j+1)) ;
    case 'det-d'
        U = W(2^(j)+1:2^(j+1),2^(j)+1:2^(j+1)) ;
    case 'det-all'
        U = [zeros(2^j,2^j) GetWcoef(W,j,'det-H') ; ...
            GetWcoef(W,j,'det-V') GetWcoef(W,j,'det-D')] ;
end