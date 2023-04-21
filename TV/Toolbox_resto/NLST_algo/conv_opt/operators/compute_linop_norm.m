 function [nxe cnt] = compute_linop_norm(sz, h, tol, maxiter)
%function [nxe cnt] = compute_linop_norm(sz, h, tol, maxiter)
%
%  Created on: 14/09/2012 - Giovanni Chierchia
%
% The function compute the norm of the linear operator F, where the latter 
% is defined by its direct and adjoint functions.
%
% COPYRIGHT:
% This function is a modification of LINOP_NORMEST in TFOCS toolbox.
%
% TFOCS v1.2 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2012 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.


% default inputs
if nargin < 3 || isempty( tol ),
    tol = 1e-8;
end
if nargin < 4 || isempty( maxiter ),
    maxiter = 100;  % should never take this many iterations. 
end
%-----%


% iterative estimation
cnt = 0;
nxe = 0;
while true,
    if nxe == 0,
        xx = randn(sz);
        nxe = norm( xx(:) );
    end
    yy = h.dir_op( xx / max(nxe,realmin) );
    nye = norm( yy(:) );
    xx = h.adj_op( yy / max(nye,realmin) );
    nxe0 = nxe;
    nxe = norm( xx(:) );
    if abs(nxe - nxe0) < tol * max(nxe0, nxe),
        break;
    end
    cnt = cnt + 1;
    if cnt >= maxiter,
        break;
    end
end





%%%%%%%%%%%%%%%%
%%% OLD CODE %%%
%%%%%%%%%%%%%%%%

 function cteDM = compute_beta(sz, f)
%function cteDM = compute_beta(sz, f)
%
%  Created on: --/--/--, Nelly Pustelnik
%
% The function compute the square norm ||F||^2 of the linear operator F,
% where the latter is defined by its direct and adjoint functions.

xn = rand(sz);

cteDM = 0;
rho_n  = 1 + 1e-6; 
rho_n1 = 1;
while abs(rho_n1 - rho_n) > 1e-8 * rho_n
    rho_n = rho_n1;
    xn1 = f.adj_op( f.dir_op(xn) );
    rho_n1 = norm( xn1(:) ) / norm( xn(:) );
    xn = xn1;
end

if ~isfinite(rho_n1)        %%% NOTA: AGGIUNTO QUESTO IF PER MIGLIORARE LA STABILITA !!!
    rho_n1 = rho_n;
end

if rho_n1 > cteDM
    cteDM = rho_n1;
end
