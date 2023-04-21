% (2015)
%
% Author :
% Giovanni Chierchia (chierchi@telecom-paristech.fr)
%
% Contributors :
% Nelly Pustelnik (nelly.pustelnik@ens-lyon.fr)
% Jean-Christophe Pesquet (jean-christophe.pesquet@univ-paris-est.fr)
% B?atrice Pesquet (beatrice.pesquet@telecom-paristech.fr)
%
% This software contains some image processing algorithms whose purpose is to be
% used primarily for research.
%
% This software is governed by the CeCILL B license under French law and
% abiding by the rules of distribution of free software. You can use,
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info".
%
% As a counterpart to the access to the source code and rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty and the software's author, the holder of the
% economic rights, and the successive licensors have only limited
% liability.
%
% In this respect, the user's attention is drawn to the risks associated
% with loading, using, modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean that it is complicated to manipulate, and that also
% therefore means that it is reserved for developers and experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or
% data to be ensured and, more generally, to use and operate it in the
% same conditions as regards security.
%
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL B license and that you accept its terms.
% January 2015 Giovanni Chierchia
% Non-local Total-Variation 
% for multicomponent (color, multispectral, hyperspectral,...)
% image restoration
%
% This toolbox implements the algorithm presented in the paper: 
% G. Chierchia, N. Pustelnik, B. Pesquet-Popescu, J.-C. Pesquet, "A Non-Local 
% Structure Tensor Based Approach for Multicomponent Image Recovery Problems",
% IEEE Trans. on Image Process., 2014
%
% Modified by N. Pustelnik
% July 8th 2021.

function [x_est, it, time, cost] = restore_TVh(data, deg, reg_mode, reg_opt, proj_mode, opt, x_inf, x_init, N_times)

% default params
if nargin < 9
    N_times = 1;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PROXIMITY OPERATORS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% range constraint
f.prox = @(y,gamma) project_range(y, [0 255]);

% smoothness constraint
[g, eta] = get_smoothness(data, x_init, reg_mode, reg_opt, proj_mode);
if strcmpi(reg_opt.norm, 'hyperbolic')
    g.grad =@(y) grad_hyperbolic();
    g.prox =@(y) prox_hyperbolic();
elseif strcmpi(reg_opt.norm, 'huber')  
    g.grad =@(y) grad_huber();
    g.prox =@(y) prox_huber();
else

% fidelity term
h.dir_op  = @(y) blur_dir(y, deg.psf, deg.valid, 0);
h.adj_op  = @(y) blur_adj(y, deg.psf, deg.valid, 0);
h.beta = compute_linop_norm( size(deg.noisy), h )^2;
%  
h.prox = @(y,gamma) deg.noisy(deg.valid) + prox_square( y - deg.noisy(deg.valid), gamma );  
%h.prox = @(y,gamma) deg.noisy(deg.valid) + prox_abs( y - deg.noisy(deg.valid), gamma );
%h.prox = @(y,gamma) deg.noisy(deg.valid) + prox_huber( y - deg.noisy(deg.valid), gamma, opt.epshuber );
%h.prox = @(y,gamma) prox_KLD( y, gamma, deg.noisy(deg.valid), 1);




%%%%%%%%%%%%%%%%%
%%% ALGORITHM %%%
%%%%%%%%%%%%%%%%%

% proj mode
if strcmpi(proj_mode, 'epi')
    algo = @(varargin) LFBF_epi(x_init, f, {g h}, [], opt, varargin{:});
else
    algo = @(varargin) LFBF_algo(x_init, f, [g h], [], opt, varargin{:});
end

for iter=1:N_times
    
    % criteria
    [fid_f, reg_f, norm_f, snr_f, epi_f] = get_criteria(data, deg, h, g, eta, reg_opt, x_inf, proj_mode);

    % algorithm
    [x_est, it, time, cost.fid, cost.reg, cost.norm, cost.snr, cost.epi] = algo(fid_f, reg_f, norm_f, snr_f, epi_f);
    
    x_inf = x_est;
end

if N_times > 1
    figure;
    plot(time, cost.norm, 'linewidth', 2);
end