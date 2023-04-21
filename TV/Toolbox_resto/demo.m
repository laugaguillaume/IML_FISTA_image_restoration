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
% Modified by N. Pustelnik
% July 8th 2021

clear all; close all;
addpath(genpath('NLST_algo'), '-end' );
addpath(genpath('build'))
addpath(genpath('images'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     LOAD DATA       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

load data_blur.mat               % Load data and blur
deg.psf   = blur ;
deg.noisy = data ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETER SELECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% regularization params (for NLST)
nlst.norm_list = {'nuclear'};  % options: 'l1' 'l2' 'linf' 'nuclear' 'frobenius' 'spectral'
nlst.tv_eta    = 0.5;          % regularization parameter TV
nlst.nltv_eta  = 0.5;          % regularization parameter NLTV
nlst.neigh_rad = 2;            % radius of the non-local support 
nlst.blk_rad   = 2;            % radius of the patch for self-similarity measure
nlst.delta     = 35;           % parameter involved in the self-similarity measure

% projection method
proj_mode = 'epi';        

% algorithmic parameters
opt.iter = 1000;                  % number of iterations
opt.tol  = 10^-5;                 % stopping criterion (you may set 10^-4 for faster, but less accurate, results)

% no decimation
deg.valid = true( size(deg.noisy) );
deg.decimation_rate = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUN THE SIMULATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st minimization
reg_opt.norm             = nlst.norm_list;
reg_opt.p_eta            = nlst.tv_eta;
[x_est1 it1 time1 cost1] = restore_NLST(deg.noisy, deg, 'TV', reg_opt, proj_mode, opt, [], deg.noisy);


% 2nd minimization
if ~isempty(nlst.nltv_eta)

    reg_opt.norm      = nlst.norm_list;
    reg_opt.p_eta     = nlst.nltv_eta;
    reg_opt.neigh_rad = nlst.neigh_rad;
    reg_opt.blk_rad   = nlst.blk_rad;
    reg_opt.delta     = nlst.delta;

    [x_est2 it2 time2 cost2] = restore_NLST(x_est1, deg, 'NL-TV', reg_opt, proj_mode, opt, [], deg.noisy);

end
       


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    PRINT RESULTS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(131);imagesc(resize_print(deg.noisy)); axis image off;  title 'Degraded'
subplot(132);imagesc(resize_print(x_est1)); axis image off; title 'TV restored'
subplot(133);imagesc(resize_print(x_est2)); axis image off; title 'NLTV restored'

