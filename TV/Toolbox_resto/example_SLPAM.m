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

clear all; %close all;
addpath( genpath('NLST_algo'), '-end' );
addpath(genpath('build'))
addpath(genpath('/Users/nellypustelnik/Dropbox/Travaux_en_cours/20XX_TOOLBOX/Optim/prox_repository/matlab'))
addpath(genpath('include/'))
addpath(genpath('vmlmb-master/'))




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETER SELECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image name
image_list = {'cameraman'};   

% noise params
deg.noise_sigma     = 10;
deg.blur            = 1;
deg.decimation_rate = 0;
deg.psf             = fspecial('average', deg.blur);
deg.border          = 'zeropad';
deg.noise           = 'Gaussian';
deg.noiseFixed      = 1;

% regularization params (for NLST)
%nlst.norm_list = {'hyperbolique'};  % options: 'l1' 'l2' 'linf' 'nuclear' 'frobenius' 'spectral' 'huber' 'hyperbolique'
%nlst.tv_eta    = 100;
%nlst.tv_beta   = 6;
%nlst.tv_lamda   = 0.015;
%nlst.nltv_eta  = 0.5;
%nlst.neigh_rad = 2;
%nlst.blk_rad   = 2;
%nlst.delta     = 35;

% projection method
proj_mode = 'epi';        

% algorithmic parameters
opt.iter = 1000;                  % number of iterations
opt.tol  = 10^-12;                 % stopping criterion (you may set 10^-4 for faster, but less accurate, results)
opt.approxl1 = 1e-4;
opt.dir = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUN THE SIMULATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx_images = 1 : length(image_list)

    % print info
    fprintf('\n %s \n', image_list{idx_images}); for mm = 1:length(image_list{idx_images})+2, fprintf('='); end, fprintf('\n');
         
    % read the image
    x_true = get_image( image_list{idx_images} );      
    
    % apply noise
    deg = degradation_blur_noise(x_true, deg);
    
    % perform restoration
    %for norm_idx = 1:length(nlst.norm_list)
        

        x_init = deg.noisy;
        
        %SLPAM
        opt.lambda= 0.6;
        opt. alpha = 50;
        
        opt.lam = opt.alpha;
        opt.bet = opt.lambda^2/2;
        opt.edges = 'distinct';
        nlst.norm_list = {'l1'};
        reg_opt.norm  = nlst.norm_list{1};
        %reg_opt.p_eta = nlst.tv_eta;     
        [x_est_SLPAM, it_SLPAM,  time_SLPAM, stepsize_SLPAM,cost_SLPAM] = restore_SLPAM(deg, 'TV', 'SLPAM', reg_opt, opt, [], x_init);

        
        % TVhyperbolique EA with BT
        %nlst.norm_list = {'hyperbolic'};
        %reg_opt.norm  = nlst.norm_list{1};
        %nlst.tv_eta    = 100;
        %reg_opt.p_eta = nlst.tv_eta;   
        opt.lam = opt.alpha;
        opt.bet = opt.lambda;   
        

        [x_est_EAwBT, it_EAwBT,  time_EAwBT, stepsize_EAwBT, cost_EAwBT] = restore_TruncatedTV(deg, 'TV', 'EAwBT', reg_opt, opt, [], x_init);

%        fprintf('TV: (SNR) = (%2.2f) --> eta: %2.2f, it: %d, time: %3.2f\n', SNR(x_true,x_est_EA), reg_opt.p_eta, it_EA, time_EA(end));
        fprintf('SLPAM: (SNR) = (%2.2f) --> lambda-beta: %2.2f - %2.2f, it: %d, time: %3.2f\n', SNR(x_true,x_est_SLPAM{1}), opt.lam, opt.bet, it_SLPAM, time_SLPAM(end));
        fprintf('TV: (SNR) = (%2.2f) --> lambda-beta: %2.2f - %2.2f, it: %d, time: %3.2f\n', SNR(x_true,x_est_EAwBT{1}), opt.lam, opt.bet, it_EAwBT, time_EAwBT(end));
            

    %end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    PRINT RESULTS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure(1);
 subplot(221);imagesc(resize_print(x_true)); axis image off; title 'Original'
 subplot(222);imagesc(resize_print(deg.noisy)); axis image off;  title 'Degraded'
 subplot(224);imshow(resize_print(x_est_EAwBT{1})); axis image off; title 'Dual with restored EAwBT'
 hold on;
plot_contours_v2((x_est_EAwBT{2}));
 subplot(223);imshow(resize_print(x_est_SLPAM{1})); axis image off; title 'Primal with SLPAM'
hold on;
plot_contours_v2((x_est_SLPAM{2}));
print -depsc result3_images

% subplot(236);imagesc(resize_print(x_est_VMLMB)); axis image off; title 'TV restored VMLMB'
% 
figure(2)
semilogy(cost_SLPAM,'k');
hold on
semilogy(cost_EAwBT.fid+cost_EAwBT.reg,'r');
%subplot(122);plot(cost_EAwBT,'r');
legend('SLPAM','EAwBT');
grid on;
print -depsc result3_cost
