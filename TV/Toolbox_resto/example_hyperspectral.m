% (2015)
%
% Author :
% Giovanni Chierchia (chierchi@telecom-paristech.fr)
%
% Contributors :
% Nelly Pustelnik (nelly.pustelnik@ens-lyon.fr)
% Jean-Christophe Pesquet (jean-christophe.pesquet@univ-paris-est.fr)
% Béatrice Pesquet (beatrice.pesquet@telecom-paristech.fr)
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

clear all; close all;
addpath( genpath('NLST_algo'), '-end' );
addpath(genpath('build'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETER SELECTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image name
image_list = {'indian'};   

% noise params
deg.noise_sigma     = 5;
deg.blur            = 1;
deg.decimation_rate = 0.9;
deg.psf             = fspecial('average', deg.blur);
deg.border          = 'zeropad';
deg.noise           = 'Gaussian';
deg.noiseFixed      = 1;

% regularization params (for NLST)
nlst.norm_list = {'nuclear'};  % options: 'l1' 'l2' 'linf' 'nuclear' 'frobenius' 'spectral'
nlst.tv_eta    = 0.4;
nlst.nltv_eta  = 0.3;
nlst.neigh_rad = 2;
nlst.blk_rad   = 2;
nlst.delta     = 35;

% projection method
proj_mode = 'epi';        

% algorithmic parameters
opt.iter = 1000;                  % number of iterations
opt.tol  = 10^-5;                 % stopping criterion (you may set 10^-4 for faster, but less accurate, results)



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
    for norm_idx = 1:length(nlst.norm_list)
        disp(''); disp('------'); disp(['Norm: ' nlst.norm_list{norm_idx}]); disp('');
        
        % 1st minimization
        reg_opt.norm  = nlst.norm_list{norm_idx};
        reg_opt.p_eta = nlst.tv_eta;
        
        [x_est1 it1 time1 cost1] = restore_NLST(deg.noisy, deg, 'TV', reg_opt, proj_mode, opt, [], x_true);

        fprintf('TV: (SNR,M-SNR) = (%2.2f, %2.2f) --> eta: %2.2f, it: %d, time: %3.2f\n', SNR(x_true,x_est1), SNR_v2(x_true,x_est1), reg_opt.p_eta, it1, time1(end));
    
        % 2nd minimization
        if ~isempty(nlst.nltv_eta)

            reg_opt.norm      = nlst.norm_list{norm_idx};
            reg_opt.p_eta     = nlst.nltv_eta;
            reg_opt.neigh_rad = nlst.neigh_rad;
            reg_opt.blk_rad   = nlst.blk_rad;
            reg_opt.delta     = nlst.delta;
                        
            [x_est2 it2 time2 cost2] = restore_NLST(x_est1, deg, 'NL-TV', reg_opt, proj_mode, opt, [], x_true);
            
            fprintf('NLTV: (SNR,M-SNR) = (%2.2f, %2.2f) --> it: %d, time: %3.2f, reg: (%2.2f, %d, %d, %2.2f)\n', SNR(x_true,x_est2), SNR_v2(x_true,x_est2), it2, time2(end), reg_opt.p_eta, reg_opt.neigh_rad, reg_opt.blk_rad, reg_opt.delta);
                     
        end
       
        disp('------');
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    PRINT RESULTS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
bande1 = 100;
bande2 = bande1+2;
subplot(221);imagesc(resize_print(x_true(:,:,bande1:bande2))); axis image off; title 'Original'
subplot(222);imagesc(resize_print(deg.noisy(:,:,bande1:bande2))); axis image off;  title 'Degraded'
subplot(223);imagesc(resize_print(x_est1(:,:,bande1:bande2))); axis image off; title 'TV restored'
subplot(224);imagesc(resize_print(x_est2(:,:,bande1:bande2))); axis image off; title 'NLTV restored'

