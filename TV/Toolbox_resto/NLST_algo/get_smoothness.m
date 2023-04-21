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

function [f, eta] = get_smoothness(data, x_init, reg_mode, reg_opt, proj_mode)

if strcmpi(reg_mode, 'TV')
    
    f = get_tv;
    
elseif strcmpi(reg_mode, 'NL-TV')
    
    % compute the "guidance" image (single band)
    x_ref = mean(x_init,3);
    
    % estimate the weights    
    w = get_fov_weights(x_ref, reg_opt.delta, reg_opt.blk_rad, reg_opt.neigh_rad);
    
    % define the linear operator
    f = get_nltv(w);
        
end

% configure the constraint
if  strcmpi(reg_opt.PenOrCons, 'penalisation')
   f.prox = @(y,gamma) prox_L2(y, reg_opt.p_eta*gamma,1);
   eta = reg_opt.p_eta;
else strcmpi(reg_opt.PenOrCons, 'constraint')
    [f, eta] = get_constraint(f, data, reg_opt.norm, reg_opt.p_eta, proj_mode);
end
%elseif strcmpi(reg_opt.PenOrCons, 'cons')
%   f.prox = @(y,gamma)  project_L12(y, reg_opt.p_eta,1);
%end
