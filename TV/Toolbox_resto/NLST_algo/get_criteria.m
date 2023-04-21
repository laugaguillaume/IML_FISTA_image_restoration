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

function [fid_f, reg_f, norm_f,  snr_f, epi_f] = get_criteria(data, deg, h, g, eta, reg_opt, x_inf, proj_mode)

if  strcmpi(deg.dataterm, 'l1')
    fid_f = @(y,xi) deg.noisecoef*sum(abs( h.dir_op(y) - deg.noisy(deg.valid) ));
elseif  strcmpi(deg.dataterm, 'l2')
    fid_f = @(y,xi) deg.noisecoef*norm( h.dir_op(y) - deg.noisy(deg.valid) )^2;
elseif  strcmpi(deg.dataterm, 'l2ball')
    fid_f = @(y,xi) max(norm( h.dir_op(y) - deg.noisy(deg.valid),'fro')-deg.noisecoef,0);
elseif  strcmpi(deg.dataterm, 'l1ball')
    fid_f = @(y,xi) max(sum(sum(abs( h.dir_op(y) - deg.noisy(deg.valid))))-deg.noisecoef,0);  
elseif  strcmpi(deg.dataterm, 'linfball')
    fid_f = @(y,xi) max(max( reshape(h.dir_op(y) - deg.noisy(deg.valid),1,numel(deg.noisy(deg.valid))))-deg.noisecoef,0);      
end

if  strcmpi(reg_opt.PenOrCons, 'penalisation')
    reg_f = @(y) reg_opt.p_eta*compute_seminorm(y, g.dir_op, reg_opt.norm);
elseif  strcmpi(reg_opt.PenOrCons, 'constraint')
    reg_f = @(y,xi) max( compute_seminorm(y, g.dir_op, reg_opt.norm) - reg_opt.p_eta, 0 );
end


if ~isempty(x_inf)
    norm_f = @(y,xi) norm( y(:) - x_inf(:) ) / norm( x_inf(:) );
else
    norm_f = @(y,xi) 0;
end

snr_f = @(y,xi) SNR( data, y );

if strcmpi(proj_mode, 'epi')
    epi_f = @(y,xi) abs( eta - sum(xi(:)) );
else
    epi_f = @(y,xi) 0;
end