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

function [semi_norm Fx] = compute_seminorm(x, dir_op, norm_type, w)


% default inputs
if nargin < 4 || isempty(w)
    w = 1;
end
%-----%


% direct operator
Fx = w .* dir_op(x);

% compute the TV amount
if strcmpi(norm_type, 'L1')

    % compute the element-wise norms
    Fx = abs(Fx);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
elseif strcmpi(norm_type, 'L2')

    % compute the block-wise norms
    Fx = Fx.^2;
    Fx = sum(Fx, 1);
    Fx = sqrt(Fx);
    Fx = shiftdim(Fx,1);
    
    % sum the norms
    semi_norm = sum( Fx(:) );

        
elseif strcmpi(norm_type, 'L12')

    % compute the block-wise norms
    Fx = Fx.^2;
    Fx = sum(Fx, 1);
    Fx = sqrt(Fx);
    Fx = shiftdim(Fx,1);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
elseif strcmpi(norm_type, 'inf')
    
    % compute the block-wise norms
    Fx = abs(Fx);
    Fx = max(Fx, [], 1);
    Fx = shiftdim(Fx,1);
    
    % sum the norms
    semi_norm = sum( Fx(:) );

elseif strcmpi(norm_type, 'nuclear')
    
    % compute the block-wise norms
    [u s v] = sv_dec(Fx);
    Fx = abs(s);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
elseif strcmpi(norm_type, 'frobenius')
   
    % compute the block-wise norms
    Fx = Fx.^2;
    Fx = sum(Fx, 1);
    Fx = sum(Fx, 2);
    Fx = sqrt(Fx);
    Fx = squeeze(Fx);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
elseif strcmpi(norm_type, 'spectral')
    
    % compute the block-wise norm
    [u s v] = sv_dec(Fx);
    Fx = max(s, [], 1);
    Fx = squeeze(Fx);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
else
    error('Norm not recognized.') 
end