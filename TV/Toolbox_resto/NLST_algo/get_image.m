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

function [x best_eta] = get_image(name)


switch name
           
    case 'firemen'
        best_eta.tv = 0.5;
        best_eta.nltv = 0.5;
        
        x = imread('images/firefighters.jpg');
        
    case 'sim'
        best_eta.tv = 0.5;
        best_eta.nltv = 0.5;
        load('sim.mat');
        x = s;
                                                       
    case 'culicoidae'
        best_eta.tv   = 0.8;
        best_eta.nltv = 0.73;
        
        x = imread('images/culicoidae.jpg');
        dim = 256; x = x(381:380+dim, 921:920+dim);
        
    case 'dot'
        best_eta.tv   = 0.8;
        best_eta.nltv = 0.73;
        
        x = imread('images/dots-256.png');
         
        
    case 'indian'
        best_eta.tv   = 0.4;
        best_eta.nltv = 0.3;
        
        % load image
        x = imread('images/indian.tif');
        x(:,:,[104:108 150:164]) = [];
        
    case 'cameraman'
        best_eta.tv   = 0.4;
        best_eta.nltv = 0.3;
        
        % load image
        x = imread('images/cameraman.tif');
        %x(:,:,[104:108 150:164]) = [];        
    
    case {'paris.lan' 'rio.lan' 'mississippi.lan' 'montana.lan' 'littlecoriver.lan' 'tokyo.lan'}
        best_eta.tv   = 0;
        best_eta.nltv = 0;
        
        x = multibandread(['images/' name '.lan'], [512, 512, 7], 'uint8=>double', 128, 'bil', 'ieee-le');
    
    otherwise
        best_eta.tv   = 0;
        best_eta.nltv = 0;
        
        x = imread(im_name);
end
x = double(x);

% set dynamics range to [0 255]
stretch_fun = @(y) 255 * ( y - min(y(:)) ) / ( max(y(:)) - min(y(:)) );
for idx_band = 1:size(x,3)
    x(:,:,idx_band) = stretch_fun( x(:,:,idx_band) );
end