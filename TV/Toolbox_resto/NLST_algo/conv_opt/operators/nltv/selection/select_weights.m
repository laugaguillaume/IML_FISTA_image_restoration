 function w = select_weights(w)
%function w = select_weights(w)
%
%  Created on: 02/11/2012 - Giovanni Chierchia
%


% weight selection
w = select_symmetric(w);
    
% reset central
idx = ceil( size(w,3) / 2 );
w(:,:,idx) = 0;

% weight normalization
avg = sum(w,3);
w = w ./ repmat( avg, [1 1 size(w,3)] );