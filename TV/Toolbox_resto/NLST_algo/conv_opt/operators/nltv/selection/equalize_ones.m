function [mask Nb_new] = equalize_ones(mask)
% 
% The function adds more 1s in the input mask, in order to get the same
% amounts of 1s along the 3rd dimension.
%

% current neigh. size
Nb = size(mask,3);

% new neigh. size
N_all = sum(mask,3);
Nb_new = max( N_all(:) );

% matrix of elements to add
N = Nb_new - N_all;

% mask of non-selected elements
sel_mask = ~mask;                           
sel_mask(:,:,ceil(Nb/2),:) = false;           % never select the central pixel
    
% mask of first N non-sel. elem.
cnt = cumsum(sel_mask,3);
cnt_mask = cnt <= repmat(N, [1 1 Nb 1]);

% update the mask
mask = mask | (sel_mask & cnt_mask);