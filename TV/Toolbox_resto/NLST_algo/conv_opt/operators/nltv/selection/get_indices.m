 function [ind mask] = get_indices(mask, opt)
%function [ind mask] = get_indices(mask, opt)
%
% The function computes the linear indices of the selected pixels.


% default input
if nargin < 2
    opt = 'zeropad';
end
%-----%


% get the radius
[Nr Nc Nb Np] = size(mask);
rad = (sqrt(Nb)-1) / 2;

% equalize the mask
[mask Ns] = equalize_ones(mask);

% reorder the mask
mask = permute(mask, [3 1 2 4]);
mask = reshape(mask, [Nb Nr*Nc Np]);

% indices of non-zero elements
idx = find(mask);
idx = mod(idx-1, Nb)+1;
idx = reshape(idx, [Ns Nr*Nc Np]);
[r c] = ind2sub(2*[rad rad]+1, idx);
r = r - (rad+1);
c = c - (rad+1);

% matrix grid
[CC RR] = meshgrid(1:Nc, 1:Nr);
CC = repmat( CC(:)', [Ns 1 Np]);
RR = repmat( RR(:)', [Ns 1 Np]);

% compute the indices
rr = RR + r;
cc = CC + c;

% handle the borders
if strcmpi(opt, 'zeropad')
    m = rr<1 | rr > Nr | cc<1 | cc > Nc;
    rr(m) = RR(m);
    cc(m) = CC(m);
elseif strcmpi(opt, 'circular')
    rr = mod(rr-1, Nr)+1;
    cc = mod(cc-1, Nc)+1;
else
    error('Option not supported');
end

% compute the final indices
ind = sub2ind([Nr Nc], rr, cc);