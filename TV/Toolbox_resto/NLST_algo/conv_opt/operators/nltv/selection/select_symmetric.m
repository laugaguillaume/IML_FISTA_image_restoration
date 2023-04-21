 function y = select_symmetric(x)
%function y = select_symmetric(x)
%
%  Created on: 4/10/2012 - Giovanni Chierchia
%
%
% This functions processes a 3D matrix of non-negative weights according to
% the "semi-local" algorithm proposed by Gilboa-Osher in "NonLocal Linear 
% Image Regularization and Supervised Segmentation".
% 
% In order to work properly, the input matrix must contain the weights 
% computed in each position "n" of a neighbourhood "N_\ell" composed by 
% the (2*r+1)x(2*r+1) positions located around the position "\ell". These
% weights don't have to be normalized.


% number of additional weights to be selected
m = 5;

% initialize
y = zeros( size(x) );

% get the radius
Nb = size(x,3);
rad = (sqrt(Nb) - 1) / 2;

% handle 3x3 neighbourhoods
if rad == 1
    y = x;
    return;
end

% select the central pixel and the 4 nearest neighbors
sz = [2*rad+1 2*rad+1];
bb(1) = sub2ind(sz,rad+1,rad);     % left neighbour
bb(2) = sub2ind(sz,rad+1,rad+2);   % right neighbour
bb(3) = sub2ind(sz,rad,  rad+1);   % up neighbour
bb(4) = sub2ind(sz,rad+2,rad+1);   % down neighbour
bb(5) = sub2ind(sz,rad+1,rad+1);   % central neighbour
y(:,:,bb) = x(:,:,bb);

% select the "m" most similar neighbors
[ng i] = sort(x-y, 3, 'descend');       % x-y is a shorthand for setting x(:,:,bb) = 0
for j=1:m
    mask = idx2mask( i(:,:,j), Nb );
    y(mask) = x(mask);
end

% ensure simmetric weights
y = ensure_symmetry(y, x, rad, i, m);





function y = idx2mask(x, Nb)

y = false([size(x) Nb]);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        y( i, j, x(i,j) ) = true;
    end
end



function y = ensure_symmetry(y, x, rad, idx, m)

% get the size
[Nr Nc Nb] = size(x);
sz = [2*rad+1 2*rad+1];

% extend the borders
xx = padarray(x, [rad rad]);
yy = padarray(y, [rad rad]);

% build the pixel grid 
[c r] = meshgrid(1:Nc, 1:Nr);

% count of added positions
count = zeros(Nr,Nc);

% scroll the non-selected positions
for j = m+1:size(idx,3)
    
    % COMMENTS:
    % for each "ell", 
    %  1. we take the weight w(ell,b), where "b" is the position of the highest non-selected weight in N_ell.
    %  2. if w(b,ell) == 0, we set w(b,ell) = w(ell,b) and we increase the counting
    %  3. we move to w(ell,b+1), the next highest non-selected weight in N_ell and we repeat from 1

    % positions of "b" in the image domain
    b = idx(:,:,j);
    [rr cc] = ind2sub(sz, b);
    rr = r + rr - (rad+1);
    cc = c + cc - (rad+1);
    
    % position of "ell" in the neighborhood N_b
    pr = (rad+1) + r-rr;
    pc = (rad+1) + c-cc;
    bb = sub2ind(sz, pr, pc);

    % position of weight w(ell,b)
    i1 = sub2ind([Nr Nc Nb]+2*[rad rad 0], r+rad, c+rad, b);
    
    % position of weight w(b,ell)
    i2 = sub2ind([Nr Nc Nb]+2*[rad rad 0], rr+rad, cc+rad, bb);
    
    % mask of w(b,ell) ~= 0
    mask = yy(i2)~=0 & count < m;
    
    % insert w(ell,b) in N_ell
    yy(i1) = ~mask .* yy(i1) + mask .* xx(i1);
    
    % increment
    count = count + mask;
    
    % exit when done
    if all( count(:) == m )
        break;
    end
end

y = yy(rad+(1:Nr), rad+(1:Nc), :);







%-------------%
% OLD CODE
%-------------%
% [c r] = meshgrid(1:Nc, 1:Nr);
% for j = 1:m
%     
%     % COMMENTS:
%     %  - We want to assign the value w(r,c,b) to its symmetric weight.
%     %  - The weights outside the image boundary have no symmetric.
% 
%     % positions of w(:,:,b) in the image domain
%     b = i(:,:,j);
%     [rr cc] = ind2sub(sz, b);
%     rr = r + rr - (rad+1);
%     cc = c + cc - (rad+1);
%     
%     % check if there is no symmetric weight
%     if any(rr(:)<1 | rr(:)>Nr) || any(cc(:)<1 | cc(:)>Nc)
%         error('Not expected');
%     end
%     
%     % position of symmetric in the neighborhood of w(r,c,b)
%     pr = (rad+1) + r-rr;
%     pc = (rad+1) + c-cc;
%     bb = sub2ind(sz, pr, pc);
% 
%     % add a symmetric weight
%     in  = sub2ind([Nr Nc Nb],  r,  c, b);
%     out = sub2ind([Nr Nc Nb], rr, cc, bb);
%     y(out) = y(in);
%     
% end