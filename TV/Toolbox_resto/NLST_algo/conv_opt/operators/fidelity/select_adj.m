function y = select_adj(x, mask)

y = zeros( size(mask) );
y(mask) = x;