 function mask = get_random_indices(sz, decimation_rate)
%function mask = get_random_indices(sz, decimation_rate)
%
%  Created on: --/--/--, Nelly Pustelnik
%
% The function computes a logical mask where the 'zero' points are in 
% random positions and their number is about the rate specified as input.


N = prod(sz);
n_indRand = bino_gen(N, decimation_rate); %random('bino',N, decimation_rate);

tmp = rand(1,N);
tmps = sort(tmp);

mask = tmp < tmps(n_indRand);
mask = reshape(mask, sz);



function x = bino_gen(n, p)

if (0 <= p && p <= 1) && (0 <= n && round(n) == n)    
    
    x = sum( rand([1 1 n]) < p, 3 );
    
else
    error('The inputs are not correct.');
end