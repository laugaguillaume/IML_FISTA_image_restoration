

% Generate a random matrix

nRows=5000;
nCols=1000;
A=randn(nRows,nCols);

% Compute the its norm

l1infNorm=sum(max(abs(A')));
disp(['Input Matrix A: norm=' sprintf('%f', l1infNorm)]);

% Project the matrix into l1,inf balls with radius that are 
% factors of the original norms

F = [0.001, 0.005, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]; 

for i=1:length(F)    
    factor = F(i);
    tic;
    [B]=projL1Inf(A,factor*l1infNorm,ones(nRows,1));
    time = toc; 
    newnorm = sum(max(abs(B')));
    disp(['Projected Matrix B: factor=' sprintf('%f', factor) ' norm=' sprintf('%f', newnorm) ' time=' sprintf('%f', time)]);
end