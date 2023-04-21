function [Phi,Phi_dir,Phi_inv] = WaveletMatrix(N,qmf,mode)

% Build the large-sparse matrix Phi whose multiplication with an image X
% corresponds to the action of the scaling wavelet function $\Phi$ on X
% Inputs:
%   N : size of the 1D signal (or the N-by-N image)
%   qmf : quadrature mirror filter coeffs obtained with MakeONFilter.m
%   mode : border effect mode
%       - 'none' : no border effect correction
%       - 'sym' : symetric border effect, as if  X was symetrized 
%                 as follows [... w x y z | y x w]
%       - 'sym-twosided' : as 'sym' but symetrization is applied to the left
%       and right parts of X [d c b | a b c d ... w x y z | y x w]
%       and similarly for the top and the bottom of an image X
%       - 'periodic' : poriodized border effect, as if X was
%       periodic as follows : [a b c ... x y z | a b c]
%       -'periodic-twosided' : as 'periodic' but applied to the left and to 
%       the right parts of X [x y z | a b c ... x y z | a b c]
%       and similarly for the top and the bottom of an image X
%
% Outputs:
%   Phi : N-by-N matrix. Phi*X is equivalent to projecting X into
%       the lower resolution subspace without decimation by a factor 2
%   Phi_dir : (N/2)-by-N matrix : row-wise decimation of Phi by a factor 2.
%       Phi_dir*X is equivalent to computing the lower resolution
%       approximation of  X 
%   Phi_inv : N-by-(N/2) matrix : Phi_inv*X is equivalent to computing the
%   expansion of X at a finer resolution. If mode is 'none' or 'periodic',
%   then Phi_inv is simply equal to Phi_dir'.

% PG : feb. 2022

L = length(qmf) ;

switch mode
    case 'none'
        Phi = toeplitz([qmf(1) zeros(1,N-1)],[qmf zeros(1,N-L)]) ;
        Phi_dir = Phi(1:2:N,:) ;
        Phi_inv = Phi_dir' ;
    case 'sym'
        Phi = toeplitz([qmf(1) zeros(1,N-1)],[qmf zeros(1,N-L)]) ;
        SymR = flipud(toeplitz([qmf(L) zeros(1,L-2)],qmf(L:-1:2))) ;
        Phi(N-L+2:N,N-L+1:N-1) = Phi(N-L+2:N,N-L+1:N-1) + SymR ;
        Phi_dir = Phi(1:2:N,:) ;
        Phi_inv = Phi(:,1:2:N) ;
    case 'mirror-twosided'
        Phi = toeplitz([qmf(L/2+1:-1:1) zeros(1,N-L/2-1)],[qmf(L/2+1:L) zeros(1,N-L/2)]) ;
        SymL = fliplr(toeplitz([qmf(1) zeros(1,L/2-1)],qmf(1:L/2))) ;
        SymR = flipud(toeplitz([qmf(L) zeros(1,L/2-2)],qmf(L:-1:L/2+2))) ;
        Phi(1:L/2,2:L/2+1) = Phi(1:L/2,2:L/2+1) + SymL ;
        Phi(N-(L/2-1)+1:N,N-(L/2-1):N-1) = Phi(N-(L/2-1)+1:N,N-(L/2-1):N-1) + SymR ;
        Phi_dir = Phi(1:2:N,:) ;
        Phi_inv = Phi(:,1:2:N) ;
    case 'periodic'
        Phi = toeplitz([qmf(1) zeros(1,N-1)],[qmf zeros(1,N-L)]) ;
        PeriodR = (toeplitz([qmf(L) zeros(1,L-2)],qmf(L:-1:2)))' ;
        Phi(N-L+2:N,1:L-1) = Phi(N-L+2:N,1:L-1) + PeriodR ;
        Phi_dir = Phi(1:2:N,:) ;
        Phi_inv = Phi_dir' ;
    case 'periodic-twosided'
        Phi = toeplitz([qmf(L/2+1:-1:1) zeros(1,N-L/2-1)],[qmf(L/2+1:L) zeros(1,N-L/2)]) ;
        PeriodL = toeplitz([qmf(1) zeros(1,L/2-1)],qmf(1:L/2)) ;
        PeriodR = (toeplitz([qmf(L) zeros(1,L/2-2)],qmf(L:-1:L/2+2)))' ;
        Phi(1:L/2,N-L/2+1:N) = Phi(1:L/2,N-L/2+1:N) + PeriodL ;
        Phi(N-L/2+2:N,1:L/2-1) = Phi(N-L/2+2:N,1:L/2-1) + PeriodR ;
        Phi_dir = Phi(1:2:N,:) ;
        Phi_inv = Phi_dir' ;
end

Phi = sparse(Phi) ;
Phi_dir = sparse(Phi_dir) ;
Phi_inv = sparse(Phi_inv) ;


