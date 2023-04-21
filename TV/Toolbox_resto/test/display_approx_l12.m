close all
clear all
clc

addpath(genpath('../include'));
addpath(genpath('/Users/nellypustelnik/Dropbox/Travaux_en_cours/20XX_TOOLBOX/Optim/prox_repository/matlab/scalar'));
nu=0.01;
gamma=1;

%% display huber
% The v2 version is associated to the implementation made in 
% L. M. Brice?o-Arias and N. Pustelnik, Proximal or gradient steps for cocoercive operators, submitted, 2020. 
% This version is easier to compare to l1 norm else in the prox repository
% toolbox a normalization of gamma as gamma/nu is required to suit the
% behavior of the l1 norm


%
[X,Y] = meshgrid(-2:.1:2);
[nx,mx] = size(X);
fh = reshape(fun_hyperbolic_modif([X(:),Y(:)], gamma, nu, 2),nx,mx);
f = reshape(fun_L2_modif([X(:),Y(:)], gamma, 2),nx,mx);
fabs = reshape(fun_abs_modif([X(:)], gamma),nx,mx) +  reshape(fun_abs_modif([Y(:)], gamma),nx,mx);

figure;
subplot(331);
mesh(fh); title 'Hyperbolic';
subplot(332);
mesh(f); title 'L21';
subplot(333);
mesh(fabs); title 'L1';

%
g = grad_hyperbolic([X(:),Y(:)], gamma, nu, 2);
gradfh1 = reshape(g(:,1),nx,mx);
gradfh2 = reshape(g(:,2),nx,mx);

subplot(334);
mesh(gradfh1); title 'Grad: Hyperbolic x';
subplot(337);
mesh(gradfh2); title 'Grad:  Hyperbolic y';

%
p = prox_L2([X(:),Y(:)], gamma, 2);
proxf1 = reshape(p(:,1),nx,mx);
proxf2 = reshape(p(:,2),nx,mx);

subplot(335);
mesh(proxf1); title 'Prox: L12 x';
subplot(338);
mesh(proxf2); title 'Prox:  L12 y';

%
pl1_1 = reshape(prox_abs(X(:), gamma),nx,mx);
pl1_2 = reshape(prox_abs(Y(:), gamma),nx,mx);

subplot(336);
mesh(pl1_1); title 'Prox: L1 x';
subplot(339);
mesh(pl1_2); title 'Prox:  L1 y';
