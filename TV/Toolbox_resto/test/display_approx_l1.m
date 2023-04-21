close all
clear all
clc

addpath(genpath('../include'));
addpath(genpath('/Users/nellypustelnik/Dropbox/Travaux_en_cours/20XX_TOOLBOX/Optim/prox_repository/matlab/scalar'));
nu=1e-1;
gamma=1;

%% display huber
% The v2 version is associated to the implementation made in 
% L. M. Brice?o-Arias and N. Pustelnik, Proximal or gradient steps for cocoercive operators, submitted, 2020. 
% This version is easier to compare to l1 norm else in the prox repository
% toolbox a normalization of gamma as gamma/nu is required to suit the
% behavior of the l1 norm



x = [-2:0.001:2];
figure;
subplot(131);
plot(x, fun_huber_v2_modif(x,gamma,nu),'Linewidth',2);
hold on;
plot(x, fun_huber_modif(x,gamma,nu),'r--','Linewidth',2);
plot(x, fun_abs_modif(x,gamma),'k--','Linewidth',2);
title 'f(x)';
grid on;
%
subplot(132);
plot(x, grad_huber_v2(x,gamma,nu),'Linewidth',2);
hold on;
plot(x, grad_huber(x,gamma,nu),'r','Linewidth',2);
title 'Gradf(x)';
grid on;
%
subplot(133);
plot(x, prox_huber_v2(x,gamma,nu),'Linewidth',2);
hold on;
plot(x, prox_huber(x,gamma,nu),'r--','Linewidth',2);
plot(x, prox_abs(x,gamma),'k--','Linewidth',2);
title 'Proxf(x)';
grid on;
legend('Huber article Luis','Huber Prox repository','L1');

