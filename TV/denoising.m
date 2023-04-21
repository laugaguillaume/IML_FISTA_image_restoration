function [hk] = denoising(data,tau,f,prox,nlst)
%DENOISING Summary of this function goes here
%   Detailed explanation goes here
% compute the "guidance" image (single band)

%% Fixing Proximal Parameters
lambda= nlst.nltv_eta ;
normD  = f.beta;
gammak = tau/(normD)^2;
iter = nlst.iter;
eps = nlst.epsilon;
a=4; % inertia parameter of FISTA


%% Initializing variables
tmp = f.dir_op(data);
persistent yk
if isempty(yk)
    yk = zeros(size(tmp));
end
byk = yk;
hk = data - f.adj_op(yk);
fista = 1;

%% Criterion of convergence
crit=[0 0];
crit0 = 1e100;
it = 1;
gapc = eps+1;

%% Algorithm
while (abs(crit(it) - crit0)/abs(crit0) > eps) &&(it<iter)
    crit0 = crit(it);
    it = it + 1;
    %Save the dual variables
    yks = yk; byks = byk;
    %Update of dual variable
    tk    = f.dir_op(hk);
    tilyk = yks + gammak * tk;
    tempk = prox(tilyk/gammak, lambda/gammak);
    byk   = tilyk - gammak*tempk;
    %Update of inertia FISTA parameter
    fistas = fista;
    fista  = (it+a)/a;
    %Update
    yk = byk + (fistas-1)/fista*(byk-byks);
    %Update primal variable
    hk = hk - f.adj_op(yk-yks);
    tmp = f.dir_op(hk);
    crit(it) = 1/2*norm(hk(:)-data(:),'fro')^2 + lambda*penalisation_cost(tmp,nlst.reg);

end

end
