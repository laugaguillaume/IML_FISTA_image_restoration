clear all
close all
%%
addpath(genpath('TV'))
addpath(genpath('Functions'))
addpath InfoTransfer/
addpath Proximity_operators/
%addpath('/Applications/MATLAB_R2021b.app/toolbox/Wavelab850')
%WavePath;
%%
image_choice = 'JWST 256 gray';  % Image choice
param.std    = 1e-2;           % Variance du bruit
choice_blur  = 'Inpainting 50';  
snr = [0]
%%
[X,param,Ar,Ac,A1,A2,A3]   = create_data(image_choice,choice_blur,param);
[param.N,param.M,channel] = size(X);
TAU = 1; %STEP 
% TOTAL VARIATION %
param_reg.norm_type = {'L2reg'};  % options: 'nuclearreg'
param_reg.reg = {'L1'}; %L1, nuclear
param_reg.TVtype= 'TV';
%%
param.lambda = 5e-3;
p = 2;
%%
% Algo
max_iteration_number = 3;    % Nb max d'itération au niveau fin  
param.p = p;
param.starting_iterate     = 'Wiener filter';              % Initialisation algo     
bk_coarse            = 0;                     % Backtracking niveau grossier     
bk_fine              = 'No';       % Avec ou sans backtracking
param.epsilon        = 1e-8;
param.nb_level       = 1;
param.d		     = 1;
param.iterations_number   = [max_iteration_number]; % Number of iterations per level in descending order 

% Construction of struct containing each level iterations number
for level=param.nb_level:-1:1
    level_number = matlab.lang.makeValidName(['maxit_level' num2str(level)]);
    param.(level_number) = param.iterations_number(param.nb_level-level+1);
end
mgm =param.nb_level;
% Penalization term
param_reg.nltv_eta  = param.lambda;
param_reg.neigh_rad = 2;
param_reg.blk_rad   = 2;
param_reg.delta     = 35;
param_reg.iter      = 2e2;
param_reg.epsilon   = param.epsilon;
param.reg_type = [param_reg.reg,'_',param_reg.TVtype];
if strcmpi(param_reg.TVtype, 'NLTV')
    x_ref = mean(param.Z,3);
    w = get_fov_weights(x_ref, param_reg.delta, param_reg.blk_rad, param_reg.neigh_rad);
    D_op = get_nltv(w);
elseif strcmpi(param_reg.TVtype, 'TV')
    D_op = get_tv;
end
prox = get_prox(param_reg.norm_type);
param.dir_op       = @(x) D_op.dir_op(x);
param.adj_op      = @(x) D_op.adj_op(x);
param.functions.gh = @(U) param.lambda.*penalisation_cost(D_op.dir_op(U),param_reg.reg);%tlv(U,param.tlv_choice);
param.proximal.proxgh = @(U,GAMMA,param_reg) denoising(U,param.lambda*GAMMA,D_op,prox,param_reg);
param.proximal.proxg = @(U,GAMMA) prox(U,param.lambda.*GAMMA);
% Config blur for optimization
switch choice_blur    
    case 'Gauss small'
        param.functions.fh = @(U) blur_L2_color_function(U,param.Z,Ac,Ar,Acolor); 
        param.gradient.gradfh = @(U) blur_L2_color_gradient(U,param.Z,Ac,Ar,Acolor); 
    case 'Gauss big'
        param.functions.fh = @(U) blur_L2_color_function(U,param.Z,Ac,Ar,Acolor); 
        param.gradient.gradfh = @(U) blur_L2_color_gradient(U,param.Z,Ac,Ar,Acolor); 
    case 'Inpainting 50'
        param.functions.fh = @(U) L2_inpainting_function(U,param.Z,A1,A2,A3); 
        param.gradient.gradfh = @(U) L2_inpainting_gradient(U,param.Z,A1,A2,A3);
    case 'Inpainting 90'
        param.functions.fh = @(U) L2_inpainting_function(U,param.Z,A1,A2,A3); 
        param.gradient.gradfh = @(U) L2_inpainting_gradient(U,param.Z,A1,A2,A3);  
end

levelAlpha = matlab.lang.makeValidName(['level' num2str(param.nb_level) 'alpha']);
param.(levelAlpha) = TAU;

% Starting Iterate choice
switch param.starting_iterate
    case 'Zero'
        X0 = zeros(param.N,param.M);
    case 'Z'
        X0 = param.Z;
    case 'Wiener filter'
        if channel == 1
            X0 = wiener2(param.Z,[5 5]);
        elseif channel == 3
            X0 = zeros(size(param.Z));
            Z1 = param.Z(:,:,1); Z2 = param.Z(:,:,2); Z3 = param.Z(:,:,3);
            X0(:,:,1) = wiener2(Z1,[5 5]);
            X0(:,:,2) = wiener2(Z2,[5 5]);
            X0(:,:,3) = wiener2(Z3,[5 5]);
        end
end

%% Convergence simulations
[xk_1f,x_iterations_1f,Fun_1f,counter_descent_1f,coarse_descent_f_1f,time_steps_1f]= ...
     IMLFista_Color(param.nb_level,param.functions.fh,param.functions.gh,param.gradient.gradfh,param.proximal.proxgh,X0,X,param_reg,param,param.nb_level,bk_fine,bk_coarse,'FISTA');
f_1f = Fun_1f(1,:); g_1f = Fun_1f(2,:);
finale_1f = f_1f(end) + g_1f(end);
temps_1f = time_steps_1f(end);
filename = [image_choice,'_',choice_blur,'_',num2str(param.std),'_',param_reg.reg{1},param_reg.TVtype,'_',num2str(param.lambda),'_',num2str(TAU)];
Z = param.Z;
k = 0
while isfile([filename,num2str(k),'.mat'])
     k = k+1;
end
snr1 = snrdeblur(xk_1f,X);
snr = [snr, snr1];
filename = [filename,num2str(k),'.mat'];
save(filename,'Z','Ar','Ac','A1','A2','A3','D_op','X0','xk_1f','x_iterations_1f','Fun_1f','time_steps_1f','finale_1f','temps_1f','snr1');

%% IML FISTA

load(filename) %load Ac Ar D_op
%%
% Algo   % Nb max d'itération au niveau fin
max_it_coarse        = 5;            % Nb max d'itération au niveau grossier
param.moreau          = 1.1;           % Paramètre de moreau
param.p               = 2;             % Nb de fois où l'on va à l'itération grossière
starting_iterate     = 'Wiener filter';              % Initialisation algo     
bk_coarse            = 0; bk_fine              = 'No';           % Backtracking niveau grossier     
param.F_star         = 1e-10;      % Critère d'arrêt
param.epsilon        = 1e-8;
param.nb_level       = 5;
param.d		     = 1;
param.iterations_number   = [max_iteration_number,5,5,5,5,5,5,5]; % Number of iterations per level in descending order 

% Construction of struct containing each level iterations number
for level=param.nb_level:-1:1
    level_number = matlab.lang.makeValidName(['maxit_level' num2str(level)]);
    param.(level_number) = param.iterations_number(param.nb_level-level+1);
end
% Penalization term
TAU = 1;
% TOTAL VARIATION %
param_reg.nltv_eta  = param.lambda;
param_reg.neigh_rad = 2;
param_reg.blk_rad   = 2;
param_reg.delta     = 35;
param_reg.iter      = 2e2;
param_reg.epsilon   = param.epsilon;
param.reg_type = [param_reg.reg,'_',param_reg.TVtype];
prox = get_prox(param_reg.norm_type);
param.dir_op       = @(x) D_op.dir_op(x);
param.adj_op      = @(x) D_op.adj_op(x);
param.functions.gh = @(U) param.lambda.*penalisation_cost(D_op.dir_op(U),param_reg.reg);%tlv(U,param.tlv_choice);
param.proximal.proxgh = @(U,GAMMA,param_reg) denoising(U,param.lambda*GAMMA,D_op,prox,param_reg);
param.proximal.proxg = @(U,GAMMA) prox(U,param.lambda.*GAMMA);
% Config blur for optimization
switch choice_blur
    case 'Gauss small'
        param.functions.fh = @(U) blur_L2_color_function(U,Z,Ac,Ar,Acolor); 
        param.gradient.gradfh = @(U) blur_L2_color_gradient(U,Z,Ac,Ar,Acolor); 
    case 'Gauss big'
        param.functions.fh = @(U) blur_L2_color_function(U,Z,Ac,Ar,Acolor); 
        param.gradient.gradfh = @(U) blur_L2_color_gradient(U,Z,Ac,Ar,Acolor); 
    case 'Inpainting 50'
        param.functions.fh = @(U) L2_inpainting_function(U,Z,A1,A2,A3); 
        param.gradient.gradfh = @(U) L2_inpainting_gradient(U,Z,A1,A2,A3);
    case 'Inpainting 90'
        param.functions.fh = @(U) L2_inpainting_function(U,Z,A1,A2,A3); 
        param.gradient.gradfh = @(U) L2_inpainting_gradient(U,Z,A1,A2,A3); 
end

% Information transfer choice
param.wavelet_choice = 'Symmlet';
param.v_m = 10; param.TypeBorder = 'periodic' ;
param.qmfphi  = MakeONFilter(param.wavelet_choice,param.v_m);
levelAlpha = matlab.lang.makeValidName(['level' num2str(param.nb_level) 'alpha']);
param.(levelAlpha) = TAU;
Ac_H = Ac; Ar_H = Ar; ZH = Z; A1_H = A1; A2_H = A2; A3_H = A3;
param.alpha_moreau = 1/(1+param.moreau*param.lambda);
lambda_reg = param.lambda;
N = param.N;
M  = param.M;
for level = param.nb_level-1:-1:1
    N = N/2; M = M/2;
    lambda_reg = lambda_reg/4;
    levelR = matlab.lang.makeValidName(['level' num2str(level) 'R']);
    levelP = matlab.lang.makeValidName(['level' num2str(level) 'P']);
    levelAlpha = matlab.lang.makeValidName(['level' num2str(level) 'alpha']);
    levelfH = matlab.lang.makeValidName(['level' num2str(level) 'fH']);
    levelgH = matlab.lang.makeValidName(['level' num2str(level) 'gH']);
    levelgradfH = matlab.lang.makeValidName(['level' num2str(level) 'gradfH']);
    levelgradgH = matlab.lang.makeValidName(['level' num2str(level) 'gradgH']);
    levelproxgH = matlab.lang.makeValidName(['level' num2str(level) 'proxgH']);
    levelopdirectH = matlab.lang.makeValidName(['level' num2str(level) 'dir_opH']);
    levelopadjointH = matlab.lang.makeValidName(['level' num2str(level) 'adj_opH']);
    % Information transfer definition at coarser levels
    [Phi,Phisym_dir,Phisym_inv] = WaveletMatrix(2*N,param.qmfphi,param.TypeBorder) ;
    R =Phisym_dir; P = Phisym_inv;
    ZH = D3RP(ZH,R);%reshape(R*ZH(:),N,M);
    param.blur.(levelR) = R; 
    param.blur.(levelP) = P;
    switch choice_blur
        case 'Gauss small'
            Ac_H = R*Ac_H*R'; Ar_H = R*Ar_H*R';
            Ar_H(Ar_H<1e-10) = 0; Ac_H(Ac_H<1e-10) = 0;
            Ar_H = sparse(Ar_H); Ac_H = sparse(Ac_H);
            param.functions.(levelfH) = @(U) blur_L2_color_function(U,ZH,Ac_H,Ar_H,Acolor);
            param.gradient.(levelgradfH) = @(U) blur_L2_color_gradient(U,ZH,Ac_H,Ar_H,Acolor);
        case 'Gauss big'
            Ac_H = R*Ac_H*R'; Ar_H = R*Ar_H*R';
            Ar_H(Ar_H<1e-10) = 0; Ac_H(Ac_H<1e-10) = 0;
            Ar_H = sparse(Ar_H); Ac_H = sparse(Ac_H);
            param.functions.(levelfH) = @(U) blur_L2_color_function(U,ZH,Ac_H,Ar_H,Acolor);
            param.gradient.(levelgradfH) = @(U) blur_L2_color_gradient(U,ZH,Ac_H,Ar_H,Acolor);
        case 'Inpainting 50'
            if channel == 3
                A1_H = A1_H(1:2:end,1:2:end); A2_H = A2_H(1:2:end,1:2:end); A3_H = A3_H(1:2:end,1:2:end);
            elseif channel == 1
                A1_H = A1_H(1:2:end,1:2:end); A2_H = []; A3_H = [];
            end
            param.functions.(levelfH) = @(U) L2_inpainting_function(U,ZH,A1_H,A2_H,A3_H); 
            param.gradient.(levelgradfH) = @(U) L2_inpainting_gradient(U,ZH,A1_H,A2_H,A3_H);
        case 'Inpainting 90'
            if channel == 3
                A1_H = A1_H(1:2:end,1:2:end); A2_H = A2_H(1:2:end,1:2:end); A3_H = A3_H(1:2:end,1:2:end);
            elseif channel == 1
                A1_H = A1_H(1:2:end,1:2:end); A2_H = []; A3_H = [];
            end
            param.functions.(levelfH) = @(U) L2_inpainting_function(U,ZH,A1_H,A2_H,A3_H); 
            param.gradient.(levelgradfH) = @(U) L2_inpainting_gradient(U,ZH,A1_H,A2_H,A3_H); 
    end
    if strcmpi(param_reg.TVtype, 'NLTV')
        x_ref = mean(denoising_TV(ZH,param.lambda),3);
        w = get_fov_weights(x_ref, param_reg.delta, param_reg.blk_rad, param_reg.neigh_rad);
        D_opH = get_nltv(w);
    elseif strcmpi(param_reg.TVtype, 'TV')
        D_opH = get_tv;
    end
    param.op.(levelopdirectH)       = @(U) D_opH.dir_op(U);
    param.op.(levelopadjointH)     = @(U) D_opH.adj_op(U);
    
    param.functions.(levelgH) = @(U) lambda_reg.*penalisation_cost(U,param_reg.reg);
    param.proximal.(levelproxgH) = @(U,GAMMA) prox(U,lambda_reg.*GAMMA); % U is here in H,V format
end

%% Convergence simulations
[xk,x_iterations,Fun,counter_descent,coarse_descent_f,time_steps,snr_steps]= ...
    IMLFista_Color(param.nb_level,param.functions.fh,param.functions.gh,param.gradient.gradfh,param.proximal.proxgh,X0,X,param_reg,param,param.nb_level,bk_fine,bk_coarse);
%%
filename_IML = ['IMLFISTA_',filename];
k = 0;
while isfile([filename_IML,num2str(k),'.mat'])
     k = k+1;
end
filename_IML = [filename_IML,num2str(k),'.mat'];
save(filename_IML,'xk','x_iterations','Fun','time_steps','snr_steps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(filename,'Fun_1f','time_steps_1f','Z','finale_1f','xk_1f','x_iterations_1f')
load(filename_IML,'Fun','time_steps','xk','x_iterations','snr_steps')
%% NORMALIZE FUNCTION VALUES
N = max_iteration_number
F = sum(Fun,1);
F_normalized = (F(1:N+1))./F(1);
F_1f = sum(Fun_1f,1); F_1f = F_1f(1:N+1);
F_normalized = max(F_normalized,0);
F_1f_normalized = ( F_1f)./F(1);
F_1f_normalized = max(F_1f_normalized,0);
%% PLOT OBJECTIVE FUNCTION EVOLUTION
figure(1)
f = gcf;
f.Units = 'normalized'
f.InnerPosition = [0 0.0756 0.4528 0.4122];
semilogy(time_steps(1:N+1),F_1f_normalized,'r-+','LineWidth',2.5,'MarkerSize',10)
hold on 
semilogy(time_steps(1:N+1),F_normalized,'b-+','LineWidth',2.5,'MarkerSize',10)
xlabel('CPU time (s)')
ylabel('F_h')
xlim([0,max(time_steps(N+1))])
legend('FISTA','IML FISTA','FontSize',25)
title('F_h evolution with CPU time')
set(gca,'FontSize',20)
%% PLOT IMAGES
figure(2)
t = tiledlayout(2,3,'TileSpacing','compact');
ax1 = nexttile
imagesc(X); title('X'); colormap gray
set(ax1,'xtick',[],'ytick',[],'ylabel',[])
set(ax1,'FontSize',20)
ax3 = nexttile;
imagesc(x_iterations_1f); title(['FISTA ',num2str(p),'it']); colormap gray
set(ax3,'xtick',[],'ytick',[],'ylabel',[])
set(ax3,'FontSize',20)
ax5 = nexttile;
imagesc(xk_1f); title('FISTA'); colormap gray
set(ax5,'xtick',[],'ytick',[],'ylabel',[])
set(ax5,'FontSize',20)
ax2 = nexttile;
imagesc(Z); xlabel('Z','fontweight','bold'); colormap gray
set(ax2,'xtick',[],'ytick',[],'title',[],'ylabel',[])
set(ax2,'FontSize',20)
ax4 = nexttile;
imagesc(x_iterations); xlabel(['IML FISTA ',num2str(p),'it'],'fontweight','bold'); colormap gray
set(ax4,'xtick',[],'ytick',[],'title',[],'ylabel',[])
set(ax4,'FontSize',20)
ax6 = nexttile;
imagesc(xk); xlabel('IML FISTA','fontweight','bold'); colormap gray
set(ax6,'xtick',[],'ytick',[],'title',[],'ylabel',[])
set(ax6,'FontSize',20)
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'xy')
