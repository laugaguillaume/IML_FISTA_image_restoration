 function [x, y, it, time, varargout] = CP_uncertainties_LJ(x_init, f, g, h, opt, varargin)
%function [x it time varargout] = LFBF_epi(x_init, f, g, h, opt, varargin)
%
%  Created on: 19/01/12 - Giovanni Chierchia
%
% Modified by N. Pustelnik
% July 8th 2021.

% default parameters
if nargin < 4
    error('Too few parameters');
end
if nargin < 5 || isempty(opt)
    opt = [];
end
if ~isfield(opt, 'gamma') || isempty(opt.gamma)
    opt.gamma = 1;
end
if ~isfield(opt, 'tol') || isempty(opt.tol)
    opt.tol = 10^-4;
end
if ~isfield(opt, 'iter') || isempty(opt.iter)
    opt.iter = 500;
end
if ~isfield(opt, 'show')
    opt.show = [];
end

% initialize the time vector
time  = zeros(1, opt.iter);

% initialize the criterion vectors
n_out = nargout - 4;
n_in  = nargin  - 5;
n_fun = min(n_in, n_out);
varargout = cell(1,n_fun);
for i = 1:n_fun
    varargout{i} = zeros(1, opt.iter);
end

% open a figure
if ~isempty(opt.show)
    fig = figure;
end

% open a wait bar
hdl = wait_bar(0, 'Running CP (epi)...');


% check 'g' input
if ~iscell(g)
    if isfield(g, 'project_epi')
        g = {g []};
    else
        error('project_epi not found in ''g''');
    end
elseif ~isfield(g{1}, 'project_epi')
    error('project_epi not found in ''g''');
end
    


%---------------------------------%
%----- ALGORITHM STARTS HERE -----%
%---------------------------------%
% default functions
if isempty(f)
    f.prox = @(y,gamma) y;
end
if length(g)==1
    g{2}.prox   = @(y,gamma) y;
    g{2}.dir_op = @(y) y;
    g{2}.adj_op = @(y) y;
    g{2}.beta = 0;    
end
if isempty(h) 
    h.grad = @(y) 0;
    h.beta = 0;
end


% select the step-size ---------------------------------------------------%
if ~isfield(opt, 'gamma') || isempty(opt.gamma)
    opt.gamma = 1;
end
gamma = opt.gamma;
beta = g{1}.beta +  g{2}.beta +1
tau = 0.99/(h.beta/2 + gamma*beta);
fprintf('gamma=%3.2f et tau =%3.2f\n',gamma,tau)

%-------------------------------------------------------------------------%


% initialize the solution
x = x_init;
y = x;
% initialize the temporary variables
w1 = zeros(size(g{1}.dir_op(x)));
w2 = zeros(size(g{2}.dir_op(x)));
u = zeros(size(w1));
v = zeros(size(w2));
[~,m,l,h] = size(u);
zeta1 = zeros(m,l,h);
zeta2 = zeros(size(v));
d_zeta1 = zeta1;
d_zeta2 = zeta2;
d_y = y;
% perform the M+LFBF algorithm
for it = 1:opt.iter
    
    tic;
   
    x_hat = x + gamma * (g{1}.adj_op(w1) + g{2}.adj_op(w2)) ;
    u_hat = u - gamma * w1;    
    v_hat = v - gamma * w2;  
    zeta1_hat = zeta1 - gamma * d_zeta1;     
    zeta2_hat = zeta2 - gamma * d_zeta2;        
    y_hat = y - gamma * d_y;   
    
    %p_x = project_S(x_hat, opt);
    val_inv = inv([(1+gamma), -1;-1 (1+gamma)]);
    p_x = val_inv(1,1)*x_hat + val_inv(1,2)*y_hat;
    p_y = val_inv(2,1)*x_hat + val_inv(2,2)*y_hat ; 
    
    [p_u, p_epi1] = g{1}.project_epi(u_hat, zeta1_hat);
    [p_v, p_epi2] = f(2).project(v_hat, zeta2_hat); 


   
   
    % 2/5: function 'f'

    
    w1 = w1  + tau * (2*p_u - u - g{1}.dir_op(2*p_x-x));
    w2 = w2  + tau * (2*p_v - v - g{2}.dir_op(2*p_x-x));
    d_zeta1 = d_zeta1 + tau*(2*p_epi1 - zeta1);
    d_zeta2 = d_zeta2 + tau*(2*p_epi2 - zeta2);    
    d_y = d_y + tau*(2*p_y - y);    
    
    p_w1 = zeros(size(w1));
    p_w2 = zeros(size(w2));  
    

    [p_eta1, p_eta2] = f(1).project(d_zeta1/tau,d_zeta2/tau);

    
    w1 = w1 - tau*p_w1;
    w2 = w2 - tau*p_w2;  
    

    
    d_zeta1 = d_zeta1 - tau*p_eta1;
    d_zeta2 = d_zeta2 - tau*p_eta2;    
    
    d_y     = d_y - tau * project_S(d_y/tau, opt);    
    

    % 4/4: primal forward step
    x_old = x;
    x = p_x;
    y = p_y;
    u = p_u;
    v = p_v;
    zeta1 = p_epi1;
    zeta2 = p_epi2;
    

    time(it) = toc;

Fx = shiftdim(sqrt(sum(g{1}.dir_op(x).^2, 1)),1);
eval_crit = opt.noisecoef* sum(abs( g{2}.dir_op(x) - opt.noisy(opt.valid) ))+ opt.lambda*sum( Fx(:) ) ;
fprintf('norm = %3.2f\t crit_real = %3.2f\n',norm(x-y,'fro'),eval_crit);
   
    %Fx = shiftdim(sqrt(sum(u.^2, 1)),1);
    %dEpi1 =sum(  max(Fx(:) - zeta1(:),0));
    %fprintf('dEpi1 = %3.2f\n',dEpi1);     
    
    % dEpi2 = sum(sum(max(abs(v - opt.noisy(:)) - zeta2,0)));
    %  fprintf('dEpi2 = %3.2f\n',dEpi2);
    
% 

%      Fx = shiftdim(sqrt(sum(u_hat.^2, 1)),1);
%      dEpi2 = max(sum( Fx(:) ) - zeta2_hat,0);
%     
%     [tmp1,tmp2] =   f(1).project(zeta1,zeta2);
%     dHalfSpace = norm(zeta1(:)-tmp1(:),'fro')^2 + norm(zeta2(:)-tmp2(:),'fro')^2;
%     fprintf('dHS = %3.2f\n',dHalfSpace);
%     tmp1 = g{1}.dir_op(x);
%     tmp2 = g{2}.dir_op(x) ;
%     dL1 = norm(tmp1(:) - u(:));
%     dL2 = norm(tmp2(:) - v(:)) ;   
%     fprintf('dL1 = %3.2f - dL2 = %3.2f\n',dL1,dL2);
     
%    fprintf('Epi1 = %3.2f \t Epi2 = %3.2f \t dHalfSpace = %3.2f \t dL1 = %3.2f \t dL2 = %3.2f\n\n',dEpi1,dEpi2,dHalfSpace,dL1,dL2);
    %--- OPERATIONS NOT RELATED TO THE ALGORITHM ---%
    
    
    % check the stop criterion
    if norm(x(:) - x_old(:)) < opt.tol * norm(x_old(:))  % || all(x(:) == x_old(:))
        break;
    end
    
    
    % compute the criterion
    varargout{1}(it) = varargin{1}(x);
    varargout{2}(it) = varargin{2}(y);
    varargout{n_fun}(it) = varargin{n_fun}(x,y);
    %fprintf('dCalpha = %3.2f \t dS = %3.2f\n',varargout{1}(it),varargout{2}(it));
    
    % show the current solution
    if ~isempty(opt.show) && mod(it,10) == 0
        set(0,'CurrentFigure',fig);
        opt.show(x);
        title( sprintf('Iteration %d', it) );
    end
    
    
    % update the wait bar
    wait_bar( it/opt.iter, hdl );
end



%-------------------------------%
%----- ALGORITHM ENDS HERE -----%
%-------------------------------%



% close the wait bar
close(hdl);

% close the figure
if ~isempty(opt.show)
   close(fig);
end

% compute the times
time = cumsum(time(1:it));

% compute the criterion
    varargout{1}(it) = varargin{1}(x);
    varargout{2}(it) = varargin{2}(y);
    varargout{n_fun}(it) = varargin{n_fun}(x,y);

% cut the zeros
for i = 1:n_fun
    varargout{i} = varargout{i}(1:it);
end