 function [x, it, time, varargout] = CP_uncertainties_POCS(x_init, f, g, h, opt, varargin)
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
n_out = nargout - 3;
n_in  = nargin  - 2;
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
hdl = wait_bar(0, 'Running M-LFBF (epi)...');


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
beta = g{1}.beta +  g{2}.beta;
tau = 0.99/(h.beta/2 + gamma*beta);
%fprintf('gamma=%3.2f et tau =%3.2f\n',gamma,tau)

%-------------------------------------------------------------------------%

%
% range constraint
%
fprojS = @(y,gamma) (y + gamma*x_init)/(1+gamma);
% initialize the solution
x = x_init;

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

% perform the M+LFBF algorithm
for it = 1:opt.iter
    
    tic;
   
    x_hat = x + gamma * (g{1}.adj_op(w1) + g{2}.adj_op(w2)) ;
    u_hat = u - gamma * w1;    
    v_hat = v - gamma * w2;  
    zeta1_hat = zeta1 - gamma * d_zeta1;     
    zeta2_hat = zeta2 - gamma * d_zeta2;        
    
    p_x = fprojS(x_hat, gamma);
    [p_u, p_epi1] = g{1}.project_epi(u_hat, zeta1_hat);
    [p_v, p_epi2] = f(2).project(v_hat, zeta2_hat); 


    % 2/5: function 'f'

    
    w1 = w1  + tau * (2*p_u - u - g{1}.dir_op(2*p_x-x));
    w2 = w2  + tau * (2*p_v - v - g{2}.dir_op(2*p_x-x));
    d_zeta1 = d_zeta1 + tau*(2*p_epi1 - d_zeta1);
    d_zeta2 = d_zeta2 + tau*(2*p_epi2 - d_zeta2);    
    
    p_w1 = zeros(size(w1));
    p_w2 = zeros(size(w2));  
    

    [p_eta1, p_eta2] = f(1).project(d_zeta1/tau,d_zeta2/tau);

    
    w1 = w1 - tau*p_w1;
    w2 = w2 - tau*p_w2;  
    

    
    d_zeta1 = d_zeta1 - tau*p_eta1;
    d_zeta2 = d_zeta2 - tau*p_eta2;    
    
    

    % 4/4: primal forward step
    x_old = x;
    x = p_x;
    u = p_u;
    v = p_v;
    zeta1 = p_epi1;
    zeta2 = p_epi2;
    

    time(it) = toc;
    
    
    %--- OPERATIONS NOT RELATED TO THE ALGORITHM ---%

Fx = shiftdim(sqrt(sum(g{1}.dir_op(x).^2, 1)),1);
eval_crit = sqrt(2)/10*sum(abs( g{2}.dir_op(x) - opt.noisy(:) )) + 0.1*sum( Fx(:) ) ;
fprintf('crit_real = %3.2f\n',eval_crit);    
    % check the stop criterion
    %if norm(x(:) - x_old(:)) < opt.tol * norm(x_old(:))  % || all(x(:) == x_old(:))
    %   break;
    %end
    
    
    % compute the criterion
    for i = 1:n_fun
        varargout{i}(it) = varargin{i}(x);
    end
    %    fprintf('dCalpha = %3.2f \t dS = %3.2f\n',varargout{1}(it),varargout{2}(it));
    
    
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
for i = 1:n_fun
    varargout{i}(it) = varargin{i}(x);
end

% cut the zeros
for i = 1:n_fun
    varargout{i} = varargout{i}(1:it);
end