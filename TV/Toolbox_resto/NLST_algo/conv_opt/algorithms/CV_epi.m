 function [x, it, time, varargout] = CV_epi(x_init, f, g, h, opt, varargin)
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

%-----%


% get the number of functions
Ng = length(g{2});

% select the step-size ---------------------------------------------------%
if ~isfield(opt, 'gamma') || isempty(opt.gamma)
    opt.gamma = 1;
end
gamma = opt.gamma;
beta = g{1}.beta;
for m = 1:Ng
    beta = beta + g{2}(m).beta;
end
tau = 0.99/(h.beta/2 + gamma*beta);
fprintf('gamma=%3.2f et tau =%3.2f\n',gamma,tau)

%-------------------------------------------------------------------------%


% initialize the solution
x = x_init;
x2 = g{1}.init;

% initialize the temporary variables
v = cell(1,Ng);
for m = 1:Ng
    v{m} = g{2}(m).dir_op(x);
end
v1 = g{1}.dir_op(x);
v2 = x2;

% perform the M+LFBF algorithm
for it = 1:opt.iter
    
    tic;
   
    
    
    % 1/5: adjoint operators of 'v'
    x_hat  = x  - gamma * ( h.grad(x) + g{1}.adj_op(v1) + compute_adjoint(g{2},v) );
    x2_hat = x2 - gamma * v2;

    % 2/5: function 'f'
    p  = f.prox(x_hat, gamma);
    p2 = g{1}.project_V(x2_hat);
    
    % 3/5: 'v' updating (prox & direct operators of 'v')

    for m = 1:Ng
        v_hat   = v{m}  + tau * g{2}(m).dir_op(2*p-x);
        v{m}    = v_hat - tau * g{2}(m).prox(v_hat/tau, 1/tau);
    end

    v1_hat  = v1 + tau * g{1}.dir_op(2*p - x);
    v2_hat = v2 + tau * (2*p2-x2);
    [epi_1, epi_2] = g{1}.project_epi(v1_hat/tau, v2_hat/tau);
    v1  = v1_hat - tau * epi_1;
    v2 = v2_hat - tau * epi_2;
 
    % 4/4: primal forward step
    x_old = x;
    x2_old = x2;
    x = p;
    x2 = p2;
    

    time(it) = toc;
    
    
    %--- OPERATIONS NOT RELATED TO THE ALGORITHM ---%
    
    
    % check the stop criterion
    if norm(x(:) - x_old(:)) < opt.tol * norm(x_old(:))  % || all(x(:) == x_old(:))
        break;
    end
    
    
    % compute the criterion
    for i = 1:n_fun
        varargout{i}(it) = varargin{i}(x,x2);
    end
    
    
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
    varargout{i}(it) = varargin{i}(x,x2);
end

% cut the zeros
for i = 1:n_fun
    varargout{i} = varargout{i}(1:it);
end