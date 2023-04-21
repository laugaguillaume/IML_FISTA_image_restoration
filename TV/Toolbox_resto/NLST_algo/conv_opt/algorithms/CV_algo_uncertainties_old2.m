 function [x, it, time, varargout] = CV_algo_uncertainties(x_init, f, g, h, opt, varargin)
%function [x it time varargout] = LFBF(x_init, f, g, h, opt, varargin)
%
%  Created on: 19/01/12 - Giovanni Chierchia


% check inputs
if nargin < 4
    error('Too few parameters');
end
if nargin < 5 || isempty(opt)
    opt = [];
end
if ~isfield(opt, 'tol') || isempty(opt.tol)
    opt.tol = 10^-4;
end
if ~isfield(opt, 'gamma') || isempty(opt.gamma)
    opt.gamma = 1;
end
if ~isfield(opt, 'iter') || isempty(opt.iter)
    opt.iter = 500;
end
if ~isfield(opt, 'show') || isempty(opt.show)
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
hdl = wait_bar(0, 'Running M+LFBF...');



%---------------------------------%
%----- ALGORITHM STARTS HERE -----%
%---------------------------------%

% default functions
if isempty(f)
    f.prox = @(y,gamma) y;
end
if isempty(g)
    g.prox   = @(y,gamma) y;
    g.dir_op = @(y) y;
    g.adj_op = @(y) y;
    g.beta = 1;
end
if isempty(h) 
    h.grad = @(y) 0;
    h.beta = 0;
end
%-----%


% get the number of functions
Ng = length(g);

% select the step-size ---------------------------------------------------%
if ~isfield(opt, 'gamma') || isempty(opt.gamma)
    opt.gamma = 1;
end
gamma = opt.gamma;
g_beta = 0;
for m = 1:Ng
    g_beta = g_beta + g(m).beta;
end
tau = 0.99/(h.beta/2 + gamma*g_beta);
fprintf('gamma=%3.2f et tau =%3.2f\n',gamma,tau)
%-------------------------------------------------------------------------%


% initialize the solution
x = x_init;
x2 = x;
% initialize the temporary variables
v = cell(1,Ng);
for m = 1:Ng
    v{m} = g(m).dir_op(x);
end
v2 = x2;
t = cell(1,Ng);

% perform the M+LFBF algorithm
for it = 1:opt.iter
    
    tic;
    
    % 1/4: primal forward step
    x_hat = x - tau * ( h.grad(x) + compute_adjoint(g,v) );
    x2_hat = x2 - tau * v2;
    
    % 2/4: primal backward step
    %p = f.prox(x_hat, tau);
    val_inv = inv([(1+tau), -tau;-tau (1+tau)]);
    p = val_inv(1,1)*x_hat + val_inv(1,2)*x2_hat;
    p2 =  val_inv(2,1)*x_hat + val_inv(2,2)*x2_hat ;  
    % 3/4: dual forward-backward-forward step
    for m = 1:Ng
        v_hat   = v{m}  + gamma * g(m).dir_op(2*p - x);
        v{m}    = v_hat - gamma * g(m).prox(v_hat/gamma, 1/gamma);
    end
    v2_hat   = v2  + gamma * (2*p2 - x2);
    v2    = v2_hat - gamma * f.prox(v2_hat/gamma, 1/gamma);
    % 4/4: primal forward step
    x_old = x;
    x = p;
    x2_old = x2;
    x2 = p2;    
    
    time(it) = toc;
    
    %disp(norm(x-x2,'fro'));
    %--- OPERATIONS NOT RELATED TO THE ALGORITHM ---%


    % check the stop criterion
    %if norm(x(:) - x_old(:)) < opt.tol * norm(x_old(:))  % || all(x(:) == x_old(:))
    %    break;
    %end
    
    % compute the criterion
    for i = 1:n_fun
        varargout{i}(it) = varargin{i}(x);
        
    end
    varargout{n_fun+1}(it) = norm(x-x2,'fro');
    %fprintf('crit_real_data = %3.2f\t crit_real_pen = %3.2f\n',varargout{1}(it),varargout{2}(it));
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
varargout{n_fun+1}(it) = norm(x-x2,'fro');

% cut the zeros
for i = 1:n_fun
    varargout{i} = varargout{i}(1:it);
end