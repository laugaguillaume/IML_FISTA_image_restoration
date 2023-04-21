 function [x, it, time, stepsize, varargout] = EA_with_BTBB(x_init, f, g, h, opt, varargin)
%function [x it time varargout] = EA_with_BT(x_init, f, g, h, opt, varargin)
%
% Created on: 19/01/12 - Giovanni Chierchia
%
% Descente de gradient avec backtracking
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


% open a wait bar
hdl = wait_bar(0, 'Running EA...');

x = x_init;

% perform the EA algorithm with backwtracking
fcost = @(x) varargin{1}(x) + varargin{2}(x);

for it = 1:opt.iter
    
    tic;
    
    gradfx = h.adj_op(h.grad(h.dir_op(x),1)) + g.adj_op(g.grad(g.dir_op(x),1));
    s = -gradfx;
    
    if it==1
        gradfxm1 = gradfx;
        [alpha,ind]=backtracking(1,x,s,fcost,gradfx);
    else
        % Barzilai Borwein
        alpha=alpham1*(sum(gradfxm1(:).*sm1(:)))/(sum(gradfx(:).*s(:)));
        [alpha,ind]=backtracking(alpha,x,s,fcost,gradfx);
        gradfxm1=gradfx;
    end
    

    if ind==1
        x_old = x;
        x=x+alpha*s;
    else
        break
    end
    
    sm1     = s; 
    alpham1 = alpha;
    
    time(it) = toc;
    
    
    %--- OPERATIONS NOT RELATED TO THE ALGORITHM ---%
    
    % check the stop criterion
    if norm(x(:) - x_old(:)) < opt.tol * norm(x_old(:))  % || all(x(:) == x_old(:))
        break;
    end
    
    
    % compute the criterion

    for i = 1:opt.nbf
        varargout{i}(it) = varargin{i}(x);
    end
    stepsize(it) = alpha;
    
    
    
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
% 
% 
% 
% close the wait bar
 close(hdl);

% close the figure
if ~isempty(opt.show)
   close(fig);
end

% compute the times
time = cumsum(time(1:it));

% compute the criterion
for i = 1:opt.nbf
    varargout{i}(it) = varargin{i}(x);
end

% cut the zeros
for i = 1:opt.nbf
    varargout{i} = varargout{i}(1:it);
end