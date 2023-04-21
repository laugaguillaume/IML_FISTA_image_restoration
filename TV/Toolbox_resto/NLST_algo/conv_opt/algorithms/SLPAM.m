 function [u, it, time, stepsize,cost] = SLPAM(x_init, f, g, h, opt)
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
hdl = wait_bar(0, 'Running SLPAM...');

x = x_init;

% descent parameters
norm_D     = g.beta;
lipch_S_du = opt.bet*norm_D + 1e-16;

ck = 1.01*lipch_S_du;
dk = ck/1000;

switch opt.edges
    case 'similar',  S.D_sq = @(u) sum(g.dir_op(u).^2,2);dim_e = 1; e = sum(g.dir_op(x),2);
    case 'distinct', S.D_sq = @(u) g.dir_op(u).^2;dim_e = size(x_init,3);e = g.dir_op(x);
end

S.du = @(u,e)  2*g.adj_op(repmat((1-e).^2,[1 size(x_init,3)/dim_e 1 1]).*g.dir_op(u));
S.de = @(u,e) -2*(1-e).*S.D_sq(u);

S.expr = @(u,e) sum(sum(sum(sum((repmat(1-e,[1 size(x_init,3)/dim_e 1 1]).*g.dir_op(u)).^2))));

J = @(u,e) f.func(u) + opt.bet*S.expr(u,e) + opt.lam*sum(abs(e(:)));

% perform the EA algorithm
for it = 1:opt.iter
    
    tic;
    tmp = x - (opt.bet/ck)*S.du(x,e);
    x2 = reshape(h.prox(tmp(:),1/ck),size(x));
     
    % update of e
    num = opt.bet*S.D_sq(x2) + dk/2.*e; 
    den = opt.bet*S.D_sq(x2) + dk/2.; 
    e  = g.prox(num./den,opt.lam./(2*den));    
    
    x_old = x;
    x = x2;
 
    
    time(it) = toc;
    stepsize(1,it) = ck;
    stepsize(2,it) = dk;
    cost(it) = J(x,e);
    
    %--- OPERATIONS NOT RELATED TO THE ALGORITHM ---%
    
    % check the stop criterion
    if norm(x(:) - x_old(:)) < opt.tol * norm(x_old(:))  % || all(x(:) == x_old(:))
        break;
    end
    
    
    % compute the criterion
% 
%     for i = 1:opt.nbf
%         varargout{i}(it) = varargin{i}(x);
%     end
    
    
    % show the current solution
    if ~isempty(opt.show) && mod(it,10) == 0
        set(0,'CurrentFigure',fig);
        opt.show(x);
        title( sprintf('Iteration %d', it) );
    end
    
    
    % update the wait bar
    wait_bar( it/opt.iter, hdl );
end

u{1} = x;
for k=1:size(x,3);
    u{2}(:,:,k,1) = e(1,k,:,:);
    u{2}(:,:,k,2) = e(2,k,:,:);
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

% % compute the criterion
% for i = 1:opt.nbf
%     varargout{i}(it) = varargin{i}(x);
% end

% % cut the zeros
% for i = 1:opt.nbf
%     varargout{i} = varargout{i}(1:it);
% end