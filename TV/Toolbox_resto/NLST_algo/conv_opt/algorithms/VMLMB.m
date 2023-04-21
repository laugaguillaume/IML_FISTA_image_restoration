 function [x, it, time, stepsize, varargout] = VMLMB(x_init, f, g, h, opt, varargin)
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
hdl = wait_bar(0, 'Running EA...');

x = x_init;
gamma = 1.99/(h.beta+ g.beta);

% perform the EA algorithm
for it = 1:opt.iter
    
    tic;
    
    
    
    time(it) = toc;
    stepsize(it) = gamma;
    
    %--- OPERATIONS NOT RELATED TO THE ALGORITHM ---%
    
    % check the stop criterion
    if norm(x(:) - x_old(:)) < opt.tol * norm(x_old(:))  % || all(x(:) == x_old(:))
        break;
    end
    
    
    % compute the criterion

    for i = 1:opt.nbf
        varargout{i}(it) = varargin{i}(x);
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


