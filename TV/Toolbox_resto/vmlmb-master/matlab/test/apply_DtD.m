%%     r = apply_DtD(x);
%%
%% Apply the D'*D operator where D is finite difference along the dimensions
%% of the 2-D array `x`.
%%
function r = apply_DtD(x)
    dims = size(x);
    ndims = numel(dims);
    r = zeros(dims);
    if ndims == 2
        %% compute D1'.D1.x
        n = dims(1);
        if n > 1
            dx = x(2:n,:) - x(1:n-1,:);   % dx <-- D1.x
            r(2:n,:) = dx;
            r(1:n-1,:) = r(1:n-1,:) - dx;  % t <-- D1'.D1.x
        end
        %% compute D2'.D2.x
        n = dims(2);
        if n > 1
            dx = x(:,2:n) - x(:,1:n-1);   % dx <-- D2.x
            t = zeros(size(x));
            t(:,2:n) = dx;
            t(:,1:n-1) = t(:,1:n-1) - dx;
            r = r + t;
        end
    else
        error('only 2-D arrays are supported for now');
    end
end
