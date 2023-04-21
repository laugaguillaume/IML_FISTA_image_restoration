%%    z = zeropad(x, shape)
%%
%% Zero-pad array `x` to dimensions given by `shape`.  The contents of `x` will
%% be approximately centered in the result.
%%
%% Examples:
%%   zeropad(x, size(y));
%%   zeropad(x, 512);
%%   zeropad(x, [200, 300]);
function z = zeropad(x, shape)
    xshape = size(x);
    rank = numel(xshape);
    if numel(shape) ~= rank
        error('bad number of dimensions');
    end
    diff = shape - xshape;
    if min(diff) < 0
        error('output dimensions must be larger or equal input ones');
    end
    offset = floor(diff/2);
    first = offset + 1;
    last = offset + xshape;
    z = zeros(shape);
    if rank == 1
        i1 = first(1); j1 = last(1);
        z(i1:j1) = x;
    elseif rank == 2
        i1 = first(1); j1 = last(1);
        i2 = first(2); j2 = last(2);
        z(i1:j1,i2:j2) = x;
    elseif rank == 3
        i1 = first(1); j1 = last(1);
        i2 = first(2); j2 = last(2);
        i3 = first(3); j3 = last(3);
        z(i1:j1,i2:j2,i3:j3) = x;
    elseif rank == 4
        i1 = first(1); j1 = last(1);
        i2 = first(2); j2 = last(2);
        i3 = first(3); j3 = last(3);
        i4 = first(4); j4 = last(4);
        z(i1:j1,i2:j2,i3:j3,i4:j4) = x;
    elseif rank == 5
        i1 = first(1); j1 = last(1);
        i2 = first(2); j2 = last(2);
        i3 = first(3); j3 = last(3);
        i4 = first(4); j4 = last(4);
        i5 = first(5); j5 = last(5);
        z(i1:j1,i2:j2,i3:j3,i4:j4,i5:j5) = x;
    elseif rank == 6
        i1 = first(1); j1 = last(1);
        i2 = first(2); j2 = last(2);
        i3 = first(3); j3 = last(3);
        i4 = first(4); j4 = last(4);
        i5 = first(5); j5 = last(5);
        i6 = first(6); j6 = last(6);
        z(i1:j1,i2:j2,i3:j3,i4:j4,i5:j5,i6:j6) = x;
    else
        error('too many dimensions');
    end
end
