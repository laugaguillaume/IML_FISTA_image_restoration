function y = project_S(x,opt)

y = x;
ind = find(opt.patch == 1);
if ~isempty(ind)
    y(ind)=255;
end