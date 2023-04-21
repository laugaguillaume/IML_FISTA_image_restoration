function [f, g] = define_gf(x,fid_f,reg_f,h,gg)

f =   fid_f(x) + reg_f(x);
g =   h.adj_op(h.grad(h.dir_op(x),1)) + gg.adj_op(gg.grad(gg.dir_op(x),1));
