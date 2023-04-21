function [alpha,ind]=backtracking(alpha,x,p,f,gradfx)
bmax=20;
c1=1e-4;%norm(p(:))^2;%1e-4;
gamma=0.5;
fx=f(x);

for b=1:bmax
    if f(x+alpha.*p)<=fx+alpha.*c1.*sum(gradfx(:).*p(:))
        ind=1;
        break
    else
        ind=-1;
        alpha=gamma*alpha;
    end
   
end