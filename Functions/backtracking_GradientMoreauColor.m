function [alpha,ind,Y,proxgY]=backtracking_GradientMoreauColor(alpha,moreau,X,direction,gradFx,Fx,f,g,proxg,D,D_adjoint,bk)
iter_max=20;
c1=1e-4;
gamma=0.5;
Y = X;
if bk == 1
    scalar = c1.*gradFx(:)'*direction(:);
    for b=1:iter_max
        X = X+alpha.*direction;
        DX = D(X);
        proxgk = proxg(DX,moreau);
        fxalpha = f(X) + g(proxgk) + 1/(2*moreau)*norm(DX(:)-proxgk(:))^2;
        if  fxalpha <= Fx+alpha*scalar 
            ind=1;
            Y = X;
            proxgY = proxgk;
            break
        else
            ind=-1;
            alpha=gamma*alpha;
        end
    end
else
    ind = 1;
    Y = X + alpha.*direction;
    proxgY = proxg(D(Y),moreau); 
end