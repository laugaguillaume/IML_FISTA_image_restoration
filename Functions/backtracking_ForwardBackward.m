function [alpha,ind,X_updated]=backtracking_ForwardBackward(alpha,X,direction,f,gradf,proxg,method)
iter_max1=20;
% iter_max2 = 5;
delta = 1e-4;
epsilon = 1e-4;
gamma1=0.5;
fx = f(X);
gradfx = gradf(X);
% alpha = 2.5;
switch method
    case 'Method 1'
        for b=1:iter_max1
            prox_g_X = proxg(X+alpha.*direction,alpha);
            DL = prox_g_X-X;
            if f(prox_g_X) <=fx + epsilon*(gradfx(:)'*DL(:) + 1/(2*alpha)*norm(DL)^2)
                ind=1;
                X_updated = prox_g_X;
                break
            else
                ind=-1;
                alpha=gamma1*alpha;
                X_updated = X;
            end
        end
    case 'Method 2'
        for b=1:iter_max1
            prox_g_X = proxg(X+alpha.*direction,alpha);
            DL = prox_g_X-X;
            gradDL = gradf(prox_g_X)-gradfx;
            if alpha*norm(gradDL(:)) <= delta*norm(DL(:))
                ind=1;
                X_updated = prox_g_X;
                break
            else
                ind=-1;
                alpha=gamma1*alpha;
                if b==iter_max1
                    X_updated = X;
                end
            end
        end
    case 'No'
        X_updated = proxg(X+alpha.*direction,alpha); ind = 1;
end