function [alpha,ind,X_updated,PQ]=backtracking_ForwardBackward_Color(alpha,X,direction,gradfx,f,gradf,proxg,method,parameters_reg)
iter_max1=20;
% iter_max2 = 5;
delta = 1e-4;
epsilon = 1e-4;
gamma1=0.5;
%gradfx = gradf(X);
% alpha = 2.5;
switch method
    case 'Method 1'
        fx = f(X);
        for b=1:iter_max1
            [prox_g_X] = proxg(X+alpha.*direction,alpha,parameters_reg);
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
        fx = f(X);
        for b=1:iter_max1
            [prox_g_X] = proxg(X+alpha.*direction,alpha,parameters_reg);
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
        [prox_g_X] = proxg(X+alpha.*direction,alpha,parameters_reg);
        X_updated = prox_g_X; ind = 1;
end