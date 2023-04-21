function [X,X_iteration,fun,time_steps,snr_steps]=StandardForwardbackward_Color(f,g,gradf,proxg,parameters_reg,parameters,X0,Xref,level,method)

X=X0;
gradfx=gradf(X); 
% vec_gradfxm1=vec_gradfx; 
it_number = matlab.lang.makeValidName(['maxit_level' num2str(level)]); 
max_iteration = parameters.(it_number);
time_steps = NaN(1,max_iteration);
time_steps(1) = 0;
snr_steps = NaN(1,max_iteration);
snr_steps(1) = snrdeblur(X0,Xref);
fun = NaN(2,max_iteration);
fun(:,1)=[f(X);g(X)];
alpha_index = matlab.lang.makeValidName(['level' num2str(level) 'alpha']);
alpha = parameters.(alpha_index);
for iteration=1:max_iteration
    %pars.MAXITER = parameters.init_subit + 1;
    tic
    [alpha,ind,X]=backtracking_ForwardBackward_Color(alpha,X,-gradfx,gradfx,f,gradf,proxg,method,parameters_reg);
    gradfx=gradf(X);
    fun(:,iteration+1)=[f(X);g(X)];
%     if abs(sum(fun(:,iteration))-sum(fun(:,iteration+1)))/sum(fun(:,iteration+1))<pars.epsilon
%             pars.epsilon = pars.epsilon/10;
%     end
    if sum(fun(:,iteration))<sum(fun(:,iteration+1))
            parameters_reg.epsilon = parameters_reg.epsilon/10;
    end
    time_steps(iteration+1) = toc + time_steps(iteration);
    %snr_steps(iteration+1) = snrdeblur(X,Xref);
    snr_steps(iteration+1) = norm(X(:)-Xref(:));
    %stopping criterions
    if iteration == parameters.p
        X_iteration = X;
    end
    if iteration==max_iteration
        disp('stop per maxit')
        break
    end
end

end
