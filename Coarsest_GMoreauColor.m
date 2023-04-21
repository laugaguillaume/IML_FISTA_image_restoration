function [X]=Coarsest_GMoreauColor(f,g,gradf,proxg,D,D_adjoint,parameters,X,level,bk)
DX = D(X);
proxgk = proxg(DX,parameters.moreau);
Fx = f(X) + g(proxgk) + 1/(2*parameters.moreau)*norm(DX(:)-proxgk(:))^2;
gradFx=gradf(X)+parameters.moreau^(-1).*D_adjoint(DX-proxgk); 
it_number = matlab.lang.makeValidName(['maxit_level' num2str(level)]); 
max_iteration = parameters.(it_number);
alpha_lip = parameters.alpha_moreau;
for iteration=1:max_iteration
    [alpha,ind,X,proxgk]=backtracking_GradientMoreauColor(alpha_lip,parameters.moreau,X,-gradFx,...
        gradFx,Fx,f,g,proxg,D,D_adjoint,bk);
    DX = D(X);
    if bk == 1
        Fx = f(X) + g(proxgk) + 1/(2*parameters.moreau)*norm(DX(:)-proxgk(:))^2;
    end
    gradFx=gradf(X)+parameters.moreau^(-1).*D_adjoint(DX-proxgk);
    if iteration==max_iteration
        disp('stop per maxit')
        break
    end
end

