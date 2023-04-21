function y=edge_function(x,alpha)
% ym=sqrt((1e-4+x.^2));
% y=sum(sum(sum(ym)));
n=sqrt(length(x)/2);
x=reshape(x,n,n,2);
[n,m,~]=size(x); y=0;
for i=1:n
    for j=1:m
        z=[x(i,j,1),x(i,j,2)];
        y=y+sqrt((alpha+norm(z)^2));
    end
end
