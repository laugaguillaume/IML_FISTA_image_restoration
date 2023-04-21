function snr=snrdeblur(y,x)
% x: ground truth
% y: estimate
%x=x(:)./norm(x,'inf');
%y=y(:)./norm(y,'inf');
[n,p,q] = size(x);
if q == 1
    snr=20*log10(norm(x,'fro')/norm(y-x,'fro'));
else
    snr=20*log10(norm(x(:))/norm(y(:)-x(:)));
end
