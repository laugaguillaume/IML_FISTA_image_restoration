function [psnr var] = PSNR(x_true, x_est)

var = sum( (x_est(:)-x_true(:)).^2 ) / (numel(x_est)-1);
psnr = 10 * log10( 255^2 / var );