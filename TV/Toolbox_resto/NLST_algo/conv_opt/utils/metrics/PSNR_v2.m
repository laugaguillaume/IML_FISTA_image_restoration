function [psnr var] = PSNR_v2(x_true, x_est)

if size(x_true,3) == 3
    x_true = color2gray(x_true);
    x_est  = color2gray(x_est);
end

var = sum( (x_est(:)-x_true(:)).^2 ) / (numel(x_est)-1);
psnr = 10 * log10( 255^2 / var );