function snr = SNR_v2(x_true, x_est)

snr = 0;
for k = 1:size(x_true,3)
    x1 = x_true(:,:,k);
    x2 =  x_est(:,:,k);
    snr = snr + 20 * log10( norm( x1(:) ) / norm( x2(:)-x1(:) ) );
end
snr = snr / size(x_true,3);