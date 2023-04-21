%% Function to compute the objective function and its gradient.
function [f, g] = fg(x)
    global mu dat mtf wgt;
    H = @(x) real(ifft2(mtf.*fft2(x)));
    Ht = @(x) real(ifft2(conj(mtf).*fft2(x)));
    W = @(x) wgt.*x;
    lhs = @(x) Ht(W(H(x))) + mu*apply_DtD(x);

    r = H(x) - dat;
    Wr = W(r);
    muDtDx = mu*apply_DtD(x);
    g = 2*(Ht(Wr) + muDtDx);
    f = optm_inner(r, Wr) + optm_inner(x, muDtDx);
end
