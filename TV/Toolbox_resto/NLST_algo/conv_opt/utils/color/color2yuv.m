function yuv = color2yuv(rgb,idx)

if nargin == 1 || isempty(idx)
    idx = [1 2 3];
end

if size(rgb,3) == 1
    yuv = rgb;
else
    T = [65.481 128.553 24.966; -37.797 -74.203 112; 112 -93.786 -18.214] / 255;
    offset = [16; 128; 128];

    for p = idx
        yuv(:,:,p) = T(p,1) * rgb(:,:,1) + T(p,2) * rgb(:,:,2) + T(p,3) * rgb(:,:,3) + offset(p);
    end
end