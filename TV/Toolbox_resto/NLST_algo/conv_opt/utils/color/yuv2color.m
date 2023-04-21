function rgb = yuv2color(yuv,idx)

if nargin == 1 || isempty(idx)
    idx = [1 2 3];
end

if size(yuv,3) == 1
    rgb = yuv;
else

    T = [65.481 128.553 24.966; -37.797 -74.203 112; 112 -93.786 -18.214] / 255;
    T = T^-1;
    offset = T * [16; 128; 128];

    for p = idx
        rgb(:,:,p) = T(p,1) * yuv(:,:,1) + T(p,2) * yuv(:,:,2) + T(p,3) * yuv(:,:,3) - offset(p);
    end
    
    rgb = min(max(rgb,0.0),255);
end