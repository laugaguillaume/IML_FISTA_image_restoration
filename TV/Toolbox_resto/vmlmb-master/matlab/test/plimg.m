% A simple function to plot an image.
function plimg(img, fig)
    if nargin >= 2
        figure(fig);
    end
    clf;
    imagesc(img, [0,max(img(:))])
    set(gca, 'YDir', 'normal');
end
