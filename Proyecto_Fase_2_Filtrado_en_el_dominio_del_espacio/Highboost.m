function [gxy_out] = Highboost(img,k)
J = double(img);
w = [1 1 1; 1 1 1; 1 1 1];
fxyb = imfilter(J, w, 'replicate', 'conv') / 9;
mask = J - fxyb;
gxy = J + k * mask;
gxy_norm = (gxy - min(gxy(:))) / (max(gxy(:)) - min(gxy(:)));
gxy_out = uint8(gxy_norm * 255);
end