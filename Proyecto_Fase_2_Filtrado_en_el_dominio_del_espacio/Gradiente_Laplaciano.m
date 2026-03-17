function [IF] = Gradiente_Laplaciano(img,gamma)

J = double(img);



pfx=Gradientex(J);
pfy=Gradientey(J);

mag = sqrt(pfx.^2 + pfy.^2);


w=[1 1 1; 1 1 1; 1 1 1];
mags=imfilter(mag,w,'replicate','conv')/9;

mags= mags / max(mags(:));

fxyb=Laplaciano(J);

c = 1;
Rxy = J + c*fxyb;
mask = Rxy.*mags;
gxy = J + mask;

gxy_norm = (gxy - min(gxy(:))) / (max(gxy(:)) - min(gxy(:)));
gxyg = gxy_norm.^gamma;


IF = uint8(gxyg * 255);








end