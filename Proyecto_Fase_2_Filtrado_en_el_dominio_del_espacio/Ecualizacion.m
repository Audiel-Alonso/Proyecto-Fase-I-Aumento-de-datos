function [img_e] = Ecualizacion(img,L)
% Ecualizacion de una imagen.bmp
%   Detailed explanation goes here
[M,N]=size(img);
[nk1,rk1]=imhist(img);
pk1 = nk1/(M*N);

sk1 = cumsum(pk1) * (L - 1);

T1 = round(sk1);

img_e = uint8(T1(img+1));

end