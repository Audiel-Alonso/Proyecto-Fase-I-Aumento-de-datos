function [img_out] = AdaptativoLocal(img,n)


[M,N] = size(img);

gxy = padarray(img, [(n-1)/2 (n-1)/2], 'replicate');
gxy = double(gxy);
fxy = zeros(M,N);
SigN = var(gxy, 0, 'all');
for x = (n-1)/2+1:M+(n-1)/2
    for y = (n-1)/2+1:N+(n-1)/2
        pixVecl = gxy(x-(n-1)/2:1:x+(n-1)/2,y-(n-1)/2:1:y+(n-1)/2);
        SigSxy = var(pixVecl, 0, 'all');
        if SigN>SigSxy
            gan = 1;
        else 
            gan = (SigN)/(SigSxy);
        end
        fxy(x-(n-1)/2,y-(n-1)/2) = gxy(x,y) -gan*(gxy(x,y)-mean(mean(pixVecl)));
    end
end

img_out = uint8(fxy);

end