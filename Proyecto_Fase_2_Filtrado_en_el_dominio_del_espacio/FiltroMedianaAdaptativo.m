function [img_out] = FiltroMedianaAdaptativo(img, Smax)
[M,N] = size(img);
img_out = double(img);
padM = (Smax-1)/2;
img1A = double(padarray(img, [padM padM], 'replicate'));

for x = padM+1:M+padM
    for y = padM+1:N+padM
        S = 3;
        listo = 0;
        
        while listo == 0 && S <= Smax
            w = (S-1)/2;
            pixVecl = img1A(x-w:1:x+w, y-w:1:y+w);
            
            zmin = min(min(pixVecl));
            zmax = max(max(pixVecl));
            zmed = median(pixVecl(:));
            zxy = img1A(x,y);
            
            if zmin < zmed && zmed < zmax
                if zmin < zxy && zxy < zmax
                    img_out(x-padM, y-padM) = zxy;
                else
                    img_out(x-padM, y-padM) = zmed;
                end
                listo = 1;
            else
                S = S + 2;
                if S > Smax
                    img_out(x-padM, y-padM) = zmed;
                    listo = 1;
                end
            end
        end
    end
end
img_out = uint8(img_out);
end