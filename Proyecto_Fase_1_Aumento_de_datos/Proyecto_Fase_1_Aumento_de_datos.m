clear;clc;
imo = imread('86_png.rf.484b890790f01cf4800a8a5059de22fd.jpg');
im = double(imo);
%% Transformacion a escala de grises
img = 0.299*im(:,:,1) + 0.587*im(:,:,2) + 0.114*im(:,:,3);
[M,N] = size(img);

%% a) Volteado
imgV = zeros(M,N);

for i=1:N
    imgV(:,i) = img(:,N+1-i);
end

%% b) Rotacion aplicando interpolacion bilineal para definir la intensidad
%     de la imagen resultante

Desplx = N/2;
Desply = M/2;

ang = 45;
th = ang*pi/180;

matrizTras = [1, 0, Desplx;
        0, 1, Desply;
        0, 0, 1];

matrizInvTras = [1, 0, -Desplx;
        0, 1, -Desply;
        0, 0, 1];
matrizRotacion= [cos(th), -sin(th), 0;
       sin(th), cos(th), 0;
       0, 0, 1];
matrizTF = matrizTras*matrizRotacion*matrizInvTras;

imR = zeros(N, M);
D = zeros(N, M);

for iaster=1:N
    for jaster=1:M
        nuevasCoords=inv(matrizTF)*[iaster jaster 1]';
        ip=fix(nuevasCoords(1)+0.5);
        jp=fix(nuevasCoords(2)+0.5);
        if (ip>0 && ip <N && jp>0 && jp<M)
            D(iaster, jaster) = img(ip,jp);
        end
    end
end

for iaster=1:N
    for jaster=1:M
        nuevasCoords=inv(matrizTF)*[iaster jaster 1]';
        ip=floor(nuevasCoords(1));
        jp=floor(nuevasCoords(2));
        Bv = 1-(nuevasCoords(1)-ip);
        Av = nuevasCoords(1)-ip;
        Bh = 1-(nuevasCoords(2)-jp);
        Ah = nuevasCoords(2)-jp;
        if (ip>0 && ip <N && jp>0 && jp<M)
            imR(iaster, jaster) = Bh*Bv*img(ip,jp) + Ah*Bv*img(ip,jp+1) +Av*Bh*img(ip+1,jp) + Av*Ah*img(ip+1,jp+1);
        end
    end
end

%% c.  Traslación
imT = zeros(N, M);

Desplx1 = N/8;
Desply1 = M/8;

matrizTras1 = [1, 0, Desplx1;
        0, 1, Desply1;
        0, 0, 1];

for iaster=1:N
    for jaster=1:M
        nuevasCoords=matrizTras1*[iaster jaster 1]';
        ip=fix(nuevasCoords(1)+0.5);
        jp=fix(nuevasCoords(2)+0.5);
        if (ip>0 && ip <N && jp>0 && jp<M)
            imT(iaster, jaster) = img(ip,jp);
        end
    end
end

%% d.  Escalamiento aplicando interpolación bilineal para definir la intensidad de la imagen resultante.
sx = 1;
sy = 2;

Ns = ceil(N*sx);
Ms = ceil(M*sy);

imS = zeros(Ns, Ms);

for iaster=1:Ns
    for jaster=1:Ms
        nuevasCoords=[iaster/sx, jaster/sy];
        ip=floor(nuevasCoords(1));
        jp=floor(nuevasCoords(2));
        Bv = 1-(nuevasCoords(1)-ip);
        Av = nuevasCoords(1)-ip;
        Bh = 1-(nuevasCoords(2)-jp);
        Ah = nuevasCoords(2)-jp;
        if (ip>0 && ip <N && jp>0 && jp<M)
            imS(iaster, jaster) = Bh*Bv*img(ip,jp) + Ah*Bv*img(ip,jp+1) +Av*Bh*img(ip+1,jp) + Av*Ah*img(ip+1,jp+1);
        end
    end
end

%% e. Borrado Aleatorio (Random Erase)
imE = img;
p = 0.1;
Ae = M*N*p;
rat = rand()*2 + 0.3;
He = round(sqrt(Ae*rat));
We = round(sqrt(Ae/rat));
if He>=M, He=M-1; end
if We>=N, We=N-1; end
Te = randi(M-He);
Le = randi(N-We);
imE(Te:Te+He, Le:Le+We) = 0;


%% f.   Mezclado de regiones (cutmix)

imM = img;
Hm = round(M*0.3);
Wm = round(N*0.3);
Tm = randi(M-Hm);
Lm = randi(N-Wm);
Ts = randi(M-Hm);
Ls = randi(N-Wm);
imM(Tm:Tm+Hm, Lm:Lm+Wm) = imgV(Ts:Ts+Hm, Ls:Ls+Wm);

%%
figure; 
subplot(7, 3, 1); 
imshow(imo, []);
subplot(7, 3, 2); 
imshow(img, []);

subplot(7, 3, 4); 
imshow(img, []);
subplot(7, 3, 5); 
imshow(imgV, []);

subplot(7, 3, 7); 
imshow(img, []);
subplot(7, 3, 8); 
imshow(D, []);
subplot(7, 3, 9); 
imshow(imR, []);

subplot(7, 3, 10); 
imshow(img, []);
subplot(7, 3, 11); 
imshow(imT, []);

subplot(7, 3, 13); 
imshow(img, []);
subplot(7, 3, 14); 
imshow(imS, []);

subplot(7, 3, 16); 
imshow(img, []);
subplot(7, 3, 17); 
imshow(imE, []);


subplot(7, 3, 19); 
imshow(img, []);
subplot(7, 3, 20); 
imshow(imM, []);


