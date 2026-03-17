clear;clc;
imo = imread('42_png.rf.01629717e50bd7318cebd6424442ee2d.jpg');
im = double(imo);
%% Transformacion a escala de grises
imgg = 0.299*im(:,:,1) + 0.587*im(:,:,2) + 0.114*im(:,:,3);
img = uint8(imgg);
img_r = imnoise(img, "salt & pepper", 0.03);
L = 2^8;
%% 1)Ecualizacion del histograma

img_E = Ecualizacion(img,L);

%% 2)  CLAHE
cl = 70;
img_CLAHE = CLAHE(img,L,cl);

%% 3)  Highboost
k = 2;
img_HG = Highboost(img,k);

%% 4)  Gradiente-Laplaciano

gamma = 0.7;
img_GL = Gradiente_Laplaciano(img,gamma);


%% 5)  Filtro Adaptativo Local


n = 5;
img_FL = AdaptativoLocal(img_r,n);

%% 5)  Adaptative Median Filter

Smax = 5;
img_AMF = FiltroMedianaAdaptativo(img_r, Smax);


figure; 
subplot(4, 4, 1); 
imshow(imo, []);
subplot(4, 4, 2); 
imshow(img, []);

subplot(4, 4, 3); 
imshow(img, []);
subplot(4, 4, 4); 
imshow(img_E, []);

subplot(4, 4, 5); 
imshow(img, []);
subplot(4, 4, 6); 
imshow(img_CLAHE, []);


subplot(4, 4, 7); 
imshow(img, []);
subplot(4, 4, 8); 
imshow(img_HG, []);

subplot(4, 4, 9); 
imshow(img, []);
subplot(4, 4, 10); 
imshow(img_GL, []);

subplot(4, 4, 11); 
imshow(img_r, []);
subplot(4, 4, 12); 
imshow(img_FL, []);


subplot(4, 4, 13); 
imshow(img_r, []);
subplot(4, 4,14); 
imshow(img_AMF, []);


