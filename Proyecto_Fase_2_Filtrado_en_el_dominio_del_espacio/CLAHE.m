function [im_out] = CLAHE(img, L, cl)

[M,N] = size(img);
tam1 = M/2;
tam2 = N/2;
PARTE1 = img(1:tam1,1:tam2);
PARTE2 = img(tam1+1:M,1:tam2);
PARTE3 = img(1:tam1,tam2+1:N);
PARTE4 = img(tam1+1:M,tam2+1:N);

ci1 = tam1 / 2;                
ci2 = tam1 + (M - tam1) / 2;  
cj1 = tam2 / 2;               
cj2 = tam2 + (N - tam2) / 2; 

T1 = CLAHE_f_aux(PARTE1, L, cl);
T2 = CLAHE_f_aux(PARTE2, L, cl);
T3 = CLAHE_f_aux(PARTE3, L, cl);
T4 = CLAHE_f_aux(PARTE4, L, cl);




i = (1:M)';  
j = (1:N);  

Av = zeros(M, 1);
Av(i >= ci2) = 1; 
idx_mid_v = (i > ci1) & (i < ci2); 
Av(idx_mid_v) = (i(idx_mid_v) - ci1) / (ci2 - ci1);
Bv = 1 - Av;

Ah = zeros(1, N);
Ah(j >= cj2) = 1;
idx_mid_h = (j > cj1) & (j < cj2);
Ah(idx_mid_h) = (j(idx_mid_h) - cj1) / (cj2 - cj1);
Bh = 1 - Ah;

W1 = Bv * Bh;
W2 = Av * Bh;
W3 = Bv * Ah;
W4 = Av * Ah;

val = double(img) + 1;
M1 = T1(val);
M2 = T2(val); 
M3 = T3(val); 
M4 = T4(val); 

im_out = W1.*M1 + W2.*M2 + W3.*M3 + W4.*M4;
im_out = uint8(im_out);


end