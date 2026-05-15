function Resultado = FLBP(img)
J = double(img);

J_pad = padarray(J, [1 1], 'replicate');


TL = J_pad(1:end-2, 1:end-2); 
TM = J_pad(1:end-2, 2:end-1); 
TR = J_pad(1:end-2, 3:end);  
R  = J_pad(2:end-1, 3:end);  
BR = J_pad(3:end,   3:end);   
BM = J_pad(3:end,   2:end-1); 
BL = J_pad(3:end,   1:end-2);
ML = J_pad(2:end-1, 1:end-2); 


b7 = TL >= J;
b6 = TM >= J;
b5 = TR >= J;
b4 = R  >= J;
b3 = BR >= J;
b2 = BM >= J;
b1 = BL >= J;
b0 = ML >= J;


LBP = b7*(128) + b6*(64) + b5*(32) + b4*(16) + b3*(8) + b2*(4) + b1*(2) + b0*(1);

Resultado = uint8(LBP);
end