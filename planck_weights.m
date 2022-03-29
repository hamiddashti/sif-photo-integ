function I_W = planck_weights(wav,T)
    h = 6.626*10^-34;
    c = 2.998*10^8;
    k = 1.38*10^-23;
    a = 2*h*c^2;
    b = exp((h*c)./(k*T.*wav))-1;
    I = a./((wav.^5).*b);
    I_W = I./sum(I) 
end 
