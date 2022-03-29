function m= main_photo()
data.Qin = 2400;   % PAR, umol PPFD m-2 s-1
data.Tin = 25;                 % Leaf temperature, C
data.Cin = 200;                % Mesophyll CO2, ubar
data.Oin = 209;                 % Atmospheric O2, mbar
v = configure_fun(data);
m = model_fun(v);
end 