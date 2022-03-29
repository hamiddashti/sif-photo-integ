function [chi2] = chi_fun(x,spectral,optipar,type,E,Fu_mes,Fd_mes)

% phys variables 
v.Abs = x(1);           % Total leaf absorptance to PAR, mol mol-1
v.beta = x(2);          % PSII fraction of total absorptance, mol mol-1
v.CB6F = x(3);
v.RUB = x(4);

% photochemical vars 
v.Kf = x(5);
v.Kd = x(6);
v.Kp1 = x(7);
v.Kn1 = x(8);
v.Kp2 = x(9);


v.Ku2 = x(3);           % Rate constant for exciton sharing at PSII, s-1
v.CB6F = x(4);      % Cyt b6f density, mol sites m-2
v.RUB = x(5);      % Rubisco density, mol sites m-2
v.eps1 = x(6);             % PS I transfer function, mol mol-1




leafbio.Cab= x(1);
leafbio.Cw= x(2);
leafbio.Cdm= x(3);
leafbio.Cca= x(4);
leafbio.N= x(5);
leafbio.Cx= x(6);
leafbio.Cant= x(7);
leafbio.fqe(1)= x(8);
leafbio.fqe(2)= x(9);
leafbio.Cs = x(10); 
leafopt = main_flu(spectral,leafbio,optipar,type,E);
chi2(1) = sqrt(sum((Fu_mes-leafopt.Fu').^2,'omitnan'));
chi2(2) = sqrt(sum((Fd_mes-leafopt.Fd').^2,'omitnan'));
end
