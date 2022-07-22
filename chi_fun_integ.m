function [chi2] = chi_fun_integ(x,v,spectral,optipar,type,E,Fu_mes,Fd_mes)

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

% Biochem constant
v.Kq = x(10);
v.nl = x(11);
v.nc = x(12);
v.Kc = x(13);
v.ko = x(14);

% fluspect params
leafbio.Cab= x(15);
leafbio.Cw= x(16);
leafbio.Cdm= x(17);
leafbio.Cca= x(18);
leafbio.N= x(19);
leafbio.Cx= x(20);
leafbio.Cant= x(21);
leafbio.Cs= 0;

m = model_fun(v);
if m.phi1F_a <0 | m.phi2F_a <0
chi2(1) = 0;
chi2(2) = 0;
return
end

if isnan(m.phi1F_a) | isnan(m.phi2F_a) 
chi2(1) = 0;
chi2(2) = 0;
return
end

leafbio.fqe(1)= m.phi1F_a;
leafbio.fqe(2)= m.phi2F_a;

leafopt = main_flu(spectral,leafbio,optipar,type,E);
chi2(1) = sqrt(sum((Fu_mes-leafopt.Fu').^2,'omitnan'));
chi2(2) = sqrt(sum((Fd_mes-leafopt.Fd').^2,'omitnan'));
end
