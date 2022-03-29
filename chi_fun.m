function [chi2] = chi_fun(x,spectral,optipar,type,E,Fu_mes,Fd_mes)
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
