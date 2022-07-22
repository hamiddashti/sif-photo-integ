
function vcmax = calc_vcmax_lu(a,b,c,d,chl)
% chl = chl*1e4   % convert from ug cm-2 to ug m-2
R = a.*chl+b      % Rubisco content
R_mol = R*1e6/550000;
kcat = c./(1+d.*R_mol.*8); % Rubisco turnover rate at 25°C
vcmax = kcat*(8/550)*R*1e-3*1e6; %umol co2 m-2 s-1
% vcmax = vcmax .* 1e-6 % mol co2 m-2 s-1 ;
end
