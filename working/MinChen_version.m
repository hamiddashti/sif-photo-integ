%% Prepare inputs
%calculation of Ci (internal CO2 concentration)
clear
clc

Afinal = NaN(11,1);
k = 0;
for irwc = 0:0.1:1
    k = k + 1;
    % RH = min(1, eb./satvap(T-273.15) ); % jak: don't allow "supersaturated" air! (esp. on T curves)
    RH = 0.50;
    data.Qin = 1450;   % PAR, umol PPFD m-2 s-1
    data.Tin = 25;                  % Leaf temperature, C
    % data.Cin = 200;                 % Atmpspheric CO2, ubar
    data.Oin = 209;                 % Atmospheric O2, mbar
    BallBerrySlope = 9.0;
    BallBerry0 = 0.0100;
    minCi = 0.3;
    Ca = 400;    %ppm
    Cs = Ca;
    p = 970;
    v = configure_fun(data);
    ppm2bar     =  1e-6 .* (p .*1E-3);
    % Cs          = Cs .* ppm2bar;
    tol = 1e-7;  % 0.1 ppm more-or-less
    rwc = irwc;
    s=1;
    rwc_t = 0.9;
    rwc_c = 0.4;

    if rwc < rwc_c
        rwc_scalar = 0;
    elseif rwc > rwc_t
        rwc_scalar = 1;
    else
        rwc_scalar = (rwc-rwc_c)/(rwc_t-rwc_c);
    end

    % Unit conversion
    Cs          = Cs .* ppm2bar;
    Ca          = Ca.*ppm2bar;
    Ci_in = Ca * 0.7;

    for i = 1:100
        
        m = model_fun_ci(Ci_in,v);
        A = m.An_a *1e6*rwc_scalar; % [mol co2 m-2 s-1 --> umol co2 m-2 s-1]
        gs_bb = max(BallBerry0, (BallBerrySlope.* (A.*ppm2bar).* (RH ./Cs)  + BallBerry0));
        Ci_in_new = Cs - 1.6/gs_bb*A*ppm2bar;
        
        if abs(Ci_in - Ci_in_new)<1e-100
            break;
        end

        Ci_in = Ci_in_new;
    end

    Afinal(k) = A;
end
figure;plot(Afinal)
% figure;plot(xxx)
%
% mm = model_fun_ci(Ci_in,v);
% A_final = mm.An_a*1e6*ppm2bar;
% gs_final = A_final/(Cs-Ci_in)*1.6;
% [gs_bb, gs_final]