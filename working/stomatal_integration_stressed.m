%% Prepare inputs
%calculation of Ci (internal CO2 concentration)
clear
clc

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
rwc = 1;
s_g0=2;
s_vcmax=1;
rwc_t = 0.9;
rwc_c = 0.4;

% Unit conversion
Cs          = Cs .* ppm2bar;
Ca          = Ca.*ppm2bar;
Ci_in = Ca.*0.7;

%%
% Apply water stress
rwc_list = 0:0.1:1;

for i =1:length(rwc_list)
    rwc = rwc_list(i);
    BallBerry0_Stressed = water_stress(rwc,BallBerry0,s_g0,rwc_t,rwc_c);
    vcmax = v.RUB.*v.kc;
    vcmax_stressed = water_stress(rwc,vcmax,s_vcmax,rwc_t,rwc_c);
    [Ci,err2] = fixedp_brent_ari(@(x) Ci_next(x, Cs, RH, minCi, BallBerrySlope, BallBerry0_Stressed, ppm2bar,v,vcmax_stressed), Ca.*0.7, [], tol); % [] in place of Gamma: it didn't make much difference

    % %Calculate photosyntesis
    m = model_fun_stressed(Ci,vcmax_stressed,v);
    A_final(i) = m.An_a*1e6;
end
plot(rwc_list,A_final,LineWidth=2)
xlabel("RWC [%]")
set(gca, 'XDir','reverse')
ylabel("An [umol CO2 m-2 s-1")
saveas(gcf,"./Figures/water_stress.png")

% A_stressed = water_stress(rwc,A_final,s,rwc_t,rwc_c);
% disp(A_stressed)

function [err,Ci_out] = Ci_next(Ci_in,Cs, RH, minCi, BallBerrySlope, BallBerry0_Stressed, ppm2bar,v,vcmax_stressed)
% persistent fcount
% Ci_in = Ci_in*ppm2bar;    % [ppm --> bar]
%Calculate photosyntesis

% if rwc < rwc_c
%     rwc_scalar = 0;
% elseif rwc > rwc_t
%     rwc_scalar = 1;
% else
%     rwc_scalar = (rwc-rwc_c)/(rwc_t-rwc_c);
% end


m = model_fun_stressed(Ci_in,vcmax_stressed,v);
A = m.An_a *1e6; % [mol co2 m-2 s-1 --> umol co2 m-2 s-1]

% BallBerry0_stressed = water_stress(rwc,BallBerry0,s,rwc_t,rwc_c);
gs = max(BallBerry0_Stressed, (BallBerrySlope.* (A.*ppm2bar).* (RH ./Cs)  + BallBerry0_Stressed));
Ci_out = max(minCi .* Cs,  Cs - 1.6 * (A.*ppm2bar)./gs);

%Calculate error
err = Ci_out-Ci_in;
% fcount = fcount + 1; % # of times we called computeA
% m.fcount = fcount;


end

%% Ball Berry Model
function [Ci, gs] = BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input)
%  Cs  : CO2 at leaf surface
%  RH  : relative humidity
%  A   : Net assimilation in the same units of CO2 as Cs/m2/s
% BallBerrySlope, BallBerry0,
% minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
% Ci_input : will calculate gs if A is specified.
if nargin > 6 && ~isempty(Ci_input)
    % Ci is given: try and compute gs
    Ci = Ci_input;
    gs = [];
    if ~isempty(A) && nargout > 1
        gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    end
elseif all(BallBerry0 == 0) || isempty(A)
    % EXPLANATION:   *at equilibrium* CO2_in = CO2_out => A = gs(Cs - Ci) [1]
    %  so Ci = Cs - A/gs (at equilibrium)                                 [2]
    %  Ball-Berry suggest: gs = m (A RH)/Cs + b   (also at equilib., see Leuning 1990)
    %  if b = 0 we can rearrange B-B for the second term in [2]:  A/gs = Cs/(m RH)
    %  Substituting into [2]
    %  Ci = Cs - Cs/(m RH) = Cs ( 1- 1/(m RH)  [ the 1.6 converts from CO2- to H2O-diffusion ]
    Ci      = max(minCi .* Cs,  Cs.*(1-1.6./(BallBerrySlope .* RH)));
    gs = [];
else
    %  if b > 0  Ci = Cs( 1 - 1/(m RH + b Cs/A) )
    % if we use Leuning 1990, Ci = Cs - (Cs - Gamma)/(m RH + b(Cs - Gamma)/A)  [see def of Gamma, above]
    % note: the original B-B units are A: umol/m2/s, ci ppm (umol/mol), RH (unitless)
    %   Cs input was ppm but was multiplied by ppm2bar above, so multiply A by ppm2bar to put them on the same scale.
    %  don't let gs go below its minimum value (i.e. when A goes negative)
    gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    Ci = max(minCi .* Cs,  Cs - 1.6 * A./gs) ;
end
end

function gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
% add in a bit just to avoid div zero. 1 ppm = 1e-6 (note since A < 0 if Cs ==0, it gives a small gs rather than maximal gs
gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./ (Cs+1e-9)  + BallBerry0);
% clean it up:
%gs( Cs == 0 ) = would need to be max gs here;  % eliminate infinities
gs( isnan(Cs) ) = NaN;  % max(NaN, X) = X  (MATLAB 2013b) so fix it here
end

