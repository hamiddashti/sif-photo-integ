%% Prepare inputs
%calculation of Ci (internal CO2 concentration)
clear 
clc

% RH = min(1, eb./satvap(T-273.15) ); % jak: don't allow "supersaturated" air! (esp. on T curves)
RH = 0.50;
data.Qin = 450;   % PAR, umol PPFD m-2 s-1
data.Tin = 25;                  % Leaf temperature, C
% data.Cin = 200;                 % Atmpspheric CO2, ubar
data.Oin = 209;                 % Atmospheric O2, mbar
BallBerrySlope = 9.0;
BallBerry0 = 0.0100;
minCi = 0.3;
Cs = 200;    %ppm
p = 970;
v = configure_fun(data);
ppm2bar     =  1e-6 .* (p .*1E-3);
% Cs          = Cs .* ppm2bar;
tol = 1e-7;  % 0.1 ppm more-or-less
%%
[Ci,err2] = fixedp_brent_ari(@(x) Ci_next(x, Cs, RH, minCi, BallBerrySlope, BallBerry0, ppm2bar,v), Cs, [], tol); % [] in place of Gamma: it didn't make much difference


function [err,Ci_out] = Ci_next(Ci_in,Cs, RH, minCi, BallBerrySlope, BallBerry0, ppm2bar,v)
% persistent fcount
Ci_in = Ci_in*ppm2bar;    % [ppm --> bar]
%Calculate photosyntesis
m = model_fun_ci(Ci_in,v);
A = m.An_a *1e6 % [mol co2 m-2 s-1 --> umol co2 m-2 s-1]

gs = max(BallBerry0,  (BallBerrySlope.* A .* (RH ./Cs)  + BallBerry0));
Ci_out = max(minCi .* Cs,  Cs - 1.6 * A./gs) ;
Ci_in = Ci_in./ppm2bar;
err = Ci_out-Ci_in
% fcount = fcount + 1; % # of times we called computeA
% m.fcount = fcount;
end
Ci = Ci*ppm2bar;    % [ppm --> bar]
%Calculate photosyntesis
m = model_fun_ci(Ci,v);
m.An_a*1e6

%% Ball Berry Model
function [Ci, gs] = BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input)
%  Cs  : CO2 at leaf surface
%  RH  : relative humidity
%  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
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

end % function

function gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
% add in a bit just to avoid div zero. 1 ppm = 1e-6 (note since A < 0 if Cs ==0, it gives a small gs rather than maximal gs
gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./ (Cs+1e-9)  + BallBerry0);
% clean it up:
%gs( Cs == 0 ) = would need to be max gs here;  % eliminate infinities
gs( isnan(Cs) ) = NaN;  % max(NaN, X) = X  (MATLAB 2013b) so fix it here
end
