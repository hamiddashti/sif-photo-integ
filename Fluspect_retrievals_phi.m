function [optipar]= Fluspect_retrievals_phi(leafbio_all,target,measurement,wl,optipar)

%% get constants
global constants
[constants] = define_constants();

%% the input for FLUSPECT
%load 08-Nov-2016OptiparKcaVZ

%% Define spectral regions
spectral        = define_bands;
nwlP            = length(spectral.wlP);
nwlT            = length(spectral.wlT);
spectral.IwlP   = 1 : nwlP;
spectral.IwlT   = nwlP+1 : nwlP+nwlT;
spectral.wlM    = wl;

P.Z = spectral.wlP>650&spectral.wlP<820;
spectral.Fqe = spectral.wlP(P.Z);
measurement.E      = interp1(spectral.wlM,measurement.E,spectral.wlE);
[~,P.intersect,~] = intersect(spectral.wlF',spectral.Fqe','rows');

%% initial parameter values, and boundary boxes
P.PHI               = optipar.phiII(P.Z) ;
P.CoeffsKcaLength   = spectral.Fqe(1):2:spectral.Fqe(end);
P.Fqe               = length(P.CoeffsKcaLength);
FunFQE              = csaps(spectral.Fqe,P.PHI,1,P.CoeffsKcaLength); %Cubic smoothing spline
params0             = FunFQE';

LB                  = zeros(length(P.CoeffsKcaLength),1);
UB                  = 0.5*ones(length(P.CoeffsKcaLength),1);

%% do the job: fit the model to the data
input               = {leafbio_all,optipar,spectral,target};
f                   = @(params)COST_4Fluspect_Phi(params,measurement,input,P);
tic
params1             = lsqnonlin(f,params0,LB,UB);
toc

phi                 = params1;
wl0 = [640, P.CoeffsKcaLength,850]';
phi = [0;phi;0];

[p] = polyfitB(wl0(1:15)-wl0(1),phi(1:15),3,0);
y = polyval(p,wl0(1:15)-wl0(1));
phi(1:15) = y;

[p] = polyfitB(wl0(end-15:end)-wl0(end),phi(end-15:end),3,0);
y = polyval(p,wl0(end-15:end)-wl0(end));
phi(end-15:end) = y;

phii                = interp1(wl0,phi,spectral.wlF,'spline');
phii                = phii'/sum(phii);
phii                = interp1(spectral.wlF,phii,spectral.wlP,'nearest');
phii(isnan(phii))   = 0;
optipar.phi         = phii';
