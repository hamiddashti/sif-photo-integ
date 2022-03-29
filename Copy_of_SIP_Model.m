% _______________________________________________________________________
% SIP model
% Version 1.0 (November, 5th 2020)
% _______________________________________________________________________
% for any question or request, please contact:
%
% Dr. Yelu Zeng & Dr. Min Chen
% Joint Global Change Research Institute, Pacific Northwest National Laboratory, College Park, MD 20740, USA
% E-mail: zengyelu@163.com & chenminbnu@gmail.com
%
% Dr. Shengbiao Wu & Dr. Dalei Hao
% E-mail: wushengbiao90@163.com & dalei.hao.93@gmail.com
%
% https://github.com/chenminbnu/SIP-RTM-Leaf
% _______________________________________________________________________
% Function: Simulate leaf spectra properties from 400 nm to 2500 nm with 1 nm interval
%    Input:
%           - Cab = chlorophyll a+b content in ug/cm?
%           - Car = carotenoids content in ug/cm?
%           - Anth = Anthocyanin content in ug/cm?
%           - Cbrown= brown pigments content in arbitrary units
%           - Cw  = equivalent water thickness in g/cm? or cm
%           - Cm  = dry matter content in g/cm?
%   Output:leaf single scattering albedo, reflectance and transmittance 
%
% reference: Wu S.,Zeng Y.,Hao D.,Liu Q.,Li J., Chen X.,R.Asrar G.,Yin G.,
%            Wen J., Yang B., Zhu P.,Chen M.,2020.Quantifying leaf optical
%            properties with spectral invariants theory.
%            DOI: https://doi.org/10.1016/j.rse.2020.112131
%
% Acknowledgement: The authors thank the PROPSECT team for providing the
%                  specific absorption coefficient in calctav.m, dataSpec_PDB.m


% function LRT=SIP_Model(Cab,Car,Ant,Brown,Cw,Cm)
clear all
alpha=600; % constant for the the optimal size of the leaf scattering element 

% input specific absorption coefficient
data    = dataSpec_PDB;
lambda  = data(:,1);    nr      = data(:,2);
Kab     = data(:,3);    Kca    = data(:,4);
Kant    = data(:,5);    Ks  = data(:,6);
Kw      = data(:,7);    Kdm      = data(:,8);

% Some temporary input parameter; lets get it from prospect CX model
datdir = '/home/hamid/SIF/prospect-cx/original/data/output/fluspect_output/2015/2019-06-19-1452/'
load([datdir 'spectral.mat'])
load([datdir 'leafbio.mat'])
load([datdir 'measured.mat']);
load([datdir 'optipar.mat']);


k=4
Cab = leafbio_all.Cab(k);
Car = leafbio_all.Cca(k);
Cm = leafbio_all.Cdm(k);
Cw = leafbio_all.Cw(k);
Brown = leafbio_all.Cs(k);
N = leafbio_all.N(k);
Ant = leafbio_all.Cant(k);

leafbio.Cab = leafbio_all.Cab(k);
leafbio.Cca = leafbio_all.Cca(k);
leafbio.Cdm = leafbio_all.Cdm(k);
leafbio.Cw = leafbio_all.Cw(k);
leafbio.Cs = leafbio_all.Cs(k);
leafbio.N = leafbio_all.N(k);
leafbio.Cx = leafbio_all.Cx(k);
leafbio.Cant = leafbio_all.Cant(k);


Kall    = (leafbio.Cab*Kab+leafbio.Cca*Kca+leafbio.Cant*Kant+...
    leafbio.Cs*Ks+leafbio.Cw*Kw+leafbio.Cdm*Kdm)/(leafbio.Cdm*alpha);
w0=exp(-Kall);


% spectral invariant parameters
fLMA=2765.0*leafbio.Cdm;
gLMA=102.8*leafbio.Cdm;
p=1-(1-exp(-fLMA))./(fLMA+eps);
q=1-2*exp(-gLMA);
qabs=(q^2)^0.5;

% leaf single scattering albedo
w=w0.*(1-p)./(1-p.*w0+eps); 
a = w0.*(1-p)
b = (1-p.*w0+eps)

plot(p.*w0)
% leaf reflectance and leaf transmittance
refl=w.*(1/2+q/2.*(1-p.*w0)./(1-qabs.*p.*w0));
tran=w.*(1/2-q/2.*(1-p.*w0)./(1-qabs.*p.*w0));

LRT     = [lambda w refl tran];

plot(LRT(:,1),LRT(:,3)); hold on 
plot(spectral.wlM(51:end),measured.refl(51:end,k),'r')

%% run fluspect 
leafbio.fqe = 0.01;

[leafopt] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);
E = (measured.E(:,k))';
leafopt_all.Fu(:,k) = leafopt.Mb * E(51:401)';
leafopt_all.Fd(:,k)         = leafopt.Mf * E(51:401)';
leafopt_all.F(:,k) = leafopt_all.Fu(:,k) + leafopt_all.Fd(:,k)  
hold on
plot(leafopt_all.Fu(:,k),'black')
plot(leafopt_all.Fd(:,k),'r')
plot(leafopt_all.F(:,k),'b')
hold off

%% incorporate F into leaf-sip 
wle         = spectral.wlE';    % excitation wavelengths, transpose to column
wlf         = spectral.wlF';    % fluorescence wavelengths, transpose to column
wlp         = spectral.wlP;     % PROSPECT wavelengths, kept as a row vector

minwle      = min(wle);
maxwle      = max(wle);
minwlf      = min(wlf);
maxwlf      = max(wlf);

% indices of wle and wlf within wlp

Iwle        = find(wlp>=minwle & wlp<=maxwle);
Iwlf        = find(wlp>=minwlf & wlp<=maxwlf);

Kall        = (leafbio.Cab*optipar.Kab + leafbio.Cca*optipar.Kca +...
    leafbio.Cdm*optipar.Kdm + leafbio.Cw*optipar.Kw  + leafbio.Cs*optipar.Ks + ... 
    leafbio.Cant*optipar.Kant)/leafbio.N;
j           = find(Kall>0);  
kChlrel(j)  = leafbio.Cab*optipar.Kab(j)./(Kall(j)*N);

sigmoid     = 1./(1+exp(-wlf/10)*exp(wle'/10));  % matrix computed as an outproduct
[Mf, Mb] = deal(leafbio.fqe(1) * ((optipar.phi(Iwlf))*eps) * kChlrel(Iwle).*sigmoid);

size(Mf)
E = (measured.E(51:end,k))';
leafopt_all.Fu (:,k)         = leafopt.Mb * E(51:401)';
leafopt_all.Fd (:,k)         = leafopt.Mf * E(51:401)';


% 
% E = measured.E(50:400,k);
% Etot_absorbed = 1E-3*( (1-leafopt.refl(1:351)-leafopt.tran(1:351)) .* leafopt.kChlrel(1:351))'*E;
% Fem = 1E3*optipar.phi((640:850)-399)*Etot_absorbed*leafbio.fqe; % without considering any absorption in leaf

size(E)
size(wlp)

