% Sanity check: The main goal is to check if our conversion of PAR from [umol m-2 s-1] 
% to [w m-2 nm-1] based on planck was a OK: to do so I ran the model with
% measured E from Vilfan 2019 data and our converted PAR 
%% 
% Clean up working environment
clear all; % variables
close all; % figures
clc; % command window

% Set up subdirectories
currdir = pwd;
workdir = fullfile(currdir);
resultsdir = fullfile(currdir,'/outputs/');
if ~exist(resultsdir, 'dir')
    mkdir(resultsdir)
end
cd(resultsdir);
outputname = 'Sanity_check';
mkdir(outputname);
Exdir = strcat(resultsdir,outputname);
cd(workdir);
%% 
% Read in data file and create output directory


% Specify environmental conditions
n = 2400;                                   % Steps in vector
data.Qin = transpose(linspace(1,2400,n));   % PAR, umol PPFD m-2 s-1
data.Tin = repmat(25,n,1);                  % Leaf temperature, C
data.Cin = repmat(200,n,1);                 % Mesophyll CO2, ubar
data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
%% Configure model simulations

% Load default set of parameters
v = configure_fun(data);

% Adjust parameters of interest from default values
v.Abs = 0.85;           % Total leaf absorptance to PAR, mol mol-1
v.beta = 0.52;          % PSII fraction of total absorptance, mol mol-1
v.Ku2 = 0e09;           % Rate constant for exciton sharing at PSII, s-1
v.CB6F = 1.2./1e6;      % Cyt b6f density, mol sites m-2
v.RUB = 27.8./1e6;      % Rubisco density, mol sites m-2
v.eps1 = 0;             % PS I transfer function, mol mol-1

%% Run photosynthesis simulation 

% Run simulation
m = model_fun(v);

%% Run Fluspect-CX with simulated phi-F from JB model
datdir = './data/output/fluspect_output/2015/2019-06-19-1452/';
spectral = define_bands()
% load([datdir 'spectral.mat']);
load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);
load([datdir 'measured.mat']);
iwle = (400:700)-349;

% Randomly select a sample from Vilfan dataset
k=17
leafbio.Cab = leafbio_all.Cab(k);
leafbio.Cca = leafbio_all.Cca(k);
leafbio.Cdm = leafbio_all.Cdm(k);
leafbio.Cw = leafbio_all.Cw(k);
leafbio.Cs = leafbio_all.Cs(k);
leafbio.N = leafbio_all.N(k);
leafbio.Cx = leafbio_all.Cx(k);
leafbio.Cant = leafbio_all.Cant(k);

leafbio.fqe = 0.01;
%Select one leaf and do the prospect simulations
[leafopt_test] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);

E = (measured.E(iwle,k)); % Measured PAR 
Fu_E = leafopt_test.Mb*E;
Fd_E = leafopt_test.Mf*E;

% Now convert the PAR
T = 5777;   %Sun temp [K]
wav = (400:1:700)*10^-9; %Wavelenght [m]
weights = planck_weights(wav,T); %weights based on plancks law
alpha = 1/4.565; %converting par [umol sr-1 m-2] to [w m-2]
% Select one leaf and do the prospect simulations
[leafopt_test] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);
par = 2400; % That is maximum PAR which should be equivalent to measured E
par_w = par*alpha; 
par_nm = par_w.*weights;
par_um = par_nm*1000;

Fu = leafopt_test.Mb*par_um';
Fd = leafopt_test.Mf*par_um';
scatter(Fu_E,Fu,'red')
saveas(gcf,strcat(Exdir,"/scatter_plot_sif.png"))
clf
plot((640:850),Fu_E,'blue'); hold on
plot((640:850),Fu,'r')
saveas(gcf,strcat(Exdir,"/compare_sif.png"))

% 
% ax1 = axes('Position',[0.13 0.58 0.77 0.34]);
% ax1.PositionConstraint = 'outerposition';
