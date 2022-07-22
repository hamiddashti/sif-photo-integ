% Simulating case 9 of the Jen's paper: 
% Clean up working environment
clear all; % variables
close all; % figures
clc; % command window

% Set up subdirectories
currdir = pwd;
workdir = fullfile(currdir);
resultsdir = fullfile(currdir,'/outputs/');
mkdir(resultsdir)
cd(resultsdir);
outputname = 'Ex2';
mkdir(outputname);
Exdir = strcat(resultsdir,outputname)
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

E = (measured.E(iwle,k));
Fu_E = leafopt_test.Mb*E;
Fd_E = leafopt_test.Mf*E;

T = 5777;   %Sun temp [K]
wav = (400:1:700)*10^-9; %Wavelenght [m]
weights = planck_weights(wav,T); %weights based on plancks law
alpha = 1/4.565; %converting par [umol sr-1 m-2] to [w m-2]

map =  jet(2400);
% for i=1:length(m.phi2F_a)
fig = figure(1)
for i=1:10
    disp(i)
    leafbio.fqe = m.phi2F_a(i);
    % Select one leaf and do the prospect simulations
    [leafopt_test] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);

    par = m.Q(i)*10^6; 
    par_w = par*alpha; 
    par_nm = par_w.*weights;
    par_um = par_nm*1000;

    Fu = leafopt_test.Mb*par_um';
    Fd = leafopt_test.Mf*par_um';
    
    subplot(2,1,1)
    h = plot((640:850),Fu);
    h.Color = map(i,:);
    hold on
    
    subplot(2,1,2)
    h2 = plot((640:850),Fd);
    h2.Color = map(i,:);
    hold on
end

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 6])
subplot(2,1,1)
title("F in backward direction",'FontSize', 8)
xlabel('wl (nm)','FontSize', 10)
ylabel('L_F (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 10)
subplot(2,1,2)
title("F in forward direction",'FontSize', 8)
xlabel('wl (nm)','FontSize', 10)
ylabel('L_F (W m^{-2}\mum^{-1}sr^{-1})','FontSize', 10)
h = axes(
,'visible','off'); 
c = colorbar(h,'Position',[0.93 0.168 0.02 0.7]);  % attach colorbar to h
colormap(c,'jet')
caxis(h,[0,2400]);             
c.Label.String = 'PAR [umol PPFD (m^{-2} s^{-1})]';

c.Layout.Tile = 'east';

saveas(gcf,strcat(Exdir,"/test1.png"))

% 
% ax1 = axes('Position',[0.13 0.58 0.77 0.34]);
% ax1.PositionConstraint = 'outerposition';
