%% Step 1: set paths
clear all 
clc
my_dir = pwd;
% my_dir = './SA/safe_R1.1'; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';
datdir = './../../data/output/fluspect_output/2015/2019-06-19-1452/';
% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))
out_dir = strcat(my_dir,'/SA_outputs/')

%% Step 2: setup the model and define input ranges

% initialize the parameter space. We will change the Qin and Tin in the
% loop later

data.Qin = 2400;   % PAR, umol PPFD m-2 s-1
data.Tin = 25;                 % Leaf temperature, C
data.Cin = 200;                % Mesophyll CO2, ubar
data.Oin = 209;                 % Atmospheric O2, mbar
v = configure_fun(data);

% load fluspect data
iwle = (400:750)-349;
datdir = './../../data/output/fluspect_output/2015/2019-06-19-1452/';
load([datdir 'spectral.mat']);
load([datdir 'measured.mat']);
load([datdir 'leafopt.mat']);
load([datdir 'Tsimulated.mat']);
load([datdir 'Rsimulated.mat']);
% load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);

M  = 21 ; % number of uncertain parameters [ Sm beta alfa Rs Rf ]
DistrFun  = 'unif'  ; % Parameter distribution
x_center = [v.Abs,v.beta,v.CB6F,v.RUB,v.Kf,v.Kd,v.Kp1,v.Kn1,v.Kp2,...
    v.kq,v.nl,v.nc,v.Kc,v.ko,40,0.009,0.012,10,1.4,0.75,5]
thresh = x_center.*0.1
thresh_low = x_center-thresh;
thresh_high = x_center+thresh;
DistrPar = {[thresh_low(1) thresh_high(1)];...
    [thresh_low(2) thresh_high(2)];
    [thresh_low(3) thresh_high(3)];
    [thresh_low(4) thresh_high(4)];
    [thresh_low(5) thresh_high(5)];
    [thresh_low(6) thresh_high(6)];
    [thresh_low(7) thresh_high(7)];
    [thresh_low(8) thresh_high(8)];
    [thresh_low(9) thresh_high(9)];
    [thresh_low(10) thresh_high(10)];
    [thresh_low(11) thresh_high(11)];
    [thresh_low(12) thresh_high(12)];
    [thresh_low(13) thresh_high(13)];
    [thresh_low(14) thresh_high(14)];
    [thresh_low(15) thresh_high(15)];
    [thresh_low(16) thresh_high(16)];
    [thresh_low(17) thresh_high(17)];
    [thresh_low(18) thresh_high(18)];
    [thresh_low(19) thresh_high(19)];
    [thresh_low(20) thresh_high(20)];
    [thresh_low(21) thresh_high(21)]};

myfun = 'chi_fun_integ'; 
type = "not_combined";
k=28; % Random measured data from Vilian's dataset
E = measured.E(:,k);
Fu_mes =  interp1(spectral.wlM,measured.Fu(:,k),spectral.wlF)
Fu_mes(spectral.wlF<660) = NaN;
Fd_mes =  interp1(spectral.wlM,measured.Fd(:,k),spectral.wlF)
Fd_mes(spectral.wlF<660) = NaN;

X_Labels = {'Abs','beta','CB6F','RUB','Kf','Kd','Kp1','Kn1','Kp2','Kq',...
    'Knl','Knc','Kc','Ko','Cab','Cw','Cdm','Cca','N','Cx','Cant'} ;
%% Run Sobol
% Comment: the base sample size N is not the actual number of input 
% samples that will be evaluated. In fact, because of the resampling
% strategy, the total number of model evaluations to compute the two
% variance-based indices is equal to N*(M+2) 
SampStrategy = 'lhs'; % Latin Hypercube                                       
N =  20000; % Base sample size.

Qin = [500,1000,1500,2000,2400];
Tin = [20,25,30,35];
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N);
[ XA, XB, XC ] = vbsa_resampling(X);

counter = 1;
for i = 1:length(Qin)
    for j = 1:length(Tin)
        
%         prog = strcat("PAR = ",num2str(data.Qin(i)),"; Temp = ", ...
%             num2str(data.Tin(j)))
        data.Qin = Qin(i);   % PAR, umol PPFD m-2 s-1
        data.Tin = Tin(j);                 % Leaf temperature, C
        data.Cin = 200;                % Mesophyll CO2, ubar
        data.Oin = 209;                 % Atmospheric O2, mbar
        v = configure_fun(data);
        sprintf("PAR = %d and Temp = %d",data.Qin,data.Tin)

        YA = model_evaluation(myfun,XA,v,spectral,optipar,type,E,Fu_mes,Fd_mes) ; % size (N,P)
        YB = model_evaluation(myfun,XB,v,spectral,optipar,type,E,Fu_mes,Fd_mes) ; % size (N,P)
        YC = model_evaluation(myfun,XC,v,spectral,optipar,type,E,Fu_mes,Fd_mes) ; % size (N,P)

        Nboot=500;
        % select the j-th model output:
        j = 1 ; 
        [  Si1(counter,:), STi1(counter,:), Si_sd1(counter,:),...
            STi_sd1(counter,:) ] = vbsa_indices(YA(:,j),YB(:,j),YC(:,j),...
            Nboot);
        j = 2 ;
        [  Si2(counter,:), STi2(counter,:), Si_sd2(counter,:),...
            STi_sd2(counter,:)] = ...
            vbsa_indices(YA(:,j),YB(:,j),YC(:,j),Nboot);
        counter = counter+1;
    end 
end
disp(" ----------- All Done ------------------")
%% ploting sobol indices

% tiledlayout(6,5,'TileSpacing','tight','Padding','tight')

figure (1)
set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
tiledlayout(5,4,'TileSpacing','tight','Padding','compact')

counter = 1;
for i = 1:length(Qin)
    for j = 1:length(Tin)
        nexttile
        boxplot2([Si1(counter,:); STi1(counter,:)],X_Labels)
        legend('main effects','total effects')
        counter =counter+1;
    end
end
% set(gcf,'Position',[100 100 900 600])
% resizefigure(spfig3,5,4,0.1,0.1,0.05,0.05,0.97,0.93)
saveas(gcf,"Fu_fll_SA.png")

figure (2)
set(gcf,'units','inch','Position',[50 50 15 12],'color','w')
tiledlayout(5,4,'TileSpacing','tight','Padding','compact')

counter = 1;
for i = 1:length(Qin)
    for j = 1:length(Tin)
        nexttile
        boxplot2([Si2(counter,:); STi2(counter,:)],X_Labels)
        legend('main effects','total effects')
        counter =counter+1;
    end
end
% set(gcf,'Position',[100 100 900 600])
% resizefigure(spfig3,5,4,0.1,0.1,0.05,0.05,0.97,0.93)
saveas(gcf,"Fd_fll_SA.png")
%% 
save(strcat(out_dir,'SA_Results.mat'))
