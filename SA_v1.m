%% Step 1: set paths
clear all 
clc

my_dir = './SA/safe_R1.1'; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';
datdir = './../../data/output/fluspect_output/2015/2019-06-19-1452/';
% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2: setup the model and define input ranges

% load data

spectral = define_bands;

iwle = (400:750)-349;
datdir = './../../data/output/fluspect_output/2015/2019-06-19-1452/';
load([datdir 'spectral.mat']);
load([datdir 'measured.mat']);
load([datdir 'leafopt.mat']);
load([datdir 'Tsimulated.mat']);
load([datdir 'Rsimulated.mat']);
% load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);
% Define input distribution and ranges:
M  = 10 ; % number of uncertain parameters [ Sm beta alfa Rs Rf ]
DistrFun  = 'unif'  ; % Parameter distribution
DistrPar  = { [0 100]; [0 0.4];[0 0.5];[0 30];[1 4];[0 1.5]; [0 10];...
    [0 0.2]; [0 0.2];[0 0.6] } ; % Parameter ranges
%% Do SA
% Sample parameter space using the resampling strategy proposed by 
% (Saltelli, 2008; for reference and more details, see help of functions
% vbsa_resampling and vbsa_indices) 
SampStrategy = 'lhs' ;
N = 10500 ; % Base sample size.
myfun = 'chi_fun'; 
type = "not_combined";
k=28; % Random measured data from Vilian's dataset
E = measured.E(:,k);
Fu_mes =  interp1(spectral.wlM,measured.Fu(:,k),spectral.wlF)
Fu_mes(spectral.wlF<660) = NaN;
Fd_mes =  interp1(spectral.wlM,measured.Fd(:,k),spectral.wlF)
Fd_mes(spectral.wlF<660) = NaN;
% Comment: the base sample size N is not the actual number of input 
% samples that will be evaluated. In fact, because of the resampling
% strategy, the total number of model evaluations to compute the two
% variance-based indices is equal to N*(M+2) 
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,2*N);
[ XA, XB, XC ] = vbsa_resampling(X);

YA = model_evaluation(myfun,XA,spectral,optipar,type,E,Fu_mes,Fd_mes) ; % size (N,P)
YB = model_evaluation(myfun,XB,spectral,optipar,type,E,Fu_mes,Fd_mes) ; % size (N,P)
YC = model_evaluation(myfun,XC,spectral,optipar,type,E,Fu_mes,Fd_mes) ; % size (N,P)


% select the j-th model output:
j = 1 ; 
[ Si1, STi1 ] = vbsa_indices(YA(:,j),YB(:,j),YC(:,j));
j = 2 ;
[ Si2, STi2 ] = vbsa_indices(YA(:,j),YB(:,j),YC(:,j));

% Compare boxplots:
X_Labels = {'Cab','Cw','Cdm','Cca','N','Cx','Cant','fqeI','fqeII','Cs'} ;
figure
subplot(121)
boxplot2([Si1; STi1],X_Labels)
title('Fu')
subplot(122)
boxplot2([Si2; STi2],X_Labels)
legend('main effects','total effects')
title('Fd')

figure
subplot(121)
stackedbar([Si1; Si2],[],'main effects',[],{'Fu','Fd'})
subplot(122)
stackedbar([STi1; STi2],X_Labels,'total effects',[],{'Fd','fu'})

% Compute confidence bounds:
Nboot = 500 ;
[ Si, STi, Si_sd, STi_sd, Si_lb, STi_lb, Si_ub, STi_ub ] = vbsa_indices(YA(:,1),YB(:,1),YC(:,1),Nboot);
% Plot:
figure % plot both in one plot:
boxplot2([Si; STi],X_Labels,[ Si_lb; STi_lb ],[ Si_ub; STi_ub ])
legend('main effects','total effects')

[ Si, STi, Si_sd, STi_sd, Si_lb, STi_lb, Si_ub, STi_ub ] = vbsa_indices(YA(:,2),YB(:,2),YC(:,2),Nboot);
% Plot:
figure % plot both in one plot:
boxplot2([Si; STi],X_Labels,[ Si_lb; STi_lb ],[ Si_ub; STi_ub ])
legend('main effects','total effects')
