clear
clc
warning('off')
data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
script_dir = pwd;
symsolver_fun();
[obs, spectral] = prepare_inputs(data_dir);

%% William's Woodgate data: LRC experiment

% Prepare observations
licor_lrc = obs.aus_data.licor_lrc;

n_sample = length(obs.aus_data.names);  % Number of leaves
n = length(licor_lrc.Qin); % Number of light levels
% Prepare model params
i = 5;
an_obs = licor_lrc.A(:,i);
data.Qin = licor_lrc.Qin(:,i);
data.Tin = licor_lrc.Tleaf(:,i);
data.Cin = licor_lrc.Ci(:,i);                 % Mesophyll CO2, ubar
data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
v = configure_fun_c3(data);
v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2

%% Invert only for vcmax 
% Static
vcmax_0 =  10;
vars_0 = vcmax_0;
opts.LBounds = 10; opts.UBounds = 600; 
opts.Restarts=3;
opts.Noise.on=1;
[vars_est_statis,fmin] = cmaes('chi2_vcmax',vars_0,[],opts,v,an_obs);

% Dynamic
v.alpha_opt = 'dynamic';
v.solve_xcs = solve_xcs;
% [vars_est,fmin,counteval,stopflag,out,bestever] = cmaes('chi2P5B',vcmax_0,[],opts,v,an_obs);
[vars_est_dynamic,fmin] = cmaes('chi2_vcmax',vars_0,[],opts,v,an_obs);
[vars_est_statis, vars_est_dynamic]

m_vcmax_static = model_fun_c3_vcmax(vars_est_statis,v);
m_vcmax_dynamic = model_fun_c3_vcmax(vars_est_dynamic,v);

plot(data.Qin, m_vcmax_static.An_a*1e6,'k')
hold on
plot(data.Qin, m_vcmax_dynamic.An_a*1e6,'b')
hold on
plot(data.Qin,an_obs,'r')



%% Invert for vcmax and vqmax 
% Static
vcmax_0 =  10;
vqmax_0 = 175;
vars_0 = [vcmax_0, vqmax_0]';
opts.LBounds = [10, 50]'; opts.UBounds = [600, 400]'; 
opts.Restarts=3;
opts.Noise.on=1;
[vars_est_statis,fmin] = cmaes('chi2_vcmax_vqmax',vars_0,[],opts,v,an_obs);
vars_est_statis
m_vcmax_static = model_fun_c3_vcmax_vqmax(vars_est_statis(1),...
   vars_est_statis(2),v);

% Dynamic
v.alpha_opt = 'dynamic';
v.solve_xcs = solve_xcs;
[vars_est_dynamic,fmin] = cmaes('chi2_vcmax_vqmax',vars_0,[],opts,v,an_obs);
vars_est_dynamic
m_vcmax_dynamic = model_fun_c3_vcmax_vqmax(vars_est_dynamic(1),...
   vars_est_dynamic(2),v);

plot(data.Qin, m_vcmax_static.An_a*1e6,'k')
hold on
plot(data.Qin, m_vcmax_dynamic.An_a*1e6,'b')
hold on
plot(data.Qin,an_obs,'r')

%% Invert for vcmax and ku

% Static
vcmax_0 =  10;
ku2_0 = 2;
vars_0 = [vcmax_0, ku2_0]';
opts.LBounds = [10, 0]'; opts.UBounds = [600, 3]'; 
opts.Restarts=3;
opts.Noise.on=1;
[vars_est_statis,fmin] = cmaes('chi2_vcmax_ku2',vars_0,[],opts,v,an_obs);
vars_est_statis
m_vcmax_static = model_fun_c3_vcmax_ku2(vars_est_statis(1),...
   vars_est_statis(2),v);

% Dynamic
v.alpha_opt = 'dynamic';
v.solve_xcs = solve_xcs;
[vars_est_dynamic,fmin] = cmaes('chi2_vcmax_ku2',vars_0,[],opts,v,an_obs);
vars_est_dynamic

m_vcmax_dynamic = model_fun_c3_vcmax_ku2(vars_est_dynamic(1),...
   vars_est_dynamic(2),v);

plot(data.Qin, m_vcmax_static.An_a*1e6,'k')
hold on
plot(data.Qin, m_vcmax_dynamic.An_a*1e6,'b')
hold on
plot(data.Qin,an_obs,'r')

%% Invert for vcmax, vqmax and ku2

vcmax_0 =  10;
vqmax_0 = 175;
ku2_0 = 2;
vars_0 = [vcmax_0,vqmax_0,ku2_0]';

opts.LBounds = [10, 50, 0.0]'; opts.UBounds = [600, 400, 4]'; 
opts.Restarts=3;
opts.Noise.on=1;
OPTS.LogModulo = 0;
[vars_est_static,fmin] = cmaes('chi2_vcmax_vqmax_ku2',vars_0,[],opts,v,an_obs);
vars_est_static

m_vcmax_staitic = model_fun_c3_vcmax_vqmax_ku2(vars_est_static(1),vars_est_static(2),vars_est_static(3),v);

v.alpha_opt = 'dynamic';
v.solve_xcs = solve_xcs;
[vars_est_dynamic,fmin] = cmaes('chi2_vcmax_vqmax_ku2',vars_0,[],opts,v,an_obs);
vars_est_dynamic
m_vcmax_dynamic = model_fun_c3_vcmax_vqmax_ku2(vars_est_dynamic(1),vars_est_dynamic(2),vars_est_dynamic(3),v);


plot(data.Qin, m_vcmax_staitic.An_a*1e6,'k')
hold on
plot(data.Qin, m_vcmax_dynamic.An_a*1e6,'b')
hold on
plot(data.Qin,an_obs,'r')
