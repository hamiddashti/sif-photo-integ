clear
clc
warning('off')
data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
script_dir = pwd;
symsolver_fun();
[obs, spectral] = prepare_inputs(data_dir);
load est_vcmax_vqmax_stat_co2.mat
% William's Woodgate data: CO2 experiment

% Prepare observations
licor = obs.aus_data.licor_lrc;
n_sample = length(obs.aus_data.names);  % Number of leaves
n = length(licor.Qin); % Number of light levels
p = 970;
ppm2ubar     =  (p .*1E-3);


%%
an_obs = licor.A(:,1);
Ci = licor.Ci(:,1);
par = licor.Qin(:,1);

data.Qin = licor.Qin(:,1);
data.Tin = licor.Tleaf(:,1)
data.Cin = licor.Ci(:,1)*ppm2ubar;
data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
v = configure_fun_c3(data);
v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
v.Rds = 0.05;

m_vcmax_static1 = model_fun_c3_vcmax_vqmax(50,est_vcmax_vqmax_stat(2,1),v);
an_sim1 = m_vcmax_static1.An_a*1e6;

data.Cin = repmat(250,n,1)*ppm2ubar;
v = configure_fun_c3(data);
m_vcmax_static2 = model_fun_c3_vcmax_vqmax(50,est_vcmax_vqmax_stat(2,1),v);
an_sim2 = m_vcmax_static2.An_a*1e6;

fig = figure (1);
set(gcf,'units','inch','Position',[50 50 8 4],'color','w');
subplot(1,2,1)
plot(Ci,an_obs,'k',LineWidth=2)
set ( gca, 'xdir', 'reverse' )
xlabel("Ci [ppm]",fontsize=12)
ylabel("An [umol m-2 s1]",fontsize=12)
title("An-Ci")
subplot(1,2,2)
plot(par,an_obs,'k',LineWidth=2)
xlabel('PAR [umol m-2 s-1]',fontsize=12)
ylabel("An [umol m-2 s1]",fontsize=12)
title("An-PAR")
sgtitle("Measurements")
saveas(gcf, './Figures/A_Ci_par.png')

fig = figure (2);
set(gcf,'units','inch','Position',[50 50 8 4],'color','w');
subplot(1,2,1)
plot(par,an_obs,'k',LineWidth=2)
hold on
plot(par,an_sim1,'r',LineWidth=2)
xlabel('PAR [umol m-2 s-1]',fontsize=12)
ylabel("An [umol m-2 s1]",fontsize=12)
title("Simulated with measured Ci")
legend('Observed photosynthesis','Simulated photosynthesis')
subplot(1,2,2)
plot(par,an_obs,'k',LineWidth=2)
hold on
plot(par,an_sim2,'r',LineWidth=2)
xlabel('PAR [umol m-2 s-1]',fontsize=12)
ylabel("An [umol m-2 s1]",fontsize=12)
title("Simulated with Ci fixed at 200 ppm")
legend('Observed photosynthesis','Simulated photosynthesis')
saveas(gcf, './Figures/simulated_An.png')

