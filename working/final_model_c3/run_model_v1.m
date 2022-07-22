% Running the model in the forward mode
clear
clc
warning('off')
data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
script_dir = pwd;
% [obs, spectral] = prepare_inputs(data_dir);

load optipar.mat
load chems_inversion.mat   % Estimated chemistry from invert_chem_v1.m
load est_vcmax_vqmax_stat_co2.mat
load obs
load spectral

obs_ref = obs.aus_data.asd.refl;
obs_tran = obs.aus_data.asd.trans;
obs_abs = obs.aus_data.asd.abs;

wl_low = 400;
wl_high = 2000;
wl_obs = spectral.wl_ASD;
wlp = spectral.wlP;
I_obs = find(wl_obs>=wl_low & wl_obs<=wl_high); % data after 2000 nm is noisy
I_wlp = find(wlp>=wl_low & wlp<=wl_high); % data after 2000 nm is noisy
wlPAR = spectral.wlPAR';
IparP    = find(wlp>=400 & wlp<=700); % Indices for PAR wavelenghts within wl


%% Now run the photosynthesis model for different co2 levels

licor = obs.aus_data.licor_co2;
n_sample = length(obs.aus_data.names);  % Number of leaves
n = length(licor.Ci); % Number of Co2 levels
symsolver_fun();
    
for i = 4:4
    
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
%     v.alpha_opt = 'dynamic';
%     v.solve_xcs = solve_xcs;
    m = model_fun_c3_vcmax_vqmax(est_vcmax_vqmax_stat(1,i)...
        ,est_vcmax_vqmax_stat(2,i),v);
    fm = m.Fm_a./m.Fo_a;
    fmp = m.Fmp_a./m.Fo_a;
    npq = (fm./fmp)-1;
    
    plot(npq)
    hold on
    plot(m.PAM9_a)
end
plot(m.Kn2_a./1e9)
plot(licor.NPQ)
%% Plot Nasstasia Fs

vilfan = obs.vilfan_data.data_co2;
leaf_id = vilfan.leaf_id;
Ci = vilfan.licor_pam.Ci;
ft = vilfan.licor_pam.Ft;
fo = vilfan.licor_pam.Fo;
ft_fo = ft./fo;


for i =1:20
hold on
% plot(ci(leaf_id==i),ft_fo(leaf_id==i))
plot(ci(leaf_id==i),fo(leaf_id==i))
end

%%

fig = figure (1);
set(gcf,'units','inch','Position',[50 50 4 0.5],'color','w');
subplot(1,3,1)
plot(licor.Ci(:,i),m_dynamic.Fs_a)
xlabel("Ci",fontsize=12)
ylabel("Simulated F_s",fontsize=12)

subplot(1,3,2)
plot(licor.Ci(:,i),licor.Fs(:,i))
xlabel("Ci",fontsize=12)
ylabel("Measured F_s",fontsize=12)

subplot(1,3,3)
scatter(m_dynamic.Fs_a,licor.Fs(:,i),"*")
xlabel("Simulated",fontsize=12)
ylabel("Measured",fontsize=12)

saveas(gcf,"./Figures/test.png")






plot(m_dynamic.An_a*1e6)
hold on
plot(an_obs)




plot(licor.Ci(:,i),m_dynamic.PAM9_a)
hold on
plot(licor.Ci(:,i),obs_NPQ,'r')



phi2F = m_vcmax_dynamic.phi2F_a
NPQ = m_vcmax_dynamic.PAM9_a

licor.NPQ

Cx = 0.3187.*NPQ
    Cx(Cx>1.5) = 1.5;
    leafbio.fqe=[phi2F];
    leafbio.Cx = Cx;
%%

leafoptics = main_flu(spectral, leafbio,optipar,"not_combined",E,false);


par = obs.aus_data.licor_co2.Qin(1,i)

% define the spectral regions
wlp      = spectral.wlP;         % PROSPECT wavelengths as a column-vector
wlPAR = spectral.wlPAR';
IparP    = find(wlp>=400 & wlp<=700); % Indices for PAR wavelenghts within wl
wl_low = 400;
wl_high = 2000;
wl_obs = spectral.wl_ASD;
wlp = spectral.wlP;
I_obs = find(wl_obs>=wl_low & wl_obs<=wl_high); % data after 2000 nm is noisy
I_sim = find(wlp>=wl_low & wlp<=wl_high); % data after 2000 nm is noisy

% calculate absorptance
refl = leafopt.refl;
tran = leafopt.tran;
absorb = 1-refl-tran;


[out_opt out_photo] = main_fun_v1(spectral,leafbio,optipar,E,v)




%% William's Woodgate data: LRC experiment

% Prepare observations
licor_lrc = obs.aus_data.licor_lrc;

% Prepare photosynthesis data
n_sample = length(obs.aus_data.names);  % Number of leaves
n_exp = length(licor_lrc.Qin) % Number of light levels

i = 1;
% data.Qin = transpose(linspace(1,2400,n));   % PAR, umol PPFD m-2 s-1
data.Qin = licor_lrc.Qin(:,i);
data.Tin = licor_lrc.Tleaf(:,i);
% data.Tin = repmat(25,n_exp,1);                  % Leaf temperature, C
% data.Cin = repmat(200,n_exp,1);                 % Mesophyll CO2, ubar
data.Cin = licor_lrc.Ci(:,i);                 % Mesophyll CO2, ubar
data.Oin = repmat(209,n_exp,1);                 % Atmospheric O2, mbar
v = configure_fun(data);

m = model_fun(v);
plot(m.An_a*1e6)


% Configure model simulations

% Load default set of parameters

% Adjust parameters of interest from default values
v.Abs = 0.85;           % Total leaf absorptance to PAR, mol mol-1
v.beta = 0.52;          % PSII fraction of total absorptance, mol mol-1
v.Ku2 = 0e09;           % Rate constant for exciton sharing at PSII, s-1
v.CB6F = 1.2./1e6;      % Cyt b6f density, mol sites m-2
v.RUB = 27.8./1e6;      % Rubisco density, mol sites m-2
v.eps1 = 0;             % PS I transfer function, mol mol-1
m = model_fun(v);
plot(m.An_a)

plot(m.An_a*1e6,'k')
hold on
plot(licor_lrc.A(:,i),'r')
xlabel('PAR')
ylabel('An')
legend(["Simulation", "Observation"])


env.Qin = licor_lrc.Qin(:,i);        % PAR, umol PPFD m-2 s-1
env.Tin = repmat(25,n_exp,1);                                     % Leaf temperature, C
env.Cin = obs.aus_data.licor_lrc.Ci(:,i);                                     % Mesophyll CO2, ubar
% env.Cin = repmat(200,n_exp,1);                                    % Atmospheric O2, mbar
env.Oin = repmat(209,n_exp,1);                                    % Atmospheric O2, mbar
v = configure_fun(env);


% DYNAMIC 
symsolver_fun();
% Adjust parameters of interest from default values
v.Abs = 0.85;           % Total leaf absorptance to PAR, mol mol-1
v.beta = 0.52;          % PSII fraction of total absorptance, mol mol-1
v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
v.CB6F = 1.2./1e6;      % Cyt b6f density, mol sites m-2
v.RUB = 27.8./1e6;      % Rubisco density, mol sites m-2
v.eps1 = 0;             % PS I transfer function, mol mol-1

% Set dynamic cross-sections and assign optimization function to 'v' 
v.alpha_opt = 'dynamic';
v.solve_xcs = solve_xcs;
m = model_fun(v);
plot(licor_lrc.Qin, m.An_a.*1e6,'k')
hold on
plot(licor_lrc.Qin,licor_lrc.A(:,i),'r')
xlabel('PAR')
ylabel('An')
legend(["Simulation", "Observation"])

%%
load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);



