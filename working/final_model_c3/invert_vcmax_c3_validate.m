clear
clc
warning('off')
data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
script_dir = pwd;
symsolver_fun();
[obs, spectral] = prepare_inputs(data_dir);

load est_vcmax_vqmax_stat_co2.mat
load est_vcmax_vqmax_dynamic_co2.mat

licor = obs.aus_data.licor_lrc;
n_sample = length(obs.aus_data.names);  % Number of leaves
n = length(licor.Qin); % Number of light levels
%% Run the LRC with the estimated Vcmax
for i = 1:8

    an_obs = licor.A(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i)
    data.Cin = licor.Ci(:,i);
%     data.Cin = repmat(200,n,1);   
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;

    m_vcmax_static = model_fun_c3_vcmax_vqmax(est_vcmax_vqmax_stat(1,i)...
        ,est_vcmax_vqmax_stat(2,i),v);
   
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    m_vcmax_dynamic = model_fun_c3_vcmax_vqmax(est_vcmax_vqmax_stat(1,i)...
        ,est_vcmax_vqmax_stat(2,i),v);

    subplot(3,3,i)
    plot(data.Qin, an_obs,'k')
    hold on
    plot(data.Qin, m_vcmax_static.An_a*1e6,'b')
    hold on
    plot(data.Qin, m_vcmax_dynamic.An_a*1e6,'r')
    xlabel("PAR")
    ylabel("An")

end

