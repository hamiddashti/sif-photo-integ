clear
clc
warning('off')
data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
script_dir = pwd;
symsolver_fun();
[obs, spectral] = prepare_inputs(data_dir);

% William's Woodgate data: CO2 experiment

% Prepare observations
licor = obs.aus_data.licor_lrc;
n_sample = length(obs.aus_data.names);  % Number of leaves
n = length(licor.Qin); % Number of light levels

%% Inversion only vcmax
vcmax_0 =  50;
vars_0 = vcmax_0;
opts.LBounds = 10; opts.UBounds = 100;
opts.Restarts=3;
opts.Noise.on=1;
for i = 1:8
    disp(i)
    % Prepare model params
    an_obs = licor.A(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    %     data.Cin = licor.Ci(:,i);
    data.Cin = repmat(250,n,1);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;

    % Invert only for vcmax
    % Static

    est_vcmax_stat(i) = cmaes('chi2_vcmax',vars_0,[],opts,v,an_obs);

    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    % [vars_est,fmin,counteval,stopflag,out,bestever] = cmaes('chi2P5B',vcmax_0,[],opts,v,an_obs);
    est_vcmax_dynamic(i) = cmaes('chi2_vcmax',vars_0,[],opts,v,an_obs);

end
disp("Inversion done!")

%% Plot the vcmax results
fig = figure (1);
set(gcf,'units','inch','Position',[50 50 15 12],'color','w');

for j=1:8

    an_obs = licor.A(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    %     data.Cin = licor.Ci(:,j);
    data.Cin = repmat(250,n,1);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05
    m_vcmax_static = model_fun_c3_vcmax(est_vcmax_stat(j),v);

    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    m_vcmax_dynamic = model_fun_c3_vcmax(est_vcmax_dynamic(j),v);

    hold on
    subplot(3,3,j)
    plot(licor.Qin(:,j), m_vcmax_static.An_a*1e6,'k')
    hold on
    plot(licor.Qin(:,j), m_vcmax_dynamic.An_a*1e6,'b')
    hold on
    plot(licor.Qin(:,j),an_obs,'r')
    ylim([-5,35])
    xlabel("PAR",fontsize=12,FontWeight="bold")
    ylabel("An",fontsize=12,FontWeight="bold")
    text(200,28,strcat("vcmax-static: ",num2str(round(est_vcmax_stat(j),1))),fontsize=10)
    text(200,25,strcat("vcmax-dynamic: ",num2str(round(est_vcmax_dynamic(j),1))),fontsize=10)
    legend(["static","dynamic","obs"])
end
sgtitle(gcf,"Only Vcmax")

saveas(gcf,"./Figures/lrc_exp_vcmax.png")

%% Vcmax and Vqmax

vcmax_0 =  50;
vqmax_0 = 175;
vars_0 = [vcmax_0, vqmax_0]';
opts.LBounds = [10, 50]'; opts.UBounds = [100, 400]';
opts.Restarts=3;
opts.Noise.on=0;

for i = 1:8
    disp(i)
    % Prepare model params
    an_obs = licor.A(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    %     data.Cin = licor.Ci(:,i);
    data.Cin = repmat(250,n,1);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;

    % Static
    [est_vcmax_vqmax_stat(:,i),fmin] = cmaes('chi2_vcmax_vqmax',vars_0,[],opts,v,an_obs);


    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    % [vars_est,fmin,counteval,stopflag,out,bestever] = cmaes('chi2P5B',vcmax_0,[],opts,v,an_obs);
    est_vcmax_vqmax_dynamic(:,i) = cmaes('chi2_vcmax_vqmax',vars_0,[],opts,v,an_obs);

end
disp("Inversion done!")
save est_vcmax_vqmax_stat_lrc.mat est_vcmax_vqmax_stat
save est_vcmax_vqmax_dynamic_lrc.mat est_vcmax_vqmax_dynamic


%% Plot the vcmax vqmax results
fig = figure (2);
set(gcf,'units','inch','Position',[50 50 15 12],'color','w');

for j=1:8

    an_obs = licor.A(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    %     data.Cin = licor.Ci(:,j);
    data.Cin = repmat(250,n,1);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    m_vcmax_static = model_fun_c3_vcmax_vqmax(est_vcmax_vqmax_stat(1,j)...
        ,est_vcmax_vqmax_stat(2,j),v);

    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05;
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    m_vcmax_dynamic = model_fun_c3_vcmax_vqmax(est_vcmax_vqmax_dynamic(1,j)...
        ,est_vcmax_vqmax_dynamic(2,j),v);

    hold on
    subplot(3,3,j)
    plot(licor.Qin(:,j), m_vcmax_static.An_a*1e6,'k')
    hold on
    plot(licor.Qin(:,j), m_vcmax_dynamic.An_a*1e6,'b')
    hold on
    plot(licor.Qin(:,j),an_obs,'r')
    ylim([-5,35])
    xlabel("PAR",fontsize=12,fontweight="bold")
    ylabel("An",fontsize=12,fontweight="bold")
    text(250,28,strcat("vcmax-static: ",num2str(round(est_vcmax_vqmax_stat(1,j),1))),fontsize=10)
    text(250,25,strcat("vqmax-static: ",num2str(round(est_vcmax_vqmax_stat(2,j),1))),fontsize=10)
    text(700,7,strcat("vcmax-dynamic: ",num2str(round(est_vcmax_vqmax_dynamic(1,j),1))),fontsize=10)
    text(700,4,strcat("vqmax-dynamic: ",num2str(round(est_vcmax_vqmax_dynamic(2,j),1))),fontsize=10)

    legend(["static","dynamic","obs"])
end
sgtitle("Optimizing Vcmax and Vqmax")
saveas(gcf,"./Figures/lrc_exp_vcmax_vqmax.png")

%% Vcmax, Vqmax ku2

vcmax_0 =  50;
vqmax_0 = 175;
ku2_0 = 2;

vars_0 = [vcmax_0, vqmax_0,ku2_0]';
opts.LBounds = [10, 50,0.0]'; opts.UBounds = [60, 400,4]';
opts.Restarts=3;
opts.Noise.on=0;

for i = 1:8
    disp(i)
    % Prepare model params
    an_obs = licor.A(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = repmat(250,n,1);
    %     data.Cin = licor.Ci(:,i);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    % Static
    [est_vcmax_vqmax_ku2_stat(:,i),fmin] = cmaes('chi2_vcmax_vqmax_ku2',vars_0,[],opts,v,an_obs);

    v = configure_fun_c3(data);
    %     v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05;
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    % [vars_est,fmin,counteval,stopflag,out,bestever] = cmaes('chi2P5B',vcmax_0,[],opts,v,an_obs);
    est_vcmax_vqmax_ku2_dynamic(:,i) = cmaes('chi2_vcmax_vqmax_ku2',vars_0,[],opts,v,an_obs);

end
disp("Inversion done!")

%% Plot the vcmax vqmax ku2 results
fig = figure (3);
set(gcf,'units','inch','Position',[50 50 15 12],'color','w');

for j=1:8

    an_obs = licor.A(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    %     data.Cin = licor.Ci(:,j);
    data.Cin = repmat(250,n,1);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05
    m_vcmax_static = model_fun_c3_vcmax_vqmax_ku2(est_vcmax_vqmax_ku2_stat(1,j),...
        est_vcmax_vqmax_ku2_stat(2,j),est_vcmax_vqmax_ku2_stat(3,j),v);

    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    m_vcmax_dynamic = model_fun_c3_vcmax_vqmax_ku2(est_vcmax_vqmax_ku2_dynamic(1,j),...
        est_vcmax_vqmax_ku2_dynamic(2,j),est_vcmax_vqmax_ku2_dynamic(3,j),v);

    hold on
    subplot(3,3,j)
    plot(licor.Qin(:,j), m_vcmax_static.An_a*1e6,'k')
    hold on
    plot(licor.Qin(:,j), m_vcmax_dynamic.An_a*1e6,'b')
    hold on
    plot(licor.Qin(:,j),an_obs,'r')
    ylim([-5,35])
    xlabel("PAR",fontsize=12,fontweight="bold")
    ylabel("An",fontsize=12,fontweight="bold")
    text(200,28,strcat("vcmax-static: ",num2str(round(est_vcmax_vqmax_ku2_stat(1,j),1))),fontsize=10)
    text(200,25,strcat("vq-static: ",num2str(round(est_vcmax_vqmax_ku2_stat(2,j),1))),fontsize=10)
    text(200,22,strcat("ku2-static: ",num2str(round(est_vcmax_vqmax_ku2_stat(3,j),1))),fontsize=10)

    text(700,15,strcat("vcmax-dynamic: ",num2str(round(est_vcmax_vqmax_ku2_dynamic(1,j),1))),fontsize=10)
    text(700,12,strcat("vqmax-dynamic: ",num2str(round(est_vcmax_vqmax_ku2_dynamic(2,j),1))),fontsize=10)
    text(700,9,strcat("ku2-dynamic: ",num2str(round(est_vcmax_vqmax_ku2_dynamic(3,j),1))),fontsize=10)


    legend(["static","dynamic","obs"])
end
sgtitle("Optimizing for Vcmax, Vqmax and Ku2")
saveas(gcf,"./Figures/lrc_exp_vcmax_vqmax_k2.png")


