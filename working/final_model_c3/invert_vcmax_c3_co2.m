clear
clc
warning('off')
data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
script_dir = pwd;
symsolver_fun();

% [obs, spectral] = prepare_inputs(data_dir);
load obs.mat 
load spectral.mat
% William's Woodgate data: CO2 experiment

% Prepare observations
licor = obs.aus_data.licor_co2;
n_sample = length(obs.aus_data.names);  % Number of leaves
n = length(licor.Ci); % Number of Co2 levels
load optipar.mat
% load est_chems_co2.mat   % Estimated chemistry from invert_chem_v1.m
load chems_inversion.mat


%% Inversion only vcmax
vcmax_0 =  50;
vars_0 = vcmax_0;
opts.LBounds = 10; opts.UBounds = 100;
opts.Restarts=3;
opts.Noise.on=1;
for i = 1:1
    
    % Prepare model params
    an_obs = licor.A(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
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
    v.Abs = chems_inversion.Abs(i);
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
    data.Cin = licor.Ci(:,j);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05
    v.Abs = chems_inversion.Abs(j);
    m_vcmax_static = model_fun_c3_vcmax(est_vcmax_stat(j),v);

    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    v.Abs = chems_inversion.Abs(j);
    m_vcmax_dynamic = model_fun_c3_vcmax(est_vcmax_dynamic(j),v);

    hold on
    subplot(3,3,j)
    plot(licor.Ci(:,j), m_vcmax_static.An_a*1e6,'k')
    hold on
    plot(licor.Ci(:,j), m_vcmax_dynamic.An_a*1e6,'b')
    hold on
    plot(licor.Ci(:,j),an_obs,'r')
    ylim([-5,35])
    xlabel("Ci",fontsize=12,FontWeight="bold")
    ylabel("An",fontsize=12,FontWeight="bold")
    text(700,10,strcat("vcmax-static: ",num2str(round(est_vcmax_stat(j),1))),fontsize=10)
    text(700,7,strcat("vcmax-dynamic: ",num2str(round(est_vcmax_dynamic(j),1))),fontsize=10)
    legend(["static","dynamic","obs"])
end
sgtitle(gcf,"Only Vcmax")

saveas(gcf,"./Figures/co2_exp_vcmax.png")

%% Vcmax and Vqmax

vcmax_0 =  50;
vqmax_0 = 175;
vars_0 = [vcmax_0, vqmax_0]';
opts.LBounds = [10, 50]'; opts.UBounds = [150, 400]';
opts.Restarts=3;
opts.Noise.on=0;

for i = 1:8

    % Prepare model params
    an_obs = licor.A(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
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
    v.Abs = chems_inversion.Abs(i);
    est_vcmax_vqmax_dynamic(:,i) = cmaes('chi2_vcmax_vqmax',vars_0,[],opts,v,an_obs);

end
disp("Inversion done!")
save est_vcmax_vqmax_stat_co2.mat est_vcmax_vqmax_stat
save est_vcmax_vqmax_dynamic_co2.mat est_vcmax_vqmax_dynamic

%% Plot the vcmax vqmax results
load est_vcmax_vqmax_stat_co2.mat
load est_vcmax_vqmax_dynamic_co2.mat

fig = figure (2);
set(gcf,'units','inch','Position',[50 50 15 12],'color','w');

for j=1:8

    an_obs = licor.A(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    data.Cin = licor.Ci(:,j);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.Ku2 = 2e09;           % Rate constant for exciton sharing at PSII, s-1
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(j);
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
    v.Abs = chems_inversion.Abs(j);
    m_vcmax_dynamic = model_fun_c3_vcmax_vqmax(est_vcmax_vqmax_dynamic(1,j)...
        ,est_vcmax_vqmax_dynamic(2,j),v);

    hold on
    subplot(3,3,j)
    plot(licor.Ci(:,j), m_vcmax_static.An_a*1e6,'k')
    hold on
    plot(licor.Ci(:,j), m_vcmax_dynamic.An_a*1e6,'b')
    hold on
    plot(licor.Ci(:,j),an_obs,'r')
    ylim([-5,35])
    xlabel("Ci",fontsize=12,fontweight="bold")
    ylabel("An",fontsize=12,fontweight="bold")
    text(250,28,strcat("vcmax-static: ",num2str(round(est_vcmax_vqmax_stat(1,j),1))),fontsize=10)
    text(250,25,strcat("vqmax-static: ",num2str(round(est_vcmax_vqmax_stat(2,j),1))),fontsize=10)
    text(700,7,strcat("vcmax-dynamic: ",num2str(round(est_vcmax_vqmax_dynamic(1,j),1))),fontsize=10)
    text(700,4,strcat("vqmax-dynamic: ",num2str(round(est_vcmax_vqmax_dynamic(2,j),1))),fontsize=10)

    legend(["static","dynamic","obs"])
end
sgtitle("Optimizing Vcmax and Vqmax")
saveas(gcf,"./Figures/co2_exp_vcmax_vqmax.png")




%% Vcmax, Vqmax ku2

vcmax_0 =  50;
vqmax_0 = 175;
ku2_0 = 2;

vars_0 = [vcmax_0, vqmax_0,ku2_0]';
opts.LBounds = [10, 50,0.0]'; opts.UBounds = [100, 400,4]';
opts.Restarts=3;
opts.Noise.on=0;

for i = 1:8
    disp(i)
    % Prepare model params
    an_obs = licor.A(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
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
    v.Abs = chems_inversion.Abs(i);
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
    data.Cin = licor.Ci(:,j);
    data.Oin = repmat(209,n,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05
    v.Abs = chems_inversion.Abs(j);
    m_vcmax_static = model_fun_c3_vcmax_vqmax_ku2(est_vcmax_vqmax_ku2_stat(1,j),...
        est_vcmax_vqmax_ku2_stat(2,j),est_vcmax_vqmax_ku2_stat(3,j),v);

    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;
    v.Rds = 0.05
    % Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    v.Abs = chems_inversion.Abs(j);
    m_vcmax_dynamic = model_fun_c3_vcmax_vqmax_ku2(est_vcmax_vqmax_ku2_dynamic(1,j),...
        est_vcmax_vqmax_ku2_dynamic(2,j),est_vcmax_vqmax_ku2_dynamic(3,j),v);

    hold on
    subplot(3,3,j)
    plot(licor.Ci(:,j), m_vcmax_static.An_a*1e6,'k')
    hold on
    plot(licor.Ci(:,j), m_vcmax_dynamic.An_a*1e6,'b')
    hold on
    plot(licor.Ci(:,j),an_obs,'r')
    ylim([-5,35])
    xlabel("Ci",fontsize=12,fontweight="bold")
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
saveas(gcf,"./Figures/co2_exp_vcmax_vqmax_k2.png")


%% plot vs vcmax
fig = figure(1);
set(gcf, 'units','inch','Position',[5,5,5,5]);
vcmax = est_vcmax_vqmax_stat(1,:);
chl = chems_inversion.est_chems(1,:);
mdl = fitlm(chl,vcmax);
plotAdded(mdl,[],Marker='+',Markersize=8)
xlabel('Chl [\mug cm^2] ','interpreter','Tex')
ylabel('Vcmax [mol CO^2 m^{-2} s^{-1}]','interpreter','Tex')
eq1 = "Vcmax = " + mdl.Coefficients.Estimate(1) + "+" + ...
    mdl.Coefficients.Estimate(2) + "Chl"
eq2 = "R^2 = " + mdl.Rsquared.Ordinary + ...
    " (P_{value} = " + round(coefTest(mdl),3) +")"
text(35,60, eq1,FontWeight='Bold')
text(35,55, eq2,FontWeight='Bold')
% text(35,50, eq3,FontWeight='Bold')
legend(gca,'off')
saveas(gcf,"./Figures/chl_vcmax.png")

