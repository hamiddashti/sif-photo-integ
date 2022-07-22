clear
clc
warning('off')
data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
script_dir = pwd;
% [obs, spectral] = prepare_inputs(data_dir);
load obs.mat
% William's Woodgate data: CO2 experiment

% Prepare observations
licor = obs.aus_data.licor_co2;
n_sample = length(obs.aus_data.names);  % Number of leaves
n = length(licor.Ci); % Number of Co2 levels
pathway_opt='C3';
%% Run the inversion

vcmax_0 =  50;
vars_0 = vcmax_0;
opts.LBounds = 10; opts.UBounds = 100;
opts.Restarts=3;
opts.Noise.on=0;

for j = 1:8
    an_obs = licor.A(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    data.Cin = licor.Ci(:,j);                 % Mesophyll CO2, ubar
    data.Oin = repmat(0.209,n,1);               % Atmospheric O2, bar
    data.Pin = repmat(1,n,1);

    % Set all parameters
    v = configure_fun(data);
    if strcmp(pathway_opt,'C3') == 1
        % Case 1: Pure C3
        v.Ku2 = 2e09;
        v.CB6F = 175./v.kq.*1e-06;
        v.RUB = 50./v.kc.*1e-06;
        v.abs_frac = 0.0;       %Bundle sheath parameter
        v.vq_frac = 0.0;        %Bundle sheath parameter
        v.vc_frac = 0.0;        %Bundle sheath parameter
        v.a2_s_frac = 0.0;      %Bundle sheath parameter
        v.a2_m_frac = 0.52;
    end
    est_vcmax(j) = cmaes('chi2_vcmax',vars_0,[],opts,v,an_obs);
end
disp("Inversion done!")

%% Plot the vcmax results
fig = figure (1);
set(gcf,'units','inch','Position',[50 50 15 12],'color','w');

for j=1:8
    an_obs = licor.A(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    data.Cin = licor.Ci(:,j);                 % Mesophyll CO2, ubar
    data.Oin = repmat(0.209,n,1);               % Atmospheric O2, bar
    data.Pin = repmat(1,n,1);

    % Set all parameters
    v = configure_fun(data);
    if strcmp(pathway_opt,'C3') == 1
        % Case 1: Pure C3
        v.Ku2 = 2e09;
        v.CB6F = 175./v.kq.*1e-06;
        v.RUB = 50./v.kc.*1e-06;
        v.abs_frac = 0.0;       %Bundle sheath parameter
        v.vq_frac = 0.0;        %Bundle sheath parameter
        v.vc_frac = 0.0;        %Bundle sheath parameter
        v.a2_s_frac = 0.0;      %Bundle sheath parameter
        v.a2_m_frac = 0.52;
    end

    m_vcmax = model_fun_c3c4_vcmax(est_vcmax_stat(j),v);
    
    subplot(3,3,j)
    plot(licor.Ci(:,j), m_vcmax.An_a*1e6,'k')
    hold on
    plot(licor.Ci(:,j),an_obs,'r')
    ylim([-5,35])
    xlabel("Ci",fontsize=12,FontWeight="bold")
    ylabel("An",fontsize=12,FontWeight="bold")
    text(700,10,strcat("vcmax: ",num2str(round(est_vcmax(j),1))),fontsize=10)
    legend(["static","obs"])
end
sgtitle(gcf,"Only Vcmax")
saveas(gcf,"./Figures/co2_exp_vcmax.png")


%% Vcmax and Vqmax

vcmax_0 =  50;
vqmax_0 = 175;
vars_0 = [vcmax_0, vqmax_0]';
opts.LBounds = [10, 50]'; opts.UBounds = [100, 400]';
opts.Restarts=3;
opts.Noise.on=0;

for j = 1:8
    an_obs = licor.A(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    data.Cin = licor.Ci(:,j);                 % Mesophyll CO2, ubar
    data.Oin = repmat(0.209,n,1);               % Atmospheric O2, bar
    data.Pin = repmat(1,n,1);

    % Set all parameters
    v = configure_fun(data);
    if strcmp(pathway_opt,'C3') == 1
        % Case 1: Pure C3
        v.Ku2 = 2e09;
        v.CB6F = 175./v.kq.*1e-06;
        v.RUB = 50./v.kc.*1e-06;
        v.abs_frac = 0.0;       %Bundle sheath parameter
        v.vq_frac = 0.0;        %Bundle sheath parameter
        v.vc_frac = 0.0;        %Bundle sheath parameter
        v.a2_s_frac = 0.0;      %Bundle sheath parameter
        v.a2_m_frac = 0.52;
    end
    est_vcmax_vqmax(:,j) = cmaes('chi2_vcmax_vqmax',vars_0,[],opts,v,an_obs);
end
disp("Inversion done!")

%% Plot the vcmax results
close all
fig = figure (1);
set(gcf,'units','inch','Position',[50 50 15 12],'color','w');

for j=1:8
    an_obs = licor.A(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    data.Cin = licor.Ci(:,j);                 % Mesophyll CO2, ubar
    data.Oin = repmat(0.209,n,1);               % Atmospheric O2, bar
    data.Pin = repmat(1,n,1);

    % Set all parameters
    v = configure_fun(data);
    if strcmp(pathway_opt,'C3') == 1
        % Case 1: Pure C3
        v.Ku2 = 2e09;
        v.CB6F = 175./v.kq.*1e-06;
        v.RUB = 50./v.kc.*1e-06;
        v.abs_frac = 0.0;       %Bundle sheath parameter
        v.vq_frac = 0.0;        %Bundle sheath parameter
        v.vc_frac = 0.0;        %Bundle sheath parameter
        v.a2_s_frac = 0.0;      %Bundle sheath parameter
        v.a2_m_frac = 0.52;
    end

    m_vcmax_vqmax = model_fun_c3c4_vcmax_vqmax(est_vcmax_vqmax(1,j),...
        est_vcmax_vqmax(2,j),v);
    
    subplot(3,3,j)
    plot(licor.Ci(:,j), m_vcmax_vqmax.An_a*1e6,'k')
    hold on
    plot(licor.Ci(:,j),an_obs,'r')
    ylim([-5,35])
    xlabel("Ci",fontsize=12,FontWeight="bold")
    ylabel("An",fontsize=12,FontWeight="bold")
    text(700,10,strcat("vcmax: ",num2str(round(est_vcmax_vqmax(1,j),1))),fontsize=10)
    text(700,7,strcat("vqmax: ",num2str(round(est_vcmax_vqmax(2,j),1))),fontsize=10)
    legend(["static","obs"])
end
sgtitle(gcf,"Vcmax and Vqmax")
saveas(gcf,"./Figures/co2_exp_vcmax_vqmax.png")
%% Vcmax, Vqmax and EPS

vcmax_0 =  50;
vqmax_0 = 175;
Eps1 = 1;
Eps2=1;
vars_0 = [vcmax_0, vqmax_0, Eps1,Eps2]';
opts.LBounds = [10, 50,0,0]'; opts.UBounds = [100, 400,1,1]';
opts.Restarts=3;
opts.Noise.on=0;
wA = 1;
wP = 1;
wNPQ = 1;

for j = 2:2
    an_obs = licor.A(:,j);
%     PhiPS2_obs = licor.PhiPS2(:,j);
    NPQ_obs = licor.NPQ(:,j);
    fs_obs = licor.Fs(:,j)./licor.Fo(:,j);
    data.Qin = licor.Qin(:,j);
    data.Tin = licor.Tleaf(:,j);
    data.Cin = licor.Ci(:,j);                 % Mesophyll CO2, ubar
    data.Oin = repmat(0.209,n,1);               % Atmospheric O2, bar
    data.Pin = repmat(1,n,1);

    % Set all parameters
    v = configure_fun(data);
    if strcmp(pathway_opt,'C3') == 1
        % Case 1: Pure C3
        v.Ku2 = 2e09;
        v.CB6F = 175./v.kq.*1e-06;
        v.RUB = 50./v.kc.*1e-06;
        v.abs_frac = 0.0;       %Bundle sheath parameter
        v.vq_frac = 0.0;        %Bundle sheath parameter
        v.vc_frac = 0.0;        %Bundle sheath parameter
        v.a2_s_frac = 0.0;      %Bundle sheath parameter
        v.a2_m_frac = 0.52;
    end

    est_vcmax_vqmax_eps(:,j) = cmaes('chi2_vcmax_vqmax_eps',vars_0,[],...
        opts,v,an_obs,fs_obs,NPQ_obs,wA,wP,wNPQ);
    m1 = model_fun_c3c4_vcmax_vqmax_eps(est_vcmax_vqmax_eps(1,j),...
    est_vcmax_vqmax_eps(2,j),est_vcmax_vqmax_eps(3,j),est_vcmax_vqmax_eps(4,j),v);
    
    figure(1)
    plot(m1.An_a.*1e6)
    hold on
    plot(an_obs,'r')
    
    plot(m1.phi2F_ma)
    
    figure(4)
    plot(fs_obs,'r')
    hold on
    plot(m1.Fs_a./m1.Fo_a)
    
    figure(3)
    plot(NPQ_obs,'r')
    hold on
    plot(m1.PAM3_a,'b')

 end
disp("Inversion done!")

%%
est_vcmax_vqmax_eps(:,j)
m1 = model_fun_c3c4_vcmax_vqmax_eps(est_vcmax_vqmax_eps(1,j),...
    est_vcmax_vqmax_eps(2,j),est_vcmax_vqmax_eps(3,j),est_vcmax_vqmax_eps(4,j),v);

figure(1)
plot(m1.An_a.*1e6)
hold on
plot(an_obs,'r')

figure(2)
plot(m1.PAM4_a)
hold on
plot(licor.PhiPS2(:,j))

figure(3)
plot(licor.NPQ(:,j))
hold on
plot(m1.PAM3_a,'r')

figure(4)
plot(fs_obs)
hold on
plot(m1.Fs_a./m1.Fo_a)




