clear
clc
warning('off')
data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
script_dir = pwd;
% Prepare observations
load obs.mat % Load observation
load spectral.mat
load optipar.mat
load chems_inversion.mat

licor = obs.aus_data.licor_lrc;
n_steps = size(licor.Qin,1);
n_samples = size(licor.Qin,2);
p = 970;
ppm2bar     =  1e-6 .* (p .*1E-3);
symsolver_fun();
names = obs.aus_data.names;

%% We do the inversion on light experiement

vcmax_0 =  50;
% vqmax_0 = 175;
ku2_0 = 2;
% a1_0 = 0.5;
% a2_0 = 0.5;
% eps1_0 = 0.5;
% eps2_0 = 0.25;
% vars_0 = [vcmax_0, vqmax_0,eps1_0,eps2_0]';
vars_0 = [vcmax_0, ku2_0]';
% opts.LBounds = [10, 50,0,0]'; opts.UBounds = [150, 500,1,1]';
opts.LBounds = [1,1]'; opts.UBounds = [150,4]';
opts.Restarts=3;
opts.Noise.on=1;
weights.wAn = 1;
weights.wNPQ = 0;
weights.wFs = 0;

for i = 1:n_samples

    observed.an_obs = licor.A(:,i);
    observed.npq_obs = licor.NPQ(:,i);
    observed.fs_obs = licor.Fs(:,i)./licor.Fo(:,i);

    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps,1);                 % Atmospheric O2, mbar

    % Static
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
%     v.Ku2 = 2e09;
    est_pars_static(:,i) = cmaes('chi2_photo_lrc',vars_0,[],opts,v,observed,weights);

% %     Dynamic
%     v.alpha_opt = 'dynamic';
%     v.solve_xcs = solve_xcs;
%     est_pars_dynamic(:,i) = cmaes('chi2_photo_lrc',vars_0,[],opts,v,observed,weights);

end
% save est_pars_static_vcvq.mat est_pars_static
% save est_pars_dynamic_vcvq.mat est_pars_dynamic
disp("Inversion done!")

%% Plot
close all
Fsize = 10;
fig1 = figure(1);
set(fig1,'units','inch','Position',[5 5 10 8],'color','w');

fig2 = figure(2);
set(fig2,'units','inch','Position',[10 10 10 8],'color','w');

fig3 = figure(3);
set(fig3,'units','inch','Position',[10 10 10 8],'color','w');

for i=1:n_samples
    an_obs = licor.A(:,i);
    npq_obs = licor.NPQ(:,i);
    fs_obs = licor.Fs(:,i)./licor.Fo(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
    v.Ku2 = 2e09;
    free_pars.vcmax = est_pars_static(1,i);
    free_pars.vqmax = est_pars_static(2,i);
    %     free_pars.ku2 = est_pars_static(3,i);
    %     free_pars.eps1 = est_pars_static(1,i);
    %     free_pars.eps2 = est_pars_static(2,i);
    sim_m_static=model_fun_c3_inv_lrc(v,free_pars);
    %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    free_pars.vcmax = est_pars_dynamic(1,i);
    free_pars.vqmax = est_pars_dynamic(2,i);
    %     free_pars.ku2 = est_pars_dynamic(3,i);
    %     free_pars.eps1 = est_pars_dynamic(1,i);
    %     free_pars.eps2 = est_pars_dynamic(2,i);
    sim_m_dynamic=model_fun_c3_inv_lrc(v,free_pars);

    set(0,'CurrentFigure',fig1)
    hold on
    subplot(3,3,i)
    plot(licor.Qin(:,i), sim_m_static.An_a*1e6,'k',LineWidth=1)
    hold on
    plot(licor.Qin(:,i), sim_m_dynamic.An_a*1e6,'b',LineWidth=1)
    hold on
    plot(licor.Qin(:,i),an_obs,'r',LineWidth=1)
    ylim([-5,35])
    xlabel("PAR",fontsize=Fsize)
    title(names(i)) 

    if (i==1)|(i==4)|(i==7)
        ylabel("A_{n} [µmol m^{-2} s^{-1}]",fontsize=Fsize)
    end
    text(250,30,strcat("vcmax_s: ",num2str(round(est_pars_static(1,i),1)))...
        ,fontsize=Fsize)
    text(250,26,strcat("vqmax_s: ",num2str(round(est_pars_static(2,i),1))),...
        fontsize=Fsize)
    text(700,5,strcat("vcmax_d: ",num2str(round(est_pars_dynamic(1,i),1))),...
        fontsize=Fsize)
    text(700,1,strcat("vqmax_d: ",num2str(round(est_pars_dynamic(2,i),1))),...
        fontsize=Fsize)
    set(gca,"FontSize",12)
    sgtitle("An")


    set(0,'CurrentFigure',fig2)
    hold on
    subplot(3,3,i)
    plot(licor.Qin(:,i), sim_m_static.Kn2_a.*1e-9,'k',LineWidth=1)
    %     plot(licor.Ci(:,i), sim_m_static.PAM9_a,'k',LineWidth=1)
    hold on
    plot(licor.Qin(:,i), sim_m_dynamic.Kn2_a.*1e-9,'b',LineWidth=1)
    %     plot(licor.Ci(:,i), sim_m_dynamic.PAM9_a,'b',LineWidth=1)
    hold on
    plot(licor.Qin(:,i),npq_obs,'r',LineWidth=1)
    sgtitle("NPQ")

    set(0,'CurrentFigure',fig3)
    hold on
    subplot(3,3,i)
    plot(licor.Qin(:,i), sim_m_static.Fs_a./sim_m_static.Fo_a,'k',LineWidth=1)
    hold on
    plot(licor.Qin(:,i), sim_m_dynamic.Fs_a./sim_m_dynamic.Fo_a,'b',LineWidth=1)
    hold on
    plot(licor.Qin(:,i),fs_obs,'r',LineWidth=1)
    sgtitle("Fs")
end
% Lgnd = legend(["static","dynamic","obs"])
% Lgnd.Position(1) = 1.1;
% Lgnd.Position(2) = 0.2;
% Lgnd.FontSize=14;

set(0,'CurrentFigure',fig1)
newPosition = [0.5 0.03 0.01 0.01];
newUnits = 'normalized';
hL = legend(["static","dynamic","obs"],"orientation","horizontal"...
    ,"location","southoutside")
hL = legend(["static","obs"],"orientation","horizontal","location"...
    ,"southoutside")
set(hL,'Position', newPosition,'Units', newUnits,"FontSize",12);

% saveas(gcf,'./Figures/fit_Vcmax_Vqmax.png')


%% Plot vcmax
fig4 = figure(4);
Fsize=12;
set(fig4,'units','inch','Position',[5 5 4 4],'color','w');
% vcmax = est_pars_dynamic(1,:);
vcmax = est_pars_static(1,:);
chl = chems_inversion.est_chems(1,:);
mdl = fitlm(chl,vcmax);
scatter(chl,vcmax,Marker="+",LineWidth=2,MarkerEdgeColor ="k")
h1 = lsline;
h1.LineStyle = "--"
h1.LineWidth = 2;
box on
% plotAdded(mdl,[],Markersize=12,color="k",LineWidth=3)

xlabel('Chl [\mug cm^2] ','interpreter','Tex')
ylabel('Vcmax [mol CO^2 m^{-2} s^{-1}]','interpreter','Tex')

eq1 = "Vcmax = " + round(mdl.Coefficients.Estimate(1),2) + "+" + ...
    round(mdl.Coefficients.Estimate(2),2) + "Chl";
eq2 = "R^2 = " + round(mdl.Rsquared.Ordinary,2);
text(40,65, eq1,FontSize=Fsize)
text(40,58, eq2,FontSize=Fsize)
% text(35,50, eq3,FontWeight='Bold')
legend(gca,'off')


%% Lets run the model wiht default values
par_idx = 2:n_steps(end); % skipping the PAR~0
for i =1:8

    observed.an_obs = licor.A(par_idx,i);
    observed.npq_obs = licor.NPQ(par_idx,i);
    observed.fs_obs = licor.Fs(:,i)./licor.Fo(par_idx,i);
    observed.fmp_obs = licor.Fm_(:,i)./licor.Fo(par_idx,i);

    data.Qin = licor.Qin(2:end,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps,1);                 % Atmospheric O2, mbar

    % Static
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
    v.Ku2 = 2e09;
    %     v.eps1 = 0.1;
    %     v.eps2 = 0.2;
    %     v.alpha_opt = 'dynamic';
    %     v.solve_xcs = solve_xcs;
    m = model_fun_c3(v);

    plot(m.Fs_a./m.Fo_a)
    hold on
    plot(observed.fs_obs,'r')
    pause
end

%% We do the inversion on light experiement
vcmax_0 =  50;
vqmax_0 = 175;
% ku2_0 = 2;
% a1_0 = 0.5;
% a2_0 = 0.5;
% eps1_0 = 0.5;
% eps2_0 = 0.25;
% vars_0 = [vcmax_0, vqmax_0,eps1_0,eps2_0]';
vars_0 = [vcmax_0, vqmax_0]';
% opts.LBounds = [10, 50,0,0]'; opts.UBounds = [150, 500,1,1]';
opts.LBounds = [0,0]'; opts.UBounds = [150,500]';
opts.Restarts=3;
opts.Noise.on=0;
weights.wAn = 1;
weights.wNPQ = 0;
weights.wFs = 0;

for i = 1:8
    %   for i = 1:1

    observed.an_obs = licor.A(:,i);
    observed.npq_obs = licor.NPQ(:,i);
    observed.fs_obs = licor.Fs(:,i)./licor.Fo(:,i);

    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps,1);                 % Atmospheric O2, mbar

    % Static
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
    v.Ku2 = 2e09;
    est_pars_static(:,i) = cmaes('chi2_photo_inv',vars_0,[],opts,v,observed,weights);

    %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    est_pars_dynamic(:,i) = cmaes('chi2_photo_inv',vars_0,[],opts,v,observed,weights);

end
% save est_pars_static_vcvq.mat est_pars_static
% save est_pars_dynamic_vcvq.mat est_pars_dynamic
disp("Inversion done!")
%%
close all
Fsize = 10;
fig1 = figure(1);
set(fig1,'units','inch','Position',[5 5 10 8],'color','w');
for i=1:8
    an_obs = licor.A(:,i);
    npq_obs = licor.NPQ(:,i);
    fs_obs = licor.Fs(:,i)./licor.Fo(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
    v.Ku2 = 2e09;
    free_pars.vcmax = est_pars_static(1,i);
    free_pars.vqmax = est_pars_static(2,i);
    sim_m_static=model_fun_c3_inv(v,free_pars);
    %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    free_pars.vcmax = est_pars_dynamic(1,i);
    free_pars.vqmax = est_pars_dynamic(2,i);
    sim_m_dynamic=model_fun_c3_inv(v,free_pars);

    set(0,'CurrentFigure',fig1)
    hold on
    subplot(3,3,i)
    plot(licor.Qin(:,i), sim_m_static.An_a*1e6,'k',LineWidth=1)
    hold on
    plot(licor.Qin(:,i), sim_m_dynamic.An_a*1e6,'b',LineWidth=1)
    hold on
    plot(licor.Qin(:,i),an_obs,'r',LineWidth=1)
    ylim([-5,35])
    xlabel("Ci [ppm]",fontsize=Fsize)
end

%%
load est_pars_static_vcvq.mat
load est_pars_dynamic_vcvq.mat
% vcmax_0 =  50;
% vqmax_0 = 175;
% ku2_0 = 2;
% a1_0 = 0.5;
% a2_0 = 0.5;
eps1_0 = 0.5;
eps2_0 = 0.25;
% vars_0 = [vcmax_0, vqmax_0,eps1_0,eps2_0]';
vars_0 = [eps1_0,eps2_0]';
% opts.LBounds = [10, 50,0,0]'; opts.UBounds = [150, 500,1,1]';
opts.LBounds = [0,0]'; opts.UBounds = [1,1]';
opts.Restarts=3;
opts.Noise.on=0;
weights.wAn = 0;
weights.wNPQ = 1;
weights.wFs = 0;

for i = 1:n_samples
    %   for i = 1:1

    observed.an_obs = licor.A(:,i);
    observed.npq_obs = licor.NPQ(:,i);
    observed.fs_obs = licor.Fs(:,i)./licor.Fo(:,i);

    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps,1);                 % Atmospheric O2, mbar

    % Static
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
    v.Ku2 = 2e09;

    v.vcmax = est_pars_static(1,i);
    v.vqmax = est_pars_static(2,i);

    est_pars_static_eps(:,i) = cmaes('chi2_photo_inv',vars_0,[],opts,v,observed,weights);

    %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    est_pars_dynamic_eps(:,i) = cmaes('chi2_photo_inv',vars_0,[],opts,v,observed,weights);

end
% save est_pars_static_vcvq.mat est_pars_static
% save est_pars_dynamic_vcvq.mat est_pars_dynamic
disp("Inversion done!")


