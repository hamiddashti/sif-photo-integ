clear
clc
warning('off')

% Prepare observations
load obs.mat % Load observation
load spectral.mat
load optipar.mat
load chems_inversion.mat

licor_co2 = obs.aus_data.licor_co2;
licor_lrc = obs.aus_data.licor_lrc;
n_steps_co2 = size(licor_co2.Ci,1);
n_steps_lrc = size(licor_lrc.Qin,1);
n_samples = size(licor_co2.Ci,2);
p = 970;
ppm2bar     =  1e-6 .* (p .*1E-3);
symsolver_fun();
names = obs.aus_data.names;

%% Estimate Vcmax and Vqmax using CO2 experiment

% Initial values
vcmax_0 =  50;
vqmax_0 = 175;

vars_0 = [vcmax_0, vqmax_0]';

% Bounderies an options for CMAES inverse method
opts.LBounds = [10, 50]'; opts.UBounds = [100, 400]';
opts.Restarts=3;
opts.Noise.on=0;

% Weights (0 or 1) to determine what to include in the objective function
weights.wAn = 1;
weights.wNPQ = 0;
weights.wFs = 0;

for i = 1:n_samples
    observed.an_obs = licor_co2.A(:,i);
    observed.npq_obs = licor_co2.NPQ(:,i);
    observed.fs_obs = licor_co2.Fs(:,i)./licor_co2.Fo(:,i);

    data.Qin = licor_co2.Qin(:,i);                     % PAR PPFD
    data.Tin = licor_co2.Tleaf(:,i);
    data.Cin = licor_co2.Ci(:,i).*ppm2bar.*1e6         %Ci [ubar CO2]
    data.Oin = repmat(209,n_steps_co2,1);              % Atmospheric O2, mbar

    % Static
    v = configure_fun(data);
    % The following parameters are from Jen's paper
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);  % Absorption are estimated based on 1-Refl-trans
    v.Ku2 = 2e09;
    est_pars_static(:,i) = cmaes('chi2_photo_inv',vars_0,[],opts,v,observed,weights);
% 
%     Dynamic
%     v.alpha_opt = 'dynamic';
%     v.solve_xcs = solve_xcs;
%     est_pars_dynamic(:,i) = cmaes('chi2_photo_inv',vars_0,[],opts,v,observed,weights);

end
save est_pars_static.mat est_pars_static
% save est_pars_dynamic.mat est_pars_dynamic
disp("Inversion done!")

%% Plot goodness of fit

load est_pars_static.mat
load est_pars_dynamic.mat

close all
Fsize = 10;
fig1 = figure (1);
set(fig1,'units','inch','Position',[25 25 10 8],'color','w');
fig2 = figure (2);
set(fig2,'units','inch','Position',[25 25 10 8],'color','w');
for i=1:8
    observed.an_obs = licor_co2.A(:,i);
    observed.npq_obs = licor_co2.NPQ(:,i);
    observed.fs_obs = licor_co2.Fs(:,i)./licor_co2.Fo(:,i);

    data.Qin = licor_co2.Qin(:,i);                     % PAR PPFD
    data.Tin = licor_co2.Tleaf(:,i);
    data.Cin = licor_co2.Ci(:,i).*ppm2bar.*1e6         %Ci [ubar CO2]
    data.Oin = repmat(209,n_steps_co2,1);              % Atmospheric O2, mbar

    % Static
    v = configure_fun(data);
    % The following parameters are from Jen's paper
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);  % Absorption are estimated based on 1-Refl-trans
    v.Ku2 = 2e09;

    free_pars.vcmax = est_pars_static(1,i);
    free_pars.vqmax = est_pars_static(2,i);
    sim_m_static=model_fun_photo_inv(v,free_pars);
%     %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    free_pars.vcmax = est_pars_dynamic(1,i);
    free_pars.vqmax = est_pars_dynamic(2,i);
    sim_m_dynamic=model_fun_photo_inv(v,free_pars);
    
    set(0,'CurrentFigure',fig1)
    hold on
    subplot(3,3,i)
    plot(licor_co2.Ci(:,i), sim_m_static.An_a*1e6,'k',LineWidth=1)
    hold on
    plot(licor_co2.Ci(:,i), sim_m_dynamic.An_a*1e6,'b',LineWidth=1)
    hold on
    plot(licor_co2.Ci(:,i),observed.an_obs,'r',LineWidth=1)
    ylim([-5,35])
    xlabel("Ci [ppm]",fontsize=Fsize)
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
    sgtitle("CO2 experiment (An plot)")
    set(gca,"FontSize",12)

    set(0,'CurrentFigure',fig2)
    subplot(3,3,i)
    plot(licor_co2.Ci(:,i), sim_m_static.Fs_a./sim_m_static.Fo_a...
        ,'k',LineWidth=1)
    hold on
    plot(licor_co2.Ci(:,i), sim_m_dynamic.Fs_a./sim_m_dynamic.Fo_a...
        ,'b',LineWidth=1)
    hold on 
    plot(licor_co2.Ci(:,i),observed.fs_obs,'r')
    xlabel("Ci [ppm]",fontsize=Fsize)
    ylabel("Fs/Fo",fontsize=Fsize)
    sgtitle(" Co2 experiment, Fs/Fo")

end

newPosition = [0.5 0.03 0.01 0.01];
newUnits = 'normalized';
hL = legend(["static","dynamic","obs"],"orientation","horizontal","location","southoutside")
set(hL,'Position', newPosition,'Units', newUnits,"FontSize",12);
saveas(fig1,'./Figures/fit_quality_co2_an.png')
saveas(fig2,'./Figures/fit_quality_co2_fs.png')


%% Lets plot estimated Vcmax vs.Chl 

close all
fig3 = figure(3);
Fsize=12;
set(fig3,'units','inch','Position',[25 25 5 5],'color','w');

vcmax_static = est_pars_static(1,:);
vcmax_dynamic = est_pars_dynamic(1,:);
chl = chems_inversion.est_chems(1,1:8);
mdl_static = fitlm(chl,vcmax_static);
scatter(chl,vcmax_static,Marker="+",LineWidth=2,MarkerEdgeColor ="k")
h1 = lsline;
h1.LineStyle = "--"
h1.LineWidth = 2;
hold on
scatter(chl,vcmax_dynamic,Marker="+",LineWidth=2,MarkerEdgeColor ="r")
box on
% plotAdded(mdl,[],Markersize=12,color="k",LineWidth=3)

xlabel('Chl [\mug cm^2] ','interpreter','Tex')
ylabel('Vcmax [mol CO^2 m^{-2} s^{-1}]','interpreter','Tex')
eq1 = "Vcmax = " + round(mdl_static.Coefficients.Estimate(1),2) + "+" + ...
    round(mdl_static.Coefficients.Estimate(2),2) + "Chl";
eq2 = "R^2 = " + round(mdl_static.Rsquared.Ordinary,2);
text(35,58, eq1,FontSize=Fsize)
text(35,52, eq2,FontSize=Fsize)
legend(["Static","","Dynamic"],"Location",'southeast')
saveas(gcf,"./Figures/Vcmax_chl.png")

%% Now run the model for LRC data
close all
fig4 = figure(4);
Fsize=12;
set(fig4,'units','inch','Position',[25 25 8 8],'color','w');

fig5 = figure(5);
Fsize=12;
set(fig5,'units','inch','Position',[25 25 8 8],'color','w');

for i = 1:8

    an_obs = licor_lrc.A(:,i);
    fs_obs = licor_lrc.Fs(:,i)./licor_lrc.Fo(:,i);

    data.Qin = licor_lrc.Qin(:,i);
    data.Tin = licor_lrc.Tleaf(:,i);
    data.Cin = licor_lrc.Ci(:,i).*ppm2bar.*1e6;
    data.Oin = repmat(209,n_steps_lrc,1);          % Atmospheric O2, mbar

    % Static
    v = configure_fun(data);
    % The following parameters are from Jen's paper
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);  % Absorption are estimated based on 1-Refl-trans
    v.Ku2 = 2e09;

    free_pars.vcmax = est_pars_static(1,i);
    free_pars.vqmax = est_pars_static(2,i);
    sim_lrc_static=model_fun_photo_inv(v,free_pars);
    %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    free_pars.vcmax = est_pars_dynamic(1,i);
    free_pars.vqmax = est_pars_dynamic(2,i);
    sim_lrc_dynamic=model_fun_photo_inv(v,free_pars);
    
    set(0,'CurrentFigure',fig4)
    subplot(3,3,i)
    plot(data.Qin,an_obs,'r')
    hold on 
    plot(data.Qin,sim_lrc_static.An_a.*1e6,'b')
    hold on
    plot(data.Qin,sim_lrc_dynamic.An_a.*1e6,'k')
    xlabel("PAR")
    ylabel("An")
    sgtitle("LRC experiment (An)")

    set(0,'CurrentFigure',fig5)
    subplot(3,3,i)
    plot(data.Qin,fs_obs,'r')
    hold on
    plot(data.Qin,sim_lrc_static.Fs_a./sim_lrc_static.Fo_a,'b')
    hold on
    plot(data.Qin,sim_lrc_dynamic.Fs_a./sim_lrc_static.Fo_a,'k')
    xlabel("PAR")
    ylabel("Fs/Fo")
    sgtitle("LRC experiment (Fs)")

end




%% Run the model with estimated pars
load est_pars_static.mat
load est_pars_dynamic.mat
close all
Fsize = 10;
fig = figure (2);
set(gcf,'units','inch','Position',[25 25 10 8],'color','w');
for i=1:8
    an_obs = licor_co2.A(:,i);
    npq_obs = licor_co2.NPQ(:,i);
    data.Qin = licor_co2.Qin(:,i);
    data.Tin = licor_co2.Tleaf(:,i);
    data.Cin = licor_co2.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps_co2,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
    free_pars.vcmax = est_pars_static(1,i);
    free_pars.vqmax = est_pars_static(2,i);
    sim_lrc_static=model_fun_c3_inv(v,free_pars);
    %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    free_pars.vcmax = est_pars_dynamic(1,i);
    free_pars.vqmax = est_pars_dynamic(2,i);
    sim_m_dynamic=model_fun_c3_inv(v,free_pars);

    hold on
    subplot(3,3,i)
    plot(licor_co2.Ci(:,i), sim_lrc_static.An_a*1e6,'k',LineWidth=1)
    hold on
    plot(licor_co2.Ci(:,i), sim_m_dynamic.An_a*1e6,'b',LineWidth=1)
    hold on
    plot(licor_co2.Ci(:,i),an_obs,'r',LineWidth=1)
    ylim([-5,35])
    xlabel("Ci [ppm]",fontsize=Fsize)
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

end
% Lgnd = legend(["static","dynamic","obs"])
% Lgnd.Position(1) = 1.1;
% Lgnd.Position(2) = 0.2;
% Lgnd.FontSize=14;

newPosition = [0.5 0.03 0.01 0.01];
newUnits = 'normalized';
hL = legend(["static","dynamic","obs"],"orientation","horizontal","location","southoutside")
set(hL,'Position', newPosition,'Units', newUnits,"FontSize",12);


vcmax = est_pars_static(1,:);
chl = chems_inversion.est_chems(1,:);
mdl = fitlm(chl,vcmax);
subplot(3,3,9)
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
saveas(gcf,'./Figures/fit_Vcmax_Vqmax.png')

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
opts.LBounds = [10, 50]'; opts.UBounds = [150, 500]';
opts.Restarts=3;
opts.Noise.on=0;
weights.wAn = 1;
weights.wNPQ = 0;
weights.wFs = 0;

% for i = 1:n_samples
for i = 1:1

    observed.an_obs = licor_co2.A(:,i);
    observed.npq_obs = licor_co2.NPQ(:,i);
    observed.fs_obs = licor_co2.Fs(:,i)./licor_co2.Fo(:,i);

    data.Qin = licor_co2.Qin(:,i);
    data.Tin = licor_co2.Tleaf(:,i);
    data.Cin = licor_co2.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps_co2,1);                 % Atmospheric O2, mbar

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
% save est_pars_static.mat est_pars_static
% save est_pars_dynamic.mat est_pars_dynamic
disp("Inversion done!")