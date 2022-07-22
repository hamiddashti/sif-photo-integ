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

licor = obs.aus_data.licor_co2;
n_steps = size(licor.Ci,1);
n_samples = size(licor.Ci,2)
p = 970;
ppm2bar     =  1e-6 .* (p .*1E-3);
symsolver_fun();
names = obs.aus_data.names;

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
    v.Ku2 = 2e09;
    est_pars_static(:,i) = cmaes('chi2_photo_inv',vars_0,[],opts,v,observed,weights);
    
    %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    est_pars_dynamic(:,i) = cmaes('chi2_photo_inv',vars_0,[],opts,v,observed,weights);

end
save est_pars_static.mat est_pars_static
save est_pars_dynamic.mat est_pars_dynamic
disp("Inversion done!")


%% Run the model with estimated pars
load est_pars_static.mat
load est_pars_dynamic.mat
close all
Fsize = 10;
fig = figure (2);
set(gcf,'units','inch','Position',[25 25 10 8],'color','w');
for i=1:8
    an_obs = licor.A(:,i);
    npq_obs = licor.NPQ(:,i);
    data.Qin = licor.Qin(:,i);
    data.Tin = licor.Tleaf(:,i);
    data.Cin = licor.Ci(:,i).*ppm2bar;
    data.Oin = repmat(209,n_steps,1);                 % Atmospheric O2, mbar
    v = configure_fun_c3(data);
    v.CB6F = 175./v.kq.*1e-06;       % Cyt b6f density, mol sites m-2
    v.RUB = 50./v.kc.*1e-06;      % Rubisco density, mol sites m-2
    v.Rds = 0.05;
    v.Abs = chems_inversion.Abs(i);
    free_pars.vcmax = est_pars_static(1,i);
    free_pars.vqmax = est_pars_static(2,i);
    sim_m_static=model_fun_c3_inv(v,free_pars);
    %Dynamic
    v.alpha_opt = 'dynamic';
    v.solve_xcs = solve_xcs;
    free_pars.vcmax = est_pars_dynamic(1,i);
    free_pars.vqmax = est_pars_dynamic(2,i);
    sim_m_dynamic=model_fun_c3_inv(v,free_pars);
    
    hold on
    subplot(3,3,i)
    plot(licor.Ci(:,i), sim_m_static.An_a*1e6,'k',LineWidth=1)
    hold on
    plot(licor.Ci(:,i), sim_m_dynamic.An_a*1e6,'b',LineWidth=1)
    hold on
    plot(licor.Ci(:,i),an_obs,'r',LineWidth=1)
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
% save est_pars_static.mat est_pars_static
% save est_pars_dynamic.mat est_pars_dynamic
disp("Inversion done!")