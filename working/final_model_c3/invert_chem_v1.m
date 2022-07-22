clear
clc
warning('off')
% data_dir = '/home/hamid/SIF/sif-photo-integ/data/';
load obs
spectral = define_bands()
% [obs, spectral] = prepare_inputs(data_dir);

load optipar.mat
load est_vcmax_vqmax_stat_co2.mat
ref = obs.aus_data.asd.refl;
tran = obs.aus_data.asd.trans;
n_samples = size(ref,2);


%% Inversion
% Initial values
Cab_0 =  50;
Cca_0 = 10;
Cdm_0 = 0.012;
Cw_0 = 0.009;
N_0 = 1.4;
% Cant_0 = 0.0;
% Cx_0 = 0.0;

% chems_0 = [Cab_0, Cca_0,Cdm_0,Cw_0,N_0,Cant_0,Cx_0]';
chems_0 = [Cab_0, Cca_0,Cdm_0,Cw_0,N_0]';

% opts.LBounds = [0,0,0,0,1,0,0]'; opts.UBounds = [100,30,0.5,0.4,4,10,1.5]';
opts.LBounds = [0,0,0,0,1]'; opts.UBounds = [100,30,0.5,0.4,4]';

opts.Restarts=3;
opts.Noise.on=1;
opt.CMA.active=1;
noise = true;
wR = 1;
wT=1;

for i = 1:n_samples
    refl_obs = ref(:,i);
    tran_obs = tran(:,i);

    [est_chems(:,i), fmin] = cmaes('chi2_chem',chems_0,[],opts,refl_obs,...
        tran_obs,wR,wT,noise,spectral,optipar);
end
chems_inversion.est_chems = est_chems;
disp("Inversion done!")


%%
% Test inversion
wl_low = 400;
wl_high = 2000;
wl_obs = spectral.wl_ASD;
wl_sim = spectral.wlP;
I_obs = find(wl_obs>=wl_low & wl_obs<=wl_high); % data after 2000 nm is noisy
I_sim = find(wl_sim>=wl_low & wl_sim<=wl_high); % data after 2000 nm is noisy
IparP    = find(wl_sim>=400 & wl_sim<=700);
wlPAR = spectral.wlPAR;

for i =1:n_samples
    refl_obs = ref(:,i);
    tran_obs = tran(:,i);
    
    leafbio_std.Cab = est_chems(1,i);
    leafbio_std.Cca = est_chems(2,i);
    leafbio_std.Cdm = est_chems(3,i);
    leafbio_std.Cw = est_chems(4,i);
    leafbio_std.N = est_chems(5,i);
    leafbio_std.Cant = 0;
    leafbio_std.Cx = 0;
    leafbio_std.Cs = 0;
    leafbio_std.fqe=[0, 0];

    [leafopt_sim] = fluspect_B_CX(spectral,leafbio_std,optipar);
    fig1 = figure(1);
    set(gcf,'unit','inch','Position',[5,5,16,9],'color','w');
   
    ax = subplot(3,3,i)
    yyaxis left
    plot(wl_obs(I_obs),refl_obs(I_obs),'r')
    hold on 
    plot(wl_sim(I_sim),leafopt_sim.refl(I_sim),'r')
    ylim([0,1])
    ylabel("Reflectance")
    yyaxis right
    plot(wl_obs(I_obs),tran_obs(I_obs),'b')
    ylabel("Transmitance")
    hold on
    plot(wl_sim(I_sim),leafopt_sim.tran(I_sim),'b')
    ax.YDir = "reverse"
    ylim([0,1])
    xlabel("wl")
    xlim([400,2000])
    left_color = [1 0 0];
    right_color = [0 0 1];
    set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
    
end
saveas(ax,"./Figures/Sim_ref_tran.png")


%% Calculate leaf absorptance 
for i = 1:n_samples    
    leafbio.Cab  = est_chems(1,i)
    leafbio.Cca = est_chems(2,i);
    leafbio.Cdm = est_chems(3,i);
    leafbio.Cw  = est_chems(4,i);
    leafbio.N = est_chems(5,i);
    leafbio.Cant = 0;
    leafbio.Cx =  0;
    leafbio.Cs = 0.0;
    leafbio.fqe = [0, 0];
    
    [leafopt] = fluspect_B_CX(spectral,leafbio,optipar);
    
    % calculate absorptance
    leafopt.absorb = 1-leafopt.refl-leafopt.tran;
    
    absorbed= 1E-3*Sint(leafopt.absorb(IparP),wlPAR)*4.565;
    reflected = 1E-3*Sint(leafopt.refl(IparP),wlPAR)*4.565;
    transmitted = 1E-3*Sint(leafopt.tran(IparP),wlPAR)*4.565;
    total_par = absorbed+reflected+transmitted;
    Abs(i) =absorbed/total_par;
end
chems_inversion.Abs = Abs;
save chems_inversion.mat chems_inversion
