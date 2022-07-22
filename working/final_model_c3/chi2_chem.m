function chi2=chi2_chem(chems,refl_obs,tran_obs,wR,wT,noise,spectral,optipar)

leafbio.Cab  = chems(1);
leafbio.Cca = chems(2);
leafbio.Cdm = chems(3);
leafbio.Cw  = chems(4);
leafbio.N = chems(5);
% leafbio.Cant = chems(6);
% leafbio.Cx =  chems(7);
leafbio.Cant = 0;
leafbio.Cx =  0;

leafbio.Cs = 0.0;
leafbio.fqe = [0, 0];

wl_low = 400;
wl_high = 2000;
wl_obs = spectral.wl_ASD;
wl_sim = spectral.wlP;
I_obs = find(wl_obs>=wl_low & wl_obs<=wl_high); % data after 2000 nm is noisy
I_sim = find(wl_sim>=wl_low & wl_sim<=wl_high); % data after 2000 nm is noisy

refl_obs = refl_obs(I_obs);
tran_obs = tran_obs(I_obs);

[leafopt] = fluspect_B_CX(spectral,leafbio,optipar);
refl_no_noise = leafopt.refl;
tran_no_noise = leafopt.tran;

% In case sim and obs have different--> interpolate (in this study they
% have same wl)
refl_no_noise = refl_no_noise(I_sim); 
tran_no_noise = tran_no_noise(I_sim); 

if noise ==true  % Add 1% noise to simulations
    noiseSigma_R = 0.01 * refl_no_noise;
    noiseSigma_T = 0.01 * tran_no_noise;
    noiseR = noiseSigma_R .* (randn(1, length(refl_no_noise))');
    noiseT = noiseSigma_T .* (randn(1, length(tran_no_noise))');
    refl_sim_noisy = refl_no_noise + noiseR;
    tran_sim_noisy = tran_no_noise + noiseT;
%     chi2 = norm(wR.*(refl_obs-refl_sim_noisy)+wT.*(tran_obs-tran_sim_noisy)); 
    chi2 = sum((wR.*(refl_obs-refl_sim_noisy).^2)+(wT.*(tran_obs-tran_sim_noisy).^2));
else
%     chi2=norm(refl_obs-refl_no_noise);
    chi2 = norm(wR.*(refl_obs-refl_no_noise)+wT.*(tran_obs-tran_no_noise));  
end
