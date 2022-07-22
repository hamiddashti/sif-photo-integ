function par_um = par2wave(spectral,par)

T = 5777;   %Sun temp [K]
wav = spectral.wlPAR.*10^-9;
% wav = (400:1:700)*10^-9; %Wavelenght [m]
weights = planck_weights(wav,T); %weights based on plancks law
alpha = 1/4.565; %converting par [umol sr-1 m-2] to [w m-2] 
par_w = par.*alpha; 
par_nm = par_w.*weights;
par_um = par_nm.*1000; % convert to W m-2 um-1
par_um = par_um';
end
