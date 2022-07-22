% Prepare input data for the integrated model

function [obs, spectral] = prepare_inputs(data_dir)
%% Prepare spectral data used by FLUSPECT-CX

spectral = define_bands;

%prepare_inputs(data_dir, leafbio,leaf_opt,spectral,env)

%% Prepare observation
%Prepare William's Woodgate data from Australia

disp("Preparing William's Woodgate data from Australia")
licor_dir =  [data_dir 'william_woodgate/LI6800/'];
asd_dir =  [data_dir 'william_woodgate/ASD/'];
qep_dir =  [data_dir 'william_woodgate/QEP/'];
obs.aus_data = prepare_aus_data(licor_dir,asd_dir,qep_dir);
% Prepare Nasstasia Vilfan's data
disp("Preparing Nasstasia Vilfan's data")
vilfan_dir = [data_dir 'Nastassia_Vilfan/'];
obs.vilfan_data = prepare_vilfan_data(vilfan_dir);
save obs.mat obs
save spectral.mat spectral

%% Prepare the photsynthesis inputs

% v = configure_fun(env);
% v.vcmax = (1.30 *leafbio.Cab+3.72)*1e-6;
%
% % Define the spectral regions that will be used in Fluspect
% leaf_opt.wlp      = spectral.wlP;         % PROSPECT wavelengths as a column-vector
% leaf_opt.wlPAR = spectral.wlPAR';
% leaf_opt.IparP    = find(wlp>=400 & wlp<=700); % Indices for PAR wavelenghts within wl
% leaf_opt.wlm = spectral.wlM;                   % measured wavelengths
% leaf_opt.IparM = find(wlm>=400 & wlm<=700);

%
end
