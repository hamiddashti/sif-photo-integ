
% first run the photosynthesis model
m = model_fun(v);

% Now run the Fluspect-CX model for combined version 
[leafopt_test] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);
E = (measured.E(iwle,k));
Fu_E = leafopt_test.Mb*E;
Fd_E = leafopt_test.Mf*E;

