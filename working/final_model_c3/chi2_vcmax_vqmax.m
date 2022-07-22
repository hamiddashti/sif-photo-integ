% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2_vcmax_vqmax(vars,v,an_obs)
vcmax = vars(1);
vqmax = vars(2);
% ku2 = vars(3);
try
    sim_m=model_fun_c3_vcmax_vqmax(vcmax,vqmax,v);
    an_est = sim_m.An_a.*1e6;
    chi2=norm(an_obs-an_est);
catch
    disp('Something went wrong with this model run')
    chi2 = NaN;
end
 

