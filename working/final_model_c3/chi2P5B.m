% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2P5B(vars,v,an_obs)
vcmax = vars(1);
vqmax = vars(2);
ku2 = vars(3);
try
    sim_m=model_fun_c3_vcmax(vcmax,vqmax,ku2,v);
    an_est = sim_m.An_a.*1e6;
    chi2=norm(an_obs-an_est);
catch
    chi2 = NaN;
end
 

