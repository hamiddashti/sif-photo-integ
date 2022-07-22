% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2_vcmax_ku2(vars,v,an_obs)
vcmax = vars(1);
ku2 = vars(2);
% ku2 = vars(3);
try
    sim_m=model_fun_c3_vcmax_ku2(vcmax,ku2,v);
    an_est = sim_m.An_a.*1e6;
    chi2=norm(an_obs-an_est);
catch
    disp('Something went wrong with this model run')
    chi2 = NaN;
end
 

