% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2_vcmax(vars,v,an_obs)
vcmax = vars(1);
% vqmax = vars(2);
% ku2 = vars(3);
try
    sim_m=model_fun_c3c4_vcmax(vcmax,v);
    an_est = sim_m.An_a.*1e6;
    chi2=norm(an_obs-an_est);
catch
    disp('Something went wrong; prob vcmax too low for dynamic opt of a1&a2')
    chi2 = NaN;
end
 

