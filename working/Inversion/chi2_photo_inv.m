% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2_photo_inv(vars,v,observed,weights)
free_pars.vcmax = vars(1);
free_pars.vqmax = vars(2);

an_obs = observed.an_obs;
npq_obs = observed.npq_obs;
fs_obs = observed.fs_obs;

wAn = weights.wAn;
wNPQ = weights.wNPQ;
wFs = weights.wFs;

try
    sim_m=model_fun_photo_inv(v,free_pars);
 
    an_sim = sim_m.An_a.*1e6;
    chi2=norm(an_obs-an_sim);
    
    %     Multiple objesctive function  
    %     npq_sim = sim_m.Kn2_a.*1e-9;
    %     npq_sim = sim_m.PAM9_a; 
    %     fs_sim  = sim_m.Fs_a./sim_m.Fo_a;
    %     chi2=norm(fs_obs-fs_sim);
    
    %     a1 = (wAn.*((an_obs-an_sim)./max(an_obs)).^2);
    %     a2 = (wFs.*((fs_obs-fs_sim)./max(fs_obs)).^2);
    %     a3 = (wNPQ.*((npq_obs-npq_sim)./max(npq_obs)).^2);
    %     chi2 = sum(a1+a2+a3);
catch
    disp('Something went wrong: usually est vcmax too low')
    chi2 = NaN;
end
 

