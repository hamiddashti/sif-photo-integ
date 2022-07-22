% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2_photo_lrc(vars,v,observed,weights)
free_pars.vcmax = vars(1);
% free_pars.vqmax = vars(2);
free_pars.ku2 = vars(2);
% free_pars.a1 = vars(4);
% free_pars.a2 = vars(5);
% free_pars.eps1 = vars(1);
% free_pars.eps2 = vars(2);

an_obs = observed.an_obs;
npq_obs = observed.npq_obs;
fs_obs = observed.fs_obs;

wAn = weights.wAn;
wNPQ = weights.wNPQ;
wFs = weights.wFs;

% if (free_pars.eps1+free_pars.eps2)>1
%     disp("eps > 1");
%     chi2 = NaN;
%     return
% end
% if((free_pars.a1+free_pars.a2)>1)
%     disp("a > 1");
%     chi2 = NaN;
%     return
% end

try
    sim_m=model_fun_c3_inv_lrc(v,free_pars);

    an_sim = sim_m.An_a.*1e6;
    chi2=norm(an_obs-an_sim);

    npq_sim = sim_m.Kn2_a.*1e-9;
    npq_sim = sim_m.PAM9_a;
    fs_sim  = sim_m.Fs_a./sim_m.Fo_a;
    %     chi2=norm(fs_obs-fs_sim);

%     a1 = (wAn.*((an_obs-an_sim)./max(an_obs)).^2);
%     a2 = (wFs.*((fs_obs-fs_sim)./max(fs_obs)).^2);
%     a3 = (wNPQ.*((npq_obs-npq_sim)./max(npq_obs)).^2);
%     chi2 = sum(a1+a2+a3);

catch
    disp('Something went wrong with this model run')
    chi2 = NaN;
end


