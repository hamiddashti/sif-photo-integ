% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2_vcmax_vqmax_eps(vars,v,an_obs,fs_obs,NPQ_obs,wA,wP,wNPQ)
vcmax = vars(1);
vqmax = vars(2);
Eps1 = vars(3);
Eps2 = vars(4);
try
    sim_m=model_fun_c3c4_vcmax_vqmax_eps(vcmax,vqmax,Eps1,Eps2,v);
    an_est = sim_m.An_a.*1e6;
    fs_est = sim_m.Fs_a./sim_m.Fo_a;
%     NPQ_est = sim_m.PAM3_a;
    NPQ_est = sim_m.Kn2_ma.*1e-9;
%     chi2=norm(fs_obs-fs_est);

    chi2 = sum((wA.*(an_obs-an_est).^2)+...
        (wP.*(fs_obs-fs_est).^2)+...
        (wNPQ.*(NPQ_obs-NPQ_est).^2));
catch
    disp('Something went wrong!')
    chi2 = NaN;
end
 

