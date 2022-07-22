function [out_opt out_photo] = main_fun_v1(spectral,leafbio,optipar,E,v)

    leafbio.fqe=[0];
    leafoptics = main_flu(spectral, leafbio,optipar,"combined",E,false);
    % [leafopt] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);
    
    % define the spectral regions
    wlp      = spectral.wlP;         % PROSPECT wavelengths as a column-vector
    wlPAR = spectral.wlPAR';
    IparP    = find(wlp>=400 & wlp<=700); % Indices for PAR wavelenghts within wl
    wlm = spectral.wlM;                   % measured wavelengths 
    IparM = find(wlm>=400 & wlm<=700);
    
    % calculate absorptance
    refl = leafoptics.refl;
    tran = leafoptics.tran;
    absorb = 1-refl-tran;
    
    E_par = E(IparM);
    absorbed= 1E-3*Sint(absorb(IparP).*E_par,wlPAR)*4.565;
    reflected = 1E-3*Sint(refl(IparP).*E_par,wlPAR)*4.565;
    transmitted = 1E-3*Sint(tran(IparP).*E_par,wlPAR)*4.565;
    total_par = absorbed+reflected+transmitted;
    Abs =absorbed/total_par;
        
    % Another way that used in SCOPE:
    % absorbed = 1E6*(0.001*Sint(e2phot(wlPAR*1E-9,par_wave.*absorb,constants),wlPAR)); %umol m-2 s-1
    % reflected = 1E6*(0.001*Sint(e2phot(wlPAR*1E-9,par_wave.*refl,constants),wlPAR)); %umol m-2 s-1
    % transmitted = 1E6*(0.001*Sint(e2phot(wlPAR*1E-9,par_wave.*tran,constants),wlPAR)); %umol m-2 s-1
    % total_par = absorbed+reflected+transmitted
    
    % Run the photosynthesis model to get the phi_F and CX
%     data.Qin = total_par;   % PAR, umol PPFD m-2 s-1
%     data.Tin = photo_data.Tin;                 % Leaf temperature, C
%     data.Cin = photo_data.Cin;                % Mesophyll CO2, ubar
%     data.Oin = photo_data.Oin;                 % Atmospheric O2, mbar
%     v = configure_fun(data);
%   
    
    v.Abs = Abs;
    m = model_fun(v);
    phi2F = m.phi2F_a;

    % phi1F = m.phi1F_a;
%     a1 = m.a1;
%     a2 = m.a2;
%     An = m.An_a;
    NPQ = m.PAM9_a;
    Cx = 0.3187.*NPQ;
    Cx(Cx>1.5) = 1.5;
    leafbio.fqe=[phi2F];
    leafbio.Cx = Cx;
    out_opt = main_flu(spectral, leafbio,optipar,"combined",E,true);

    out_photo.v = v;
    out_photo.m = m;
%     out_photo.NPQ = NPQ;
    out_photo.Cx = Cx;
    out_photo.Abs = Abs;
end
