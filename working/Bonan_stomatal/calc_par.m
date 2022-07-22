function par = calc_par(E,wl)
    % calclate the par of incoming radiation 
    % E incoming radiation 
    % wl wavelengths associated to E

    % define the spectral regions
    wlPAR = 400:1:700; 
    IparM = find(wl>=400 & wl<=700);
    E_par = E(IparM);
    par = 1E-3*Sint(E_par,wlPAR)*4.565;

end
