% A stepwise function to account for water stress. 

function y=water_stress(rwc,var,s,rwc_t,rwc_c)
% rwc:     relative water content [unitless: 0,1]
% A_0:     assimilation rate (umol m-2 s-1)
% s:       slope
% rwc_t:   rwc threshold (rwc after which reduction in rwc decreases A_0)
% rwc_c:   rwc after wich there is no photsythesis
    
    if nargin<5
        rwc_c = 0.4;
    end
    if nargin < 4
        rwc_t = 0.9;
    end

    if rwc >rwc_t
        y=var;
    elseif (rwc<=rwc_t)&&(rwc>rwc_c)
        y = var-var.*s.*(rwc_t-rwc);
    elseif (rwc <= rwc_c)
        y=0;
    end

end

