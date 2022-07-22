function leafopt = main_flu(spectral, leafbio,optipar,type,E,chlf)
% type: "not_combined" or "combined". If not combined we use the fluspect
% with both photosystems if combined we use the fluspect where photosystems
% are combined. 
% E is incoming energy

if type == "not_combined"
    if length(leafbio.fqe) < 2
        disp("Provide both quantum effiencies!")
        return
    end
    [leafopt] = fluspect_B_CX(spectral,leafbio,optipar);
    leafopt.Fu = (leafopt.MbI+leafopt.MbII)* interp1(spectral.wlM,E,spectral.wlE)';
    leafopt.Fd = (leafopt.MfI+leafopt.MfII)* interp1(spectral.wlM,E,spectral.wlE)';

elseif type == "combined"
    if length(leafbio.fqe) > 1
        disp("Provide one quantum effiency!")
        return
    end
    [leafopt] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);
    
    if chlf == true
    leafopt.Fu         = leafopt.Mb * interp1(spectral.wlM,E,spectral.wlE)';
    leafopt.Fd         = leafopt.Mf * interp1(spectral.wlM,E,spectral.wlE)';
    end
 end
end



