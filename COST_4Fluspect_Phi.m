function [ertot,leafopt2] = COST_4Fluspect_Phi(params,measurement,input,P)
leafbio_all = input{1};
optipar     = input{2};
spectral    = input{3};
target      = input{4};

optipar.phi = interp1(P.CoeffsKcaLength,params',spectral.wlP)';%params';
optipar.phi(isnan(optipar.phi)) = 0;


[leafopt2.Fu,leafopt2.Fd] =deal(zeros(length(spectral.Fqe),1));

if target ==0
    er          = zeros(2*length(spectral.Fqe),size(measurement.refl,2));
else
    
    er          = zeros(length(spectral.Fqe),size(measurement.refl,2));
end

for k = 1:length(leafbio_all.Cab)
    
    leafbio.Cab = leafbio_all.Cab(k);
    leafbio.Cca = leafbio_all.Cca(k);
    leafbio.Cdm = leafbio_all.Cdm(k);
    leafbio.Cw = leafbio_all.Cw(k);
    leafbio.Cs = leafbio_all.Cs(k);
    leafbio.N = leafbio_all.N(k);
    leafbio.Cx = leafbio_all.Cx(k);
    leafbio.Cant = leafbio_all.Cant(k);
    leafbio.fqe = 0.01;
    
    [leafopt] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);
    
    E = (measurement.E(:,k))';
    jj = find(spectral.wlM>=spectral.Fqe(1) & spectral.wlM<=spectral.Fqe(end));

    switch target
        case 0
            leafopt2.Fu (:,k)         = leafopt.Mb(P.intersect,:) * E';
            leafopt2.Fd (:,k)         = leafopt.Mf(P.intersect,:) * E';
            er(:,k) = [(measurement.Fu(jj,k)-leafopt2.Fu(:,k))*(target<2);(measurement.Fd(jj,k)-leafopt2.Fd(:,k))*~(target==1)];
            %er = er(~isnan(er));
        case 1
            leafopt2.Fu (:,k)         = leafopt.Mb(P.intersect,:) * E;
            er(:,k) = (measurement.Fu(jj,k)-leafopt2.Fu)*(target<2);
            %er = er(~isnan(er));
        case 2
            leafopt2.Fd  (:,k)        = leafopt.Mf(P.intersect,:) * E;
            er(:,k) = (measurement.Fd(jj,k)-leafopt2.Fd)*~(target==1);
            %er = er(~isnan(er));
    end
%     er2 = P.PHI-optipar.FQE;
end
er_t1=er(:);
ertot = er_t1(~isnan(er_t1));
% er_t2=er2(:);
% ertot = [er_t1;er_t2];
