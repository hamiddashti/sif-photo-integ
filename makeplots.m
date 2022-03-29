% This script creates figures for the paper 'The scattering and
% re-absorption of chlorophyll fluorescence in plant leaves and canopy' by
% Van der Tol et al., submitted 10 July to RSE.

% please use the code: do_analysis.m to prepare the results for plotting.

%% some plotting options
LB = 660;       % lower wl limit of plots of measured fluorescence
UB = 800;       % upper wl limit of plots of measured fluorescence
ab = 'a':'z';   % alphabet
MS = 5;         % marker size
samples = [17,65,28,19,33]; % leaf samples for figure 3

%% load raw data for plotting. Forward simulations still need to be done.
% the folowing loads data were the retrievals are stored. These have been made with
% 'do_analysis.m
spectral = define_bands;
iwle = (400:750)-349;
datdir = '..\data\output\fluspect_output\2015\2017-11-28-1426\';

load([datdir 'spectral.mat']);
load([datdir 'measured.mat']);
load([datdir 'leafopt.mat']);
load([datdir 'Tsimulated.mat']);
load([datdir 'Rsimulated.mat']);
load([datdir 'leafbio.mat']);
load([datdir 'optipar.mat']);

spectral.iwlF = (640:850)-399;

I = find(leafbio_all.Cab>20 & leafbio_all.Cab<50);

%% carry out the forward simulations
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
    
    if k == 1
        fluspect_B_CX_plot_doubling(spectral,leafbio,optipar);
    end
    
    E = (measured.E(:,k))';
    leafopt_all.Fu (:,k)         = leafopt.Mb * E(51:401)';
    leafopt_all.Fd (:,k)         = leafopt.Mf * E(51:401)';
    
    Fmeas = interp1(spectral.wlM,(measured.Fu(:,k)+measured.Fd(:,k)),spectral.wlF(30:end))';
    %fqe_opt = (leafopt_all.Fu(30:end,k) + leafopt_all.Fd(30:end,k))\Fmeas;%.*leafbio.fqe;
    
    %leafopt_all.Fu(:,k) = leafopt_all.Fu(:,k)*fqe_opt;
    %leafopt_all.Fd(:,k) = leafopt_all.Fd(:,k)*fqe_opt;
    
    leafopt_all.refl(:,k) = leafopt.refl;
    leafopt_all.tran(:,k) = leafopt.tran;
    leafopt_all.kChlrel(:,k) = leafopt.kChlrel;
    

    leafopt_all.Etot(k) = 1E-3*Sint(measured.E(iwle,k),wl(iwle));

end

%% This produces Figure 3 in the paper. Example of spectra of 5 leaves.
f3 = figure(3); clf
set(f3,'Position',[174.6000 224.2000 967.2000 537.6000])
wl = spectral.wlM;
I3 = find(spectral.wlM>LB&spectral.wlM<UB);
spfig3 = zeros(6,length(samples));
for I = 1:length(samples)
    k = samples(I);
    iwle = (400:750)-349;
    
    Etot = 1E-3*Sint(measured.E(iwle,k),wl(iwle));
    Etot_absorbed = 1E-3*( (1-leafopt.refl(1:350)-leafopt.tran(1:350)) .* leafopt_all.kChlrel(1:350,k))'*measured.E(51:400,k); %W m-2
    
    spfig3(1,I) = subplot(6,5,I);
    z = plot(wl,[measured.refl(:,k) 1-measured.tran(:,k)],'k'); hold on
    set(z(:),'Color','r','LineWidth',2)
    set(gca,'ylim',[0 1],'xlim',[400 2450])
    plot(spectral.wlP,[leafopt_all.refl(:,k),1-leafopt_all.tran(:,k)],'k')
    title(['C_{ab} = ' num2str(leafbio_all.Cab(k),3)]);% ', C_{dm} = ' num2str(leafbio_all.Cdm(k),3)])
    
    spfig3(2,I) = subplot(6,5,I+5);
    z= plot(wl(I3),measured.Fu(I3,k)./Etot,'k'); hold on
    set(z(:),'Color','r','LineWidth',2)
    plot(spectral.wlF,leafopt_all.Fu(:,k)./Etot,'k')
    set(gca,'ylim',[0 .033])
    
    spfig3(3,I) = subplot(6,5,I+10);
    z= plot(wl(I3),measured.Fd(I3,k)./Etot,'k'); hold on
    set(z(:),'Color','r','LineWidth',2)
    plot(spectral.wlF,leafopt_all.Fd(:,k)./Etot,'k')
    set(gca,'ylim',[0 .033])
    
    spfig3(4,I) =  subplot(6,5,I+15);
    z = plot(wl,measured.refl(:,k)./(measured.tran(:,k)+measured.refl(:,k)),'k'); hold on
    set(z(:),'Color','r','LineWidth',2)
    plot(spectral.wlP,leafopt_all.refl(:,k)./(leafopt_all.refl(:,k)+leafopt_all.tran(:,k)),'k')
    set(gca,'ylim',[0.4 1])
    
    spfig3(5,I) =  subplot(6,5,I+20);
    z = plot(wl(wl>LB&wl<UB),measured.Fu(I3,k)./(measured.Fu(I3,k)+measured.Fd(I3,k)),'k'); hold on
    set(z(:),'Color','r','LineWidth',2)
    plot(spectral.wlF,leafopt_all.Fu(:,k)./(leafopt_all.Fu(:,k)+leafopt_all.Fd(:,k)),'k')
    set(gca,'xlim',[640 850])
    set(gca,'ylim',[0.4 1])
    
    spfig3(6,I) =  subplot(6,5,I+25);
    ym = measured.Fu(I3,k)./(measured.Fu(I3,k)+measured.Fd(I3,k)) ./ (measured.refl(I3,k)./(measured.tran(I3,k)+measured.refl(I3,k)));
    ys = leafopt_all.Fu(:,k)./(leafopt_all.Fu(:,k)+leafopt_all.Fd(:,k)) ./ (leafopt_all.refl(spectral.iwlF,k)./(leafopt_all.refl(spectral.iwlF,k)+leafopt_all.tran(spectral.iwlF,k)));
    z = plot(wl(I3),ym,'k'); hold on
    set(z(:),'Color','r','LineWidth',2)
    plot(spectral.wlF,ys,'k'); hold on
    set(gca,'xlim',[640 850])
    set(gca,'ylim',[0.7 1.4])
end
resizefigure(spfig3',5,6,0.1,0.1,0.05,0.05,0.97,0.93)
set(spfig3(2:end,:),'xlim',[640 850])

ylabel(spfig3(1,1),'\rho, 1-\tau')
ylabel(spfig3(2,1),'$\displaystyle\frac{F_b}{E_{tot}} |\mu m^{-1}| $','interpreter','latex')%,'(\mum^{-1})'])
ylabel(spfig3(3,1),'$\displaystyle\frac{F_f}{E_{tot}}|\mu m^{-1}| $','interpreter','latex')%,'(\mum^{-1})'])
ylabel(spfig3(4,1),'$\displaystyle\frac{\rho}{\rho + \tau}$','interpreter','latex')
ylabel(spfig3(5,1),'$\displaystyle\frac{F_b}{F_b+F_f}$','interpreter','latex')
ylabel(spfig3(6,1),'$\displaystyle\frac{F_b}{\rho}\frac{\rho+\tau}{F_b+F_f}$','interpreter','latex')%^{F_b}/_{\rho}')%\frac{\rho+\tau}{F_b+F_f}')%/(F_b+F_f) \newline / ( \rho/(\rho + \tau))'])

for k = 1:5
    xlabel(spfig3(6,k),'wl (nm)')
end

%% Model sensitivity study
% This produces Figures 6,7,8 and 9 in the paper
clear spfig2
leafbio0.Cab= 26.4173;
leafbio0.Cdm= 0.0094;
leafbio0.Cw = 0.0084;
leafbio0.Cs= 0.1723;
leafbio0.Cca= 5.7859;
leafbio0.N= 1.4727;
leafbio0.Cant = 0;
leafbio0.Cx= 0;
leafbio0.fqe= 0.0115;

E = measured.E(50:400,end);
Etot = 1E-3*Sint(E,wl(50:400));
[Fu_,Fd_] = deal(zeros(211,20));

f = zeros(4,1);
for param = 1:4
    leafbio = leafbio0;
    f(param) = figure(param+5); clf
    set(f(param),'Position',[488.2000 369 741.6000 392.8000])
    for k = 1:20
        switch param
            case 1, leafbio.Cant = 0+k/5;%%leafbio.Cca = 1+k;%leafbio.N = 1+k/10;
            case 2, leafbio.Cab = 5+k*4.75;
            case 3, leafbio.Cdm = k/1000;
            case 4, leafbio.Cs = k/16;
        end
        leafopt = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar);
        
        Etot_absorbed = 1E-3*( (1-leafopt.refl(1:351)-leafopt.tran(1:351)) .* leafopt.kChlrel(1:351))'*E;
        Fem = 1E3*optipar.phi((640:850)-399)*Etot_absorbed*leafbio.fqe;
        Fu = leafopt.Mb*E;
        Fd = leafopt.Mf*E;
        
        if param == 1
            Fu_(:,k) = Fu;
            Fd_(:,k) = Fd;
        end
        
        spfig2(1) = subplot(2,3,1);
        z(1,k) = plot((640:850),Fem/Etot); hold on
        text(720,0.2*0.95,'a')
        
        spfig2(2) = subplot(2,3,2);
        z(2,k) = plot((640:850),(Fu+Fd)/Etot); hold on
        text(720,0.095,'b')
        
        spfig2(3) = subplot(2,3,3);
        z(3,k) = plot((640:850),(Fu+Fd)./Fem); hold on
        text(720,0.3+0.7*0.95,'c')
        
        spfig2(4) = subplot(2,3,4);
        z(4,k) = plot((640:850),(leafopt.refl((640:850)-399)./(leafopt.refl((640:850)-399)+leafopt.tran((640:850)-399)))); hold on
        text(750,0.3+0.7*0.95,'d')
        
        spfig2(5) = subplot(2,3,5);
        z(5,k) = plot((640:850),(Fu)./(Fu+Fd)); hold on
        text(720,0.3+0.7*0.95,'e')
        
        spfig2(6) = subplot(2,3,6);
        z(6,k) = plot((640:850), (leafopt.refl((640:850)-399)./(leafopt.refl((640:850)-399)+leafopt.tran((640:850)-399)))./((Fu)./(Fu+Fd))); hold on
        text(720,0.6+(1.2-0.6)*.95,'f')
        
        i675 = 675-639;
        i760 = 760-639;
        %ratio(k,param) = Fu(i675)./Fu(i760) * (Fu(i760)+Fd(i760))/(Fu(i675)+Fd(i675));
        %sumFem(k,param) = sum(Fem);
    end
    
    for k = 1:20
        set(z(:,k),'Color',[k 0 20-k]/20)
    end
    set(spfig2(:),'xlim',[640 850])
    set(spfig2(1),'ylim',[0 0.2]);
    set(spfig2(2),'ylim',[0 0.1]);
    set(spfig2(3),'ylim',[0 1]);
    set(spfig2(4),'ylim',[0.3 1]);
    set(spfig2(5),'ylim',[0.3 1]);
    set(spfig2(6),'ylim',[0.6 1.2]);
    ylabel(spfig2(1),'F_{leaf} /E_{tot} (\mum^{-1})')
    ylabel(spfig2(2),'(F_b + F_f)/E_{tot} (\mum^{-1})')
    ylabel(spfig2(3),'(F_b+F_f)/F_{leaf}')
    ylabel(spfig2(4),'\rho/(\rho+\tau)')
    ylabel(spfig2(5),'F_b/(F_b+F_f)')
    ylabel(spfig2(6),'\rho/F_b  * (F_b+F_f)/(\rho+\tau)')
    for k = 4:6
        xlabel(spfig2(k),'wl (nm)')
    end
end

%% figure 5. The RMSE and the r2 of modelled versus measured fluorescence per wavelength.
Fu_meas = measured.Fu(spectral.wlM>LB & spectral.wlM<UB,:);
Fd_meas = measured.Fd(spectral.wlM>LB & spectral.wlM<UB,:);
Fu_mod = leafopt_all.Fu(spectral.wlF>LB & spectral.wlF<UB,:);
Fd_mod = leafopt_all.Fd(spectral.wlF>LB & spectral.wlF<UB,:);

r_meas = measured.refl(spectral.wlM>LB & spectral.wlM<UB,:);
t_meas = measured.tran(spectral.wlM>LB & spectral.wlM<UB,:);
r_mod = leafopt_all.refl(spectral.wlP>LB & spectral.wlP<UB,:);
t_mod = leafopt_all.tran(spectral.wlP>LB & spectral.wlP<UB,:);

for k1 = 1:2
    switch k1
        case 1
            ratiomod = Fu_mod./(Fu_mod+Fd_mod);
            ratiomeas = Fu_meas./(Fu_meas+Fd_meas);
            ratiomeas(Fu_meas+Fd_meas<1) = NaN;
        case 2
            ratiomod = r_mod./(r_mod+t_mod);
            ratiomeas = r_meas./(r_meas+t_meas);
    end
    
    RMSEratio = sqrt( mean( (ratiomod - ratiomeas).^2 ,2));
    rr = zeros(length(Fu_mod),1);
    for k = 1:length(Fu_mod)
        r = corrcoef(ratiomod(k,:),ratiomeas(k,:));
        rr(k) = r(2).^2;
    end
    
    figure(5),
    subplot(211), z(1,k1) = plot(LB+1:UB-1,RMSEratio,'k'); hold on
    ylabel('RMSE')
    
    subplot(212),z(2,k1) = plot(LB+1:UB-1,rr,'k'); hold on
    xlabel('wl (nm)')
    ylabel('r^2')
end
set(z(:,2),'Color','b')
%% figure 4. Measured versus modelled fluorescence

I745 = find(spectral.wlM==745);
I685 = find(spectral.wlM==685);

Fu745m = measured.Fu(I745,:);
Fu685m = measured.Fu(I685,:);

Fd745m = measured.Fd(I745,:);
Fd685m = measured.Fd(I685,:);

Fd745 = leafopt_all.Fd(745-639,:);
Fd685 = leafopt_all.Fd(685-639,:);

Fu745 = leafopt_all.Fu(745-639,:);
Fu685 = leafopt_all.Fu(685-639,:);
I0  = find(Fd685>0.5);

f4 = figure(4); clf
set(f4,'Position',[373 259.4000 700 420.8000])
subplot(231)
plot(Fu745m./leafopt_all.Etot,Fu745./leafopt_all.Etot,'kx','markersize', MS), hold on
%plot(Fu745m(I(:)),Fu745(I(:)),'ko','markersize', MS), hold on
plot([0 0.04],[0 0.04],'k')
set(gca,'xlim',[0 0.04],'ylim', [0 0.04])
xlabel('modelled F_b (745)/E_{tot}')
ylabel('measured F_b (745)/E_{tot}')
r = corrcoef(Fu745m./leafopt_all.Etot,Fu745./leafopt_all.Etot);
r2 = r(2)^2;
text(0.02, 0.004, ['r^2 = ' num2str(r2,2)])

subplot(232)
plot(Fd745./leafopt_all.Etot,Fd745m./leafopt_all.Etot,'kx','markersize', MS), hold on
%plot(Fd745(I(:)),Fd745m(I(:)),'ko','markersize', MS), hold on
plot([0 .03],[0 0.03],'k')
set(gca,'xlim',[0 0.03],'ylim', [0 0.03])
xlabel('modelled F_f (745)/E_{tot}')
ylabel('measured F_f (745)/E_{tot}')
r = corrcoef(Fd745./leafopt_all.Etot,Fd745m./leafopt_all.Etot);
r2 = r(2)^2;
text(0.016, 0.003, ['r^2 = ' num2str(r2,2)])


subplot(233)
plot(Fu745./Fd745, Fu745m./Fd745m, 'kx','markersize', MS), hold on
%plot(Fu745(I(:))./Fd745(I(:)),Fu745m(I(:))./Fd745m(I(:)),'ko','markersize', MS), hold on
plot([0 3],[0 3],'k')
set(gca,'xlim',[0.8 3],'ylim', [0.8 3])
xlabel('modelled F_f/F_b (745)')
ylabel('measured F_f/F_b (745)')
r = corrcoef(Fu745./Fd745, Fu745m./Fd745m);
r2 = r(2)^2;
text( 1.95,1, ['r^2 = ' num2str(r2,2)])

subplot(234)
plot(Fu685./leafopt_all.Etot,Fu685m./leafopt_all.Etot,'kx','markersize', MS), hold on
%plot(Fu685(I(:)),Fu685m(I(:)),'ko','markersize', MS), hold on
plot([0 15],[0 15],'k')
set(gca,'xlim',[0 0.02],'ylim', [0 0.02])
xlabel('modelled F_b (685)/E_{tot}')
ylabel('measured F_b (685)/E_{tot}')
r = corrcoef(Fu685./leafopt_all.Etot,Fu685m./leafopt_all.Etot);
r2 = r(2)^2;
text(0.01,0.002,['r^2 = ' num2str(r2,2)])

subplot(235)
plot(Fd685./leafopt_all.Etot,Fd685m./leafopt_all.Etot,'kx','markersize', MS), hold on
%plot(Fd685(I(:)),Fd685m(I(:)),'ko','markersize', MS), hold on
plot([0 15],[0 15],'k')
set(gca,'xlim',[0 0.02],'ylim', [0 0.02])
xlabel('modelled F_f (685)/E_{tot}')
ylabel('measured F_f (685)/E_{tot}')
r = corrcoef(Fd685./leafopt_all.Etot,Fd685m./leafopt_all.Etot);
r2 = r(2)^2;
text(0.01,0.002, ['r^2 = ' num2str(r2,2)])

subplot(236)
plot(Fu685(:)./Fd685(:), Fu685m(:)./Fd685m(:), 'kx','markersize', MS), hold on
%plot(Fu685(I0)./Fd685(I0), Fu685m(I0)./Fd685m(I0), 'kx','markersize', MS), hold on
%plot(Fu685(I(:))./Fd685(I(:)),Fu685m(I(:))./Fd685m(I(:)),'ko','markersize', MS),
plot([0 10],[0 10],'k')
set(gca,'xlim',[0 10],'ylim', [0 10])
xlabel('modelled F_b/F_f (685)')
ylabel('measured F_b/F_f (685)')
r = corrcoef(Fu685(:)./Fd685(:), Fu685m(:)./Fd685m(:));
r2 = r(2)^2;
text(5,1, ['r^2 = ' num2str(r2,2)])

%% Figure A.11, old and new Optipar fluorescence spectra, plotted for leaf 28
k = 28;

leafbio2.Cab= leafbio_all.Cab(k);
leafbio2.Cdm= leafbio_all.Cdm(k);
leafbio2.Cw= leafbio_all.Cw(k);
leafbio2.Cs= leafbio_all.Cs(k);
leafbio2.Cca= leafbio_all.Cca(k);
leafbio2.N= leafbio_all.N(k);
leafbio2.Cx= leafbio_all.Cx(k);
leafbio2.Cant= leafbio_all.Cant(k);
leafbio2.fqe(1) = leafbio_all.fqe_opt(k)*0.2;
leafbio2.fqe(2) = leafbio_all.fqe_opt(k)*0.8;

leafbio3 = leafbio2;
leafbio3.fqe = leafbio_all.fqe_opt(k);
[leafopt_old] = fluspect_B_CX(spectral,leafbio2,optipar);
[leafopt_new] = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio3,optipar);

leafopt_old.Fu = (leafopt_old.MbI+leafopt_old.MbII)* interp1(spectral.wlM,measured.E(:,k),spectral.wlE)';
leafopt_old.Fd = (leafopt_old.MfI+leafopt_old.MfII)* interp1(spectral.wlM,measured.E(:,k),spectral.wlE)';

leafopt_new.Fu = (leafopt_new.Mb)* interp1(spectral.wlM,measured.E(:,k),spectral.wlE)';
leafopt_new.Fd = (leafopt_new.Mf)* interp1(spectral.wlM,measured.E(:,k),spectral.wlE)';
iwlE = find(spectral.wlM>399 & spectral.wlM<750);
Etot = 1E-3*Sint(measured.E(iwle,k),spectral.wlM(iwle));

Fum = interp1(spectral.wlM,measured.Fu(:,k),spectral.wlF);
Fum(spectral.wlF<660) = NaN;
Fdm = interp1(spectral.wlM,measured.Fd(:,k),spectral.wlF);
Fdm(spectral.wlF<660) = NaN;

figure(11), clf
s(1) = subplot(131);
z=plot(spectral.wlP,[ optipar.phiI optipar.phiII optipar.phi],'k');
set(z(1),'color',[.3 .3 .3])
set(z(3),'linewidth',2)
legend('PSI','PSII','new')
set(gca,'xlim',[640 850],'ylim',[0 0.03])
xlabel('wl (nm)')
ylabel('\phi')
text(710,0.027,'a','FontSize',12)
set(gca,'FontSize',12)

s(2) = subplot(132);
z=plot(spectral.wlF,[ Fum'/Etot  leafopt_old.Fu/Etot leafopt_new.Fu/Etot],'k');
set(z(1),'lineStyle','none','Marker','x','MarkerSize',3)
set(z(2),'color',[.3 .3 .3])
set(z(3),'linewidth',2)
set(gca,'xlim',[640 850],'ylim',[0 0.03])
ylabel('F_b/E_{tot} (\mum^{-1})')
xlabel('wl (nm)')
legend('measured','old','new')
text(710,0.027,'b','FontSize',12)
set(gca,'FontSize',12)

s(3) = subplot(133);
z=plot(spectral.wlF,[Fdm'/Etot leafopt_old.Fd/Etot leafopt_new.Fd/Etot  ],'k');
set(z(1),'lineStyle','none','Marker','x','MarkerSize',3)
set(z(2),'color',[.3 .3 .3])
set(z(3),'linewidth',2)
set(gca,'xlim',[640 850],'ylim',[0 0.03])
ylabel('F_f/E_{tot} (\mum^{-1})')
legend('measured','old','new')
xlabel('wl (nm)')
text(710,0.027,'c','FontSize',12)
set(gca,'FontSize',12)

resizefigure(s,3,1,.09,.22,.12,.07)

%% Figure 10. sensitivity analysis of TOC fluorescence and reflectance
f10 = figure(10); clf
set(f10,'Position',[488.2000 85 773.6000 676.8000])
parnames = {'N','C_{ab}', 'C_{dm}', 'C_s', 'LAI', 'LIDF', 'soil brightness'};

direc = '..\data\output\';
subdirecs ={    'sensitivity_N_2017-11-29-1444\'...
    'sensitivity_Cab_2017-11-29-1337\',...
    'sensitivity_Cdm_2017-11-29-1448\'...
    'sensitivity_Cs_2017-12-06-0616\'...
    'sensitivity_LAI_2017-11-29-1342\' ...
    'sensitivity_LIDF_2018-07-10-2057\' ...
    'sensitivity_soil_2019-02-19-1124\'};

for param = 1:7
    params = dlmread([direc subdirecs{param} 'pars_and_input_short.dat'],'',1,0);
    F = dlmread([direc subdirecs{param} 'fluorescence.dat'],'',2,0);
    Fem = dlmread([direc subdirecs{param} 'fluorescence_emitted_by_all_photosystems.dat'],'',2,0);
    
    r = dlmread([direc subdirecs{param} 'reflectance.dat'],'',2,0);
    wl = dlmread([direc subdirecs{param} 'wl.dat'],'',2,0);
    
    fesc = pi*F./Fem;
        
    for k = 1:length(params)
        spfig2(1,param) = subplot(7,4,1+(param-1)*4);
        z(1,k) = plot(wl,r(k,:)'); hold on
        %text(450,0.60,[ab(1+(param-1)*4) ' (range of ' parnames{param} ')'])
        text(450,0.60,[num2str(1+(param-1)*4) ' (range of ' parnames{param} ')'])
          
        spfig2(2,param) = subplot(7,4,2+(param-1)*4);
        z(2,k) = plot((640:850),F(k,:)'); hold on
        %text(650,3.6,ab(2+(param-1)*4))
        text(650,3.6,num2str(2+(param-1)*4))
        
        spfig2(3,param) = subplot(7,4,3+(param-1)*4);
        z(3,k) = plot((640:850),Fem(k,:)'); hold on
        %text(650,54,ab(3+(param-1)*4))
        text(650,54,num2str(3+(param-1)*4))
        
        spfig2(4,param) = subplot(7,4,4+(param-1)*4);
        z(4,k) = plot((640:850),fesc(k,:)'); hold on
        %text(650,0.54,ab(4+(param-1)*4)) 
        text(650,0.54,num2str(4+(param-1)*4))
    end
    for k = 1:length(params)
        set(z(:,k),'Color',[k 0 length(params)-k]/(length(params)))
    end
end

for k = 1:7
    set(spfig2(1,k),'xlim',[400 2400])
    set(spfig2(2:4,k),'xlim',[640 850]);
    set(spfig2(1,k),'ylim',[0 0.7]);
    set(spfig2(2,k),'ylim',[0 5]);
    set(spfig2(3,k),'ylim',[0 60]);
    set(spfig2(4,k),'ylim',[0 0.62]);
    
    ylabel(spfig2(1,k),'R')
    ylabel(spfig2(2,k),'L_F ')
    ylabel(spfig2(3,k),'E_F')
    ylabel(spfig2(4,k),'f_{esc}')
    
end
for k = 1:4
    xlabel(spfig2(k,7),'wl (nm)')
end
resizefigure(spfig2,4,7,.1,.12,.07,.04)

%% Figure 2. A plot of TOC fluorescence contributions, for LAI = 3
figure(2), clf

f = dlmread([direc subdirecs{5} 'fluorescence.dat'],'',2,0);
fs = dlmread([direc subdirecs{5} 'fluorescence_sunlit.dat'],'',2,0);
fsc = dlmread([direc subdirecs{5} 'fluorescence_scattered.dat'],'',2,0);
fh = dlmread([direc subdirecs{5} 'fluorescence_shaded.dat'],'',2,0);

z = plot(spectral.wlF,[f(7,:)' fs(7,:)' fh(7,:)' fsc(7,:)'],'k');
xlabel('wl (nm)')
ylabel('L_F (W m^{-2}\mum^{-1}sr^{-1})')
set(z(1),'LineWidth',2)
set(z(3),'LineStyle','--')
set(z(4),'Color','r')
legend('total','sunlit','shaded','scattered')

%% histograms of parameters (not presented in the paper)
figure(12), clf
I1 = find(strcmp(measured.variety,'WT'));
I2 = find(strcmp(measured.variety,'MG'));
X1 = (0:10:90);
X2 = (1:.2:3);
X3 = (0:.005:.06);
X4 = (0:.05:1.2);
vars = {'C_{ab} (\mug cm^{-2})','N','C_{dm} (mg cm^{-2})','C_{s}' };

spfig7 = zeros(4,2);
for k = 1:2
    switch k
        case 1, I = I1;
        case 2, I = I2;
    end
    
    spfig7(1,k) = subplot(2,4,(k-1)*4+1);
    N = hist(leafbio_all.Cab(I),X1);
    bar(X1,N,'k')
    set(gca,'xlim',[0 90])
    
    spfig7(2,k) = subplot(2,4,(k-1)*4+2);
    N = hist(leafbio_all.N(I),X2);
    bar(X2,N,'k')
    set(gca,'xlim',[1 3.1])
    
    spfig7(3,k) = subplot(2,4,(k-1)*4+3);
    N = hist(leafbio_all.Cdm(I),X3);
    bar(X3,N,'k')
    set(gca,'xlim',[-0.003 .06])
    
    spfig7(4,k) = subplot(2,4,(k-1)*4+4);
    N = hist(leafbio_all.Cs(I),X4);
    bar(X4,N,'k')
    set(gca,'xlim',[-.05 0.6])   
end
for k = 1:2
    ylabel(spfig7(1,k),'Number of leaves')
end
for k = 1:4
    xlabel(spfig7(k,2),vars{k})
end
resizefigure(spfig7,4,2,0.07,0.17,0.04,0.08,0.97,0.97)
set(spfig7(:),'FontSize',13)